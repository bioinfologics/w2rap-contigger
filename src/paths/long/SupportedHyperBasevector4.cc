// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "CoreTools.h"

#include "paths/long/MakeKmerStuff.h"
#include "paths/long/SupportedHyperBasevector.h"


namespace { // open anonymous namespace

void FixPath( vec<int>& p, const vec< triple<int,int,int> >& merges,
     const vec< vec<int> >& merges_index, const HyperBasevector& hb_orig, 
     int& left_add, int& right_add )
{    left_add = 0, right_add = 0;
     while(1)
     {    Bool changed = False;
          vec<int> mids;
          for ( int l = 0; l < p.isize( ); l++ )
               if ( p[l] >= 0 ) mids.append( merges_index[ p[l] ] );
          UniqueSort(mids);
          for ( int mj = 0; mj < mids.isize( ); mj++ )
          {    int j = mids[mj];
               int e1 = merges[j].first, e2 = merges[j].second;
               ForceAssert( e1 != e2 );
               int enew = merges[j].third;
               for ( int l = 0; l < p.isize( ); l++ )
               {    if ( l < p.isize( ) - 1 && p[l] == e1 && p[l+1] == e2 )
                    {    p[l] = enew;
                         for ( int m = l+2; m < p.isize( ); m++ )
                              p[m-1] = p[m];
                         p.pop_back( );
                         changed = True;    }
                    else if ( l == p.isize( ) - 1 && p[l] == e1 )
                    {    right_add += hb_orig.EdgeLengthKmers(e2);
                         p[l] = enew;
                         changed = True;    }
                    else if ( l == 0 && p[l] == e2 )
                    {    left_add += hb_orig.EdgeLengthKmers(e1);
                         p[l] = enew;
                         changed = True;    }
                    if (changed) break;    }
               if (changed) break;    }
          if ( !changed ) break;    }    }

} // close anonymous namespace

void SupportedHyperBasevector::RemoveUnneededVertices0( 
     vec< triple<int,int,int> >& merges )
{    vec<int> to_left, to_right;
     ToLeft(to_left), ToRight(to_right);

     for ( int i = 0; i < N( ); i++ )
     {    if ( From(i).size( ) == 1 && To(i).size( ) == 1 && From(i)[0] != i )
          {    int e1 = EdgeObjectIndexByIndexTo( i, 0 );
               int e2 = EdgeObjectIndexByIndexFrom( i, 0 );
               basevector p = Cat( e1, e2 );
               int enew = EdgeObjectCount( );
               int re1 = Inv(e1), re2 = Inv(e2);
               ForceAssert( ( re1 < 0 && re2 < 0 ) || ( re1 >= 0 && re2 >= 0 ) );
               int v = To(i)[0], w = From(i)[0];
               Bool loop = ( v == w && From(v).solo( ) && To(v).solo( ) );
               // v --e1--> i --e2--> w
               merges.push( e1, e2, enew );
               JoinEdges( i, p );
               to_left.push_back(v), to_right.push_back(w);
               if ( re1 < 0 ) InvMutable( ).push_back(-1);
               else if ( re1 == e2 && re2 == e1 )
               {    // * --e1=re2--> * --e2=re1--> *
                    InvMutable( ).push_back(enew);    }
               else if ( re2 == e2 && re1 != e1 )
               {    // * --e1--> * --e2=re2--> * --re1--> *
                    int enew2 = EdgeObjectCount( );
                    merges.push( enew, re1, enew2 );
                    basevector p2 = TrimCat( K( ), p, EdgeObject(re1) );
                    to_left.push_back(v), to_right.push_back( to_right[re1] );
                    JoinEdges( w, p2 );
                    InvMutable( ).push_back( -1, enew2 );    }
               else if ( re1 == e1 && re2 != e2 )
               {    // * --re2--> * --e1=re1--> * --e2--> *
                    int enew2 = EdgeObjectCount( );
                    merges.push( re2, enew, enew2 );
                    basevector p2 = TrimCat( K( ), EdgeObject(re2), p );
                    to_left.push_back( to_left[re2] ), to_right.push_back(w);
                    JoinEdges( v, p2 );
                    InvMutable( ).push_back( -1, enew2 );    }
               else if ( re1 == e1 && re2 == e2 )
               {    if (loop) InvMutable( ).push_back(-1);
                    else
                    {    // not sure if this can happen
                         ForceAssert( 0 == 1 );    }    }
               else
               {    // e1, e2, re1, re2 all different
                    int renew = EdgeObjectCount( );
                    basevector rp = Cat( re2, re1 );
                    merges.push( re2, re1, renew );
                    int ri = to_right[re2];
                    JoinEdges( ri, rp );
                    int rv = to_left[re2], rw = to_right[re1];
                    to_left.push_back(rv), to_right.push_back(rw);
                    InvMutable( ).push_back(renew, enew);    }    }    }    }

void SupportedHyperBasevector::RemoveUnneededVertices( )
{    vec< triple<int,int,int> > merges;
     RemoveUnneededVertices0(merges);
     vec< vec<int> > merges_index( EdgeObjectCount( ) );
     for ( int i = 0; i < merges.isize( ); i++ )
     {    merges_index[ merges[i].first ].push_back(i);
          merges_index[ merges[i].second ].push_back(i);    }
     #pragma omp parallel for
     for ( int i = 0; i < NPaths( ); i++ )
     {    int left_add, right_add;
          FixPath( PathMutable(i), merges, merges_index, *this,
               left_add, right_add );    }
     #pragma omp parallel for
     for ( int i = 0; i < NPairs( ); i++ )
     {    int left_add1, right_add1, left_add2, right_add2;
          FixPath( PairLeftMutable(i), merges, merges_index, *this,
               left_add1, right_add1 );
          FixPath( PairRightMutable(i), merges, merges_index, *this,
               left_add2, right_add2 );
          AddTrim( i, right_add1 + left_add2 );    }
     UniqueOrderPaths( );
     RemoveEdgelessVertices( );    }

void SupportedHyperBasevector::DeleteUnusedPaths( )
{    vec<Bool> used, to_delete( NPaths( ), False );
     Used(used);
     for ( int i = 0; i < NPaths( ); i++ )
     {    for ( int j = 0; j < Path(i).isize( ); j++ )
               if ( Path(i,j) >= 0 && !used[ Path(i,j) ] ) to_delete[i] = True;    }
     EraseIf( PathsMutable( ), to_delete );
     EraseIf( WeightsFwMutable( ), to_delete );    
     EraseIf( WeightsRcMutable( ), to_delete );    
     // The following if is temporary - until origins fully implemented.
     if ( to_delete.size( ) == WeightsFwOrigin( ).size( ) )
     {    EraseIf( WeightsFwOriginMutable( ), to_delete );    
          EraseIf( WeightsRcOriginMutable( ), to_delete );    }
     to_delete.resize_and_set( NPairs( ), False );
     for ( int i = 0; i < NPairs( ); i++ )
     {    vec<int> &p1 = PairLeftMutable(i), &p2 = PairRightMutable(i);
          for ( int pass = 1; pass <= 2; pass++ )
          {    vec<int>& p = ( pass == 1 ? p1 : p2 );
               for ( int j = 0; j < p.isize( ); j++ )
                    if ( p[j] >= 0 && !used[ p[j] ] ) to_delete[i] = True;    }    }
     EraseIf( PairsMutable( ), to_delete );
     EraseIf( PairDataMutable( ), to_delete );    }

void SupportedHyperBasevector::RemoveDeadEdgeObjects( )
{    vec<Bool> used;
     Used(used);
     vec<int> to_new_id( used.size( ), -1 );
     {    int count = 0;
          for ( int i = 0; i < used.isize( ); i++ )
               if ( used[i] ) to_new_id[i] = count++;    }
     vec<int> inv2;
     for ( int i = 0; i < EdgeObjectCount( ); i++ )
     {    if ( !used[i] ) continue;
          if ( !InvDef(i) ) inv2.push_back(-1);
          else inv2.push_back( to_new_id[ Inv(i) ] );    }
     InvMutable( ) = inv2;

     vec<Bool> to_delete( NPaths( ), False );
     for ( int i = 0; i < NPaths( ); i++ )
     {    vec<int>& p = PathMutable(i);
          for ( int j = 0; j < Path(i).isize( ); j++ )
          {    if ( Path(i,j) >= 0 )
               {    int n = to_new_id[ Path(i,j) ];
                    if ( n < 0 ) to_delete[i] = True;
                    else PathMutable(i)[j] = n;    }    }    }
     EraseIf( PathsMutable( ), to_delete );
     // Temporary if.
     if ( WeightsFwOrigin( ).size( ) == WeightsFw( ).size( ) )
     {    EraseIf( WeightsFwOriginMutable( ), to_delete );
          EraseIf( WeightsRcOriginMutable( ), to_delete );    }
     EraseIf( WeightsFwMutable( ), to_delete );
     EraseIf( WeightsRcMutable( ), to_delete );

     to_delete.resize_and_set( NPairs( ), False );
     for ( int i = 0; i < NPairs( ); i++ )
     {    vec<int> &p1 = PairLeftMutable(i), &p2 = PairRightMutable(i);
          for ( int pass = 1; pass <= 2; pass++ )
          {    vec<int>& p = ( pass == 1 ? p1 : p2 );
               for ( int j = 0; j < p.isize( ); j++ )
               {    if ( p[j] >= 0 )
                    {    int n = to_new_id[ p[j] ];
                         if ( n < 0 ) to_delete[i] = True;
                         else p[j] = n;    }    }    }    }
     EraseIf( PairsMutable( ), to_delete );
     EraseIf( PairDataMutable( ), to_delete );

     HyperBasevector::RemoveDeadEdgeObjects( );    }

void SupportedHyperBasevector::UniqueOrderPaths( )
{    
     // Temporary if.
     if ( WeightsFwOrigin( ).size( ) == WeightsFw( ).size( ) )
     {    SortSync( PathsMutable( ), WeightsFwMutable( ), WeightsRcMutable( ),
               WeightsFwOriginMutable( ), WeightsRcOriginMutable( ) );    }
     else SortSync( PathsMutable( ), WeightsFwMutable( ), WeightsRcMutable( ) );
     vec<Bool> to_delete( NPaths( ), False );
     for ( int i = 0; i < NPaths( ); i++ )
     {    int j = Paths( ).NextDiff(i);
          fix64_6 cfw = 0.0, crc = 0.0;
          for ( int k = i; k < j; k++ )
               cfw += WeightFw(k);
          WeightFwMutable(i) = cfw;
          for ( int k = i; k < j; k++ )
               crc += WeightRc(k);
          WeightRcMutable(i) = crc;
          // Temporary if.
          if ( WeightsFwOrigin( ).size( ) == WeightsFw( ).size( ) )
          {    for ( int k = i+1; k < j; k++ )
               {    WeightFwOriginMutable(i).append(
                         WeightFwOriginMutable(k) );
                    WeightRcOriginMutable(i).append(
                         WeightRcOriginMutable(k) );    }
               Sort( WeightFwOriginMutable(i) );
               Sort( WeightRcOriginMutable(i) );    }
          for ( int k = i+1; k < j; k++ )
               to_delete[k] = True;
          if ( cfw + crc == 0 ) to_delete[i] = True;
          i = j - 1;   }
     EraseIf( PathsMutable( ), to_delete );
     // Temporary if.
     if ( WeightsFwOrigin( ).size( ) == WeightsFw( ).size( ) )
     {    EraseIf( WeightsFwOriginMutable( ), to_delete );
          EraseIf( WeightsRcOriginMutable( ), to_delete );    }
     EraseIf( WeightsFwMutable( ), to_delete );
     EraseIf( WeightsRcMutable( ), to_delete );

     // Now do pairs.

     SortSync( PairsMutable( ), PairDataMutable( ) );
     to_delete.resize_and_set( NPairs( ), False );
     for ( int i = 0; i < NPairs( ); i++ )
     {    vec<int> &p1 = PairLeftMutable(i), &p2 = PairRightMutable(i);
          int j = Pairs( ).NextDiff(i);
          vec<pair_point> x;
          for ( int k = i; k < j; k++ )
               x.append( PairData(k) );
          Sort(x);
          PairDataMutable(i) = x;
          for ( int k = i+1; k < j; k++ )
               to_delete[k] = True;
          if ( PairData(i).empty( ) ) to_delete[i] = True;
          i = j - 1;   }
     EraseIf( PairsMutable( ), to_delete );
     EraseIf( PairDataMutable( ), to_delete );    }








