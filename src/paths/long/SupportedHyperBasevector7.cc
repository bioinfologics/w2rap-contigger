// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "Basevector.h"
#include "paths/HyperBasevector.h"
#include "paths/long/LargeKDispatcher.h"
#include "paths/long/MakeKmerStuff.h"
#include "paths/long/SupportedHyperBasevector.h"

// TransformPaths: replace each occurrence of x[i] by y[i].

void SupportedHyperBasevector::TransformPaths( const vec< vec<int> >& x,
     const vec<int>& y, vec< vec< std::pair<int,int> > >& paths_index )
{    if ( paths_index.nonempty( ) ) paths_index.resize( EdgeObjectCount( ) );
     vec<int> all;
     for ( int j = 0; j < x.isize( ); j++ )
          all.append( x[j] );
     UniqueSort(all);

     vec<int> used_paths;
     if ( paths_index.nonempty( ) )
     {    for ( int i = 0; i < all.isize( ); i++ )
          {    int e = all[i];
               for ( int j = 0; j < paths_index[e].isize( ); j++ )
                    used_paths.push_back( paths_index[e][j].first );    }
          UniqueSort(used_paths);    }
     int to_use = ( paths_index.empty( ) ? NPaths( ) : used_paths.isize( ) );

     const int infty = 1000000000;
     #pragma omp parallel for
     for ( int pi = 0; pi < to_use; pi++ )
     {    int i = ( paths_index.empty( ) ? pi : used_paths[pi] );
          vec<int>& p = PathMutable(i);
          vec<int> p_old = p;
          for ( int j = 0; j < p.isize( ); j++ )
          {    if ( !BinMember( all, p[j] ) ) continue;

               // Find all the proper overlaps between p and some x[l], that use
               // the jth entry of p.

               vec< std::pair<int,int> > pos;
               for ( int l = 0; l < x.isize( ); l++ )
               {    for ( int m = 0; m < x[l].isize( ); m++ )
                    {    if ( x[l][m] != p[j] ) continue;
                         if ( Overlap( p, x[l], j - m ) ) pos.push(l,m);    }    }
               if ( !pos.solo( ) )
               {    p.clear( );
                    WeightFwMutable(i) = 0.0;
                    WeightRcMutable(i) = 0.0;
                    break;    }
               int l = pos[0].first, m = pos[0].second;
               int start = j;
               int stop = Min( p.isize( ), x[l].isize( ) + (j-m) );
               for ( int r = start; r < stop - 1; r++ )
                    p[r] = -1;
               p[stop-1] = y[l];    }
          RemoveNegatives(p);

          if ( paths_index.nonempty( ) && p != p_old )
          {    for ( int l = 0; l < p_old.isize( ); l++ )
               {    int e = p_old[l];
                    #pragma omp critical
                    {    auto low = LowerBound( paths_index[e], {i,0} );
                         auto high = UpperBound( paths_index[e], {i,infty} );
                         paths_index[e].erase( 
                              paths_index[e].begin( ) + low,
                              paths_index[e].begin( ) + high );    }    }
               #pragma omp critical
               {    for ( int l = 0; l < p.isize( ); l++ )
                    {    int e = p[l];
                         paths_index[e].push( i, l );
                         Sort( paths_index[e] );    }    }    }    }

     #pragma omp parallel for
     for ( int i = 0; i < NPairs( ); i++ )
     {    vec<int> &p1 = PairLeftMutable(i), &p2 = PairRightMutable(i);
          for ( int pass = 1; pass <= 2; pass++ )
          {    vec<int>& p = ( pass == 1 ? p1 : p2 );
               for ( int j = 0; j < p.isize( ); j++ )
               {    if ( !BinMember( all, p[j] ) ) continue;
                    vec< std::pair<int,int> > pos;
                    for ( int l = 0; l < x.isize( ); l++ )
                    {    for ( int m = 0; m < x[l].isize( ); m++ )
                         {    if ( x[l][m] != p[j] ) continue;
                              if ( Overlap( p, x[l], j - m ) ) 
                                   pos.push(l,m);    }    }
                    if ( !pos.solo( ) )
                    {    p1.clear( ), p2.clear( );
                         PairDataMutable(i).clear( );
                         break;    }
                    int l = pos[0].first, m = pos[0].second;
                    int start = j, stop = Min( p.isize( ), x[l].isize( ) + (j-m) );
                    int add = 0;
                    if ( pass == 1 && j - m + x[l].isize( ) > p.isize( ) )
                    {    for ( int r = p.isize( ); r < x[l].isize( ) + j-m; r++ )
                              add += EdgeLengthKmers( x[l][ r - (j-m) ] );    }
                    if ( pass == 2 && j - m < 0 )
                    {    for ( int r = 0; r < -(j-m); r++ )
                              add += EdgeLengthKmers( x[l][r] );    }
                    AddTrim( i, add );
                    for ( int r = start; r < stop - 1; r++ )
                         p[r] = -1;
                    p[stop-1] = y[l];    }
               RemoveNegatives(p);    }    }    }


void SupportedHyperBasevector::writeBinary( BinaryWriter& writer ) const
{    HyperBasevector::writeBinary(writer);
     writer.write( Inv( ) );
     writer.write( Paths( ) );
     writer.write( WeightsFw( ) );
     writer.write( WeightsRc( ) );
     writer.write( WeightsFwOrigin( ) );
     writer.write( WeightsRcOrigin( ) );
     writer.write( Pairs( ) );
     writer.write( PairData( ) );
     writer.write( ReadCount( ) );
     writer.write( ReadLengthDist( ) );
     writer.write( FudgeMult( ) );    }

void SupportedHyperBasevector::readBinary( BinaryReader& reader )
{    HyperBasevector::readBinary(reader);
     reader.read( &InvMutable( ) );
     reader.read( &PathsMutable( ) );
     reader.read( &WeightsFwMutable( ) );
     reader.read( &WeightsRcMutable( ) );
     reader.read( &WeightsFwOriginMutable( ) );
     reader.read( &WeightsRcOriginMutable( ) );
     reader.read( &PairsMutable( ) );
     reader.read( &PairDataMutable( ) );
     reader.read( &ReadCountMutable( ) );
     reader.read( &ReadLengthDistMutable( ) );
     reader.read( &FudgeMultMutable( ) );    }



