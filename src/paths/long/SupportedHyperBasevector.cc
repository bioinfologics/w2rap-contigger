///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "CoreTools.h"
#include "ParallelVecUtilities.h"
#include "Set.h"
#include "VecUtilities.h"
#include "fastg/FastgGraph.h"
#include "math/Functions.h"
#include "paths/BigMapTools.h"
//#include "paths/PairedPair.h"
#include "paths/HyperBasevector.h"
#include "paths/long/CleanEfasta.h"
#include "paths/long/DigraphFromWords.h"
#include "paths/long/LargeKDispatcher.h"
#include "paths/long/LongProtoTools.h"
#include "paths/long/MakeKmerStuff.h"
#include "paths/long/PairInfo.h"
#include "paths/long/SupportedHyperBasevector.h"

namespace { // open anonymous namespace

vec<int> GetNextEdges( const HyperBasevector& hb, const vec< vec<int> >& paths,
     const vec<fix64_6>& counts_fw, const vec<fix64_6>& counts_rc, 
     const vec< vec< std::pair<int,int> > >& paths_index,
     const int median_read, const vec<int>& p, const fix64_6 keep_score, 
     const int win_ratio, const int verbosity, std::ostringstream& out )
{    
     // First find the edges that follow p, from an overlap of a path with p.

     vec<int> nexts;
     int b = p.back( );
     ForceAssertGe( b, 0 );
     for ( int k = 0; k < paths_index[b].isize( ); k++ )
     {    int id2 = paths_index[b][k].first, pos2 = paths_index[b][k].second;
          if ( pos2 == paths[id2].isize( ) - 1 ) continue;
          int n = paths[id2][pos2+1];
          if ( n < 0 || Member( nexts, n ) ) continue;
          Bool mismatch = False;
          for ( int p2 = pos2 - 1; p2 >= 0; p2-- )
          {    int p1 = p.isize( ) - 1 - ( pos2 - p2 );
               if ( p1 < 0 ) break;
               if ( p[p1] != paths[id2][p2] )
               {    mismatch = True;
                    break;    }    }
          if ( !mismatch ) nexts.push_back( paths[id2][pos2+1] );    }
     Sort(nexts);
     if ( nexts.size( ) > 1 )
     {    
          // Try to eliminate some nexts using the nonexistence criterion.

          if ( verbosity >= 2 )
          {    out << "- checking " << printSeq( p, " " );
               out << " { " << printSeq( nexts, " " ) << " }\n";    }
          for ( int s = 0; s < p.isize( ); s++ )
          {    vec<fix64_6> score( nexts.size( ), 0 );
               int len = 0;
               for ( int j = s + 1; j < p.isize( ); j++ )
               {    ForceAssertGe( p[j], 0 );
                    len += hb.EdgeLengthKmers( p[j] );    }
               if ( len > median_read ) continue;
               vec<int> p1, offsets;
               p1.SetToSubOf( p, s, p.isize( ) - s );
               for ( int l = 0; l < nexts.isize( ); l++ )
               {    vec<int> p2(p1);
                    p2.push_back( nexts[l] );

                    // Compute the number of reads that contain p2.

                    vec<int> ids;
                    int b = p2.back( );
                    for ( int k = 0; k < paths_index[b].isize( ); k++ )
                    {    int id3 = paths_index[b][k].first; 
                         int pos3 = paths_index[b][k].second;
                         if ( pos3 < p2.isize( ) - 1 ) continue;
                         Bool mismatch = False;
                         for ( int p3 = pos3 - 1; p3 >= 0; p3-- )
                         {    int p2x = p2.isize( ) - 1 - ( pos3 - p3 );
                              if ( p2x < 0 ) break;
                              if ( p2[p2x] != paths[id3][p3] )
                              {    mismatch = True;
                                   break;    }    }
                         if (mismatch) continue;
                         ids.push_back(id3);    }
                    UniqueSort(ids);
                    for ( int i = 0; i < ids.isize( ); i++ )
                         score[l] += counts_fw[ ids[i] ] + counts_rc[ ids[i] ];    }
               if ( verbosity >= 3 )
               {    PRINT2_TO( out, s, len );
                    out << "score = " << printSeq( score, " " ) << "\n";    }
               vec<Bool> to_delete( nexts.size( ), False );
               for ( int l = 0; l < nexts.isize( ); l++ )
               {    if ( score[l] < keep_score 
                         && Max(score) >= win_ratio * score[l] ) 
                    {    to_delete[l] = True;     }    }
               EraseIf( nexts, to_delete );    }    }
     return nexts;    }

vec<int> Betweens( const HyperBasevector& hb, const vec<int>& to_right, 
     const vec<int>& q, const bool debug = false )
{    vec<int> d, v;
     d.push_back( q.front( ), q.back( ) );
     digraphE<basevector> G1(hb);
     G1.DeleteEdges(d);
     for ( int j = 0; j < q.isize( ) - 1; j++ )
          v.push_back( to_right[ q[j] ] );

     if ( debug ) {
	  std::cout << "BETWEENS: v=" << printSeq(v) << std::endl;
     }
     return G1.EdgesConnectedTo(v);    }

vec<int> Betweens2(const HyperBasevector& hb, const vec<int>& to_left, const vec<int>& to_right,
	  const vec<int>& q, const bool debug = false)
{
     vec<int> v;
     typedef std::pair<int, int > VertexPair;
     vec<VertexPair > v_verboten(2);		// edges between any VertexPair here will be not be explored
						// this is how we isolate the box.  Size is hardcoded.

     // all interior vertices on v
     for (int j = 0; j < q.isize() - 1; j++)
	  v.push_back(to_right[q[j]]);

     // connections from the exterior to the interior are on the "verboten" list
     // we've hard-coded the size as the tests below are hard-coded for speed.
     v_verboten[0] = VertexPair( to_left[ q.front() ], to_right[ q.front() ] );
     v_verboten[1] = VertexPair( to_left[ q.back() ], to_right[ q.back() ] );

     // find all vertices that can be reached from vertices in v, without crossing one of the
     // verboten edges.  This is a modified version of VerticesConnectedTo.
     vec<int> v_neighbors( v );
     UniqueSort(v_neighbors);
     while ( 1 ) {
	  size_t n1 = v_neighbors.size();
	  for ( size_t i = 0; i < n1; i++ ) {

	       const vec<int>& from = hb.From( v_neighbors[i] );	// from this vertex to others
	       for ( size_t j = 0; j < from.size(); ++j ) {
		   VertexPair test( v_neighbors[i], from[j] );
		   if ( test != v_verboten[0] && test != v_verboten[1] )
			 v_neighbors.push_back( from[j] );
	       }
	       const vec<int>& to = hb.To( v_neighbors[i] );		// to this vertex from others
	       for ( size_t j = 0; j < to.size(); ++j ) {
		    VertexPair test( to[j], v_neighbors[i] );
		    if ( test != v_verboten[0] && test != v_verboten[1] )
			 v_neighbors.push_back( to[j] );
	       }
	  }
	  UniqueSort(v_neighbors);
	  if (v_neighbors.size() == n1)
	       break;

     }

     // debugging cruft
     if ( debug ) {
	  std::cout << "q=" << printSeq(q) << std::endl;
	  std::cout << "to_left[q]=" << to_left[q[0]] << "," << to_left[q[1]] << std::endl;
	  std::cout << "to_right[q]=" << to_right[q[0]] << "," << to_right[q[1]] << std::endl;
	  std::cout << "v_neighbors=" << printSeq(v_neighbors) << std::endl;
	  for ( size_t i = 0; i < v_neighbors.size(); ++i ) {
	       std::cout << "vneighbors[" << i << "].From()=" << printSeq(hb.From(v_neighbors[i])) << std::endl;
	       std::cout << "vneighbors[" << i << "].To()=" << printSeq(hb.To(v_neighbors[i])) << std::endl;
	  }
	  vec<int> dels = Betweens( hb, to_right, q, true );
	  std::cout << "dels = " << printSeq(dels) << std::endl;
     }

     // Again, a modified version of EdgesConnectedTo: find all edges connected to a list of vertices, excluding
     // edges between verboten pairs.
     vec<int> e;
     for ( size_t i = 0; i < v_neighbors.size(); i++ ) {
	  for ( size_t j = 0; j < hb.From( v_neighbors[i] ).size(); j++ ) {
	       int x = hb.EdgeObjectIndexByIndexFrom( v_neighbors[i], j );
	       VertexPair test(v_neighbors[i], hb.From(v_neighbors[i])[j] );
	       if ( test != v_verboten[0] && test != v_verboten[1]  )
		    e.push_back(x);
	  }
	  for ( size_t j = 0; j < hb.To( v_neighbors[i] ).size(); j++ ) {
	       int x = hb.EdgeObjectIndexByIndexTo( v_neighbors[i], j );
	       VertexPair test( hb.To(v_neighbors[i])[j], v_neighbors[i] );
	       if ( test != v_verboten[0] && test != v_verboten[1] )
		    e.push_back(x);
	  }
     }
     UniqueSort(e);
     return e;
}

class seq_place {

     public:

     seq_place( ) { }
     seq_place( const vec<int>& x, const int g, const Bool fw, const int start,
          const int stop ) : x(x), g(g), fw(fw), start(start), stop(stop) { }

     vec<int> x;    // seq is s[x][0], s[x[1]], ...
     int g;         // placed at G[g]
     Bool fw;       // placement orientation
     int start;     // start of placement
     int stop;      // stop of placement

};

} // close anonymous namespace

Bool Overlap( const vec<int>& v, const vec<int>& w, const int o )
{    for ( int i = 0; i < v.isize( ); i++ )
     {    int j = i - o;
          if ( j < 0 || j >= w.isize( ) ) continue;
          if ( v[i] != w[j] ) return False;    }
     return True;    }

Bool PullApartProcessVertex( SupportedHyperBasevector& shb, vec<int>& to_left, 
     vec<int>& to_right, vec< vec< std::pair<int,int> > >& paths_index, 
     const int v, const int w, const int pass, const double min_weight_split, 
     const long_logging& logc )
{    int id1 = ( pass == 1 ? 0 : 1 ), id2 = ( pass == 1 ? 1 : 0 );
     int x1 = shb.EdgeObjectIndexByIndexTo( v, id1 );
     int x2 = shb.EdgeObjectIndexByIndexTo( v, id2 );
     int r = shb.EdgeObjectIndexByIndexFrom( v, 0 );
     int y1 = shb.EdgeObjectIndexByIndexFrom( w, 0 );
     int y2 = shb.EdgeObjectIndexByIndexFrom( w, 1 );
     if ( x1 == y1 || x1 == y2 || x2 == y1 || x2 == y2 ) return False;

     fix64_6 weight_11 = 0, weight_12 = 0, weight_21 = 0, weight_22 = 0;
     vec<int> p11, p12, p21, p22;
     p11.push_back(x1,r,y1), p12.push_back(x1,r,y2);
     p21.push_back(x2,r,y1), p22.push_back(x2,r,y2);

     vec<int> p11_paths, p12_paths, p21_paths, p22_paths;
     #pragma omp parallel for
     for ( int l = 0; l < paths_index[x1].isize( ); l++ )
     {    int i = paths_index[x1][l].first, j = paths_index[x1][l].second;
          const vec<int>& p = shb.Path(i);
          if ( p.Contains( p11, j ) )
          {
               #pragma omp critical
               {    p11_paths.push_back(i);    }    }
          if ( p.Contains( p12, j ) )
          {
               #pragma omp critical
               {    p12_paths.push_back(i);    }    }    }
     #pragma omp parallel for
     for ( int l = 0; l < paths_index[x2].isize( ); l++ )
     {    int i = paths_index[x2][l].first, j = paths_index[x2][l].second;
          const vec<int>& p = shb.Path(i);
          if ( p.Contains( p21, j ) )
          {
               #pragma omp critical
               {    p21_paths.push_back(i);    }    }
          if ( p.Contains( p22, j ) )
          {
               #pragma omp critical
               {    p22_paths.push_back(i);    }    }    }
     UniqueSort(p11_paths), UniqueSort(p12_paths);
     UniqueSort(p21_paths), UniqueSort(p22_paths);
     for ( int l = 0; l < p11_paths.isize( ); l++ )
     {    int i = p11_paths[l];
          weight_11 += shb.Weight(i);    }
     for ( int l = 0; l < p12_paths.isize( ); l++ )
     {    int i = p12_paths[l];
          weight_12 += shb.Weight(i);    }
     for ( int l = 0; l < p21_paths.isize( ); l++ )
     {    int i = p21_paths[l];
          weight_21 += shb.Weight(i);    }
     for ( int l = 0; l < p22_paths.isize( ); l++ )
     {    int i = p22_paths[l];
          weight_22 += shb.Weight(i);    }

     const int min_weight_split_low = 2;
     Bool OK = False;
     if ( weight_11 >= min_weight_split_low && weight_22 >= min_weight_split_low
          && weight_12 == 0 && weight_21 == 0 )
     {    OK = True;    }
     if ( weight_11 >= min_weight_split && weight_22 >= min_weight_split
          // && weight_12 <= 1 && weight_21 <= 1
          && weight_12 + weight_21 <= 2
          && shb.EdgeLengthKmers(r) <= shb.MedianCorrectedReadLengthFudge( ) )
     {    OK = True;    }
     if ( weight_11 >= min_weight_split/2
          && weight_22 >= min_weight_split/2
          && weight_12 + weight_21 < 2
          && shb.EdgeLengthKmers(r) <= shb.MedianCorrectedReadLengthFudge( ) )
     {    OK = True;    } 
     if ( !OK ) return False;

     // Find images of edges under involution.

     int rx1 = shb.Inv(x1), ry1 = shb.Inv(y1), rr = shb.Inv(r);
     int rx2 = shb.Inv(x2), ry2 = shb.Inv(y2);

     // Test for first special case.
     //
     // a1 --x1-->            --y1=ry2-->               --rx2--> c1
     //            v --r--> w             rv --rr--> rw
     // a2 --x2-->            --y2=ry1-->               --rx1--> c2
     //
     // This pulls apart to x1,r,y1,rr,rx2
     //                     x2,r,y2,rr,rx1,
     // which are rc to each other.

     Bool special1 = ( rx1 >= 0 && y1 == ry2 && r != rr
          && rx1 != x1 && rx1 != x2 && rx2 != x1 && rx2 != x2 );
     if (special1)
     {    if ( logc.PULL_APART_DEBUG ) std::cout << "\n";
          if ( logc.verb[ "PULL_APART" ] >= 1 || logc.PULL_APART_DEBUG )
          {    std::cout << "pulling part (special1) " << x1 << "," << r << "," << y1 
                    << " from " << x2 << "," << r << "," << y2 << std::endl;    }
          if ( logc.PULL_APART_DEBUG ) 
          {    PRINT8( x1, rx1, x2, rx2, r, rr, y1, y2 );
               PRINT4( v, w, id1, id2 );    }
          const basevector &X1 = shb.EdgeObject(x1), &RX1 = shb.EdgeObject(rx1); 
          const basevector &X2 = shb.EdgeObject(x2); 
          const basevector &RX2 = shb.EdgeObject(rx2);
          const basevector &R = shb.EdgeObject(r); 
          const basevector &RR = shb.EdgeObject(rr);
          const basevector &Y1 = shb.EdgeObject(y1), &Y2 = shb.EdgeObject(y2);
          basevector Z1 = shb.Cat( x1, r, y1, rr, rx2 );
          basevector Z2 = shb.Cat( x2, r, y2, rr, rx1 );
          int a1 = shb.To(v)[id1], a2 = shb.To(v)[id2];
          int rv = to_left[rr], rw = to_right[rr];
          shb.DeleteEdgesAtVertex(v), shb.DeleteEdgesAtVertex(w);
          shb.DeleteEdgesAtVertex(rw);
          int c1 = to_right[rx2], c2 = to_right[rx1];
          int z1 = shb.EdgeObjectCount( ), z2 = shb.EdgeObjectCount( ) + 1;
          shb.InvMutable( ).push_back( z2, z1 );
          shb.AddEdge(a1,c1,Z1), shb.AddEdge(a2,c2,Z2);
          to_left.push_back(a1,a2), to_right.push_back(c1,c2);
          vec< vec<int> > e(2);
          e[0].push_back( x1, r, y1, rr, rx2 );
          e[1].push_back( x2, r, y2, rr, rx1 );
          vec<int> f;
          f.push_back( z1, z2 );
          shb.TransformPaths( e, f, paths_index );
          return True;    }

     // Find images of edges under involution (different labeling).

     rx1 = shb.Inv(y1), ry1 = shb.Inv(x1), rr = shb.Inv(r),
     rx2 = shb.Inv(y2), ry2 = shb.Inv(x2);

     // Check for consistency with involution.

     Bool eq = False;
     if ( rx1 >= 0 )
     {    Bool left_eq = ( rx1 == x1 && rx2 == x2 ) || ( rx1 == x2 && rx2 == x1 );
          Bool right_eq = ( ry1 == y1 && ry2 == y2 ) || ( ry1 == y2 && ry2 == y1 );
          eq = left_eq && right_eq && rr == r;
          vec<int> s, t;
          if ( !eq )
          {    s.push_back( x1, x2, r, y1, y2 );
               t.push_back( rx1, rx2, rr, ry1, ry2 );
               UniqueSort(s), UniqueSort(t);    }
          if ( !eq && Meet(s,t) ) return False;    }

     // Two passes.

     vec< vec<int> > e;
     vec<int> f;
     if ( logc.PULL_APART_DEBUG ) 
     {    std::cout << "\n";
          PRINT6( x1, rx1, r, rr, y1, ry1 );
          PRINT4( x2, rx2, y2, ry2 );
          PRINT4( v, w, id1, id2 );
          PRINT4( to_left[x1], to_right[x1], to_left[x2], to_right[x2] );
          PRINT2( to_left[r], to_right[r] );
          PRINT4( to_left[y1], to_right[y1], to_left[y2], to_right[y2] );
          PRINT4( to_left[rx1], to_right[rx1], to_left[rx2], to_right[rx2] );
          PRINT2( to_left[rr], to_right[rr] );
          PRINT4( to_left[ry1], to_right[ry1], to_left[ry2], to_right[ry2] );    }
     int a1_1 = -1, a1_2 = -1, a2_1 = -1, a2_2 = -1;
     int b1_1 = -1, b1_2 = -1, b2_1 = -1, b2_2 = -1;
     basevector Z1_1, Z1_2, Z2_1, Z2_2;
     int v_1 = -1, v_2 = -1, w_1 = -1, w_2 = -1;
     int npasses = ( !eq ? 2 : 1 );
     for ( int xpass = 1; xpass <= npasses; xpass++ )
     {    if ( xpass == 2 && rx1 < 0 ) continue;
          int px1, px2, py1, py2, pr;
          if ( xpass == 1 )
          {    px1 = x1, px2 = x2, py1 = y1, py2 = y2, pr = r;    }
          else
          {    px1 = rx1, px2 = rx2, py1 = ry1, py2 = ry2;
               pr = rr;    }
          int v = to_left[pr], w = to_right[pr];
          if ( logc.PULL_APART_DEBUG ) PRINT3( xpass, v, w );
          int id1 = ( shb.EdgeObjectIndexByIndexTo(v,0) == px1 ? 0 : 1 );
          int id2 = 1 - id1;
          int id1b
               = ( shb.EdgeObjectIndexByIndexFrom(w,0) == py1 ? 0 : 1 );
          if ( xpass == 1 )
          {    a1_1 = shb.To(v)[id1], a2_1 = shb.To(v)[id2];    }
          else
          {    a1_2 = shb.To(v)[id1], a2_2 = shb.To(v)[id2];    }
          if ( xpass == 1 )
          {    b1_1 = shb.From(w)[id1b], b2_1 = shb.From(w)[1-id1b];    }
          else
          {    b1_2 = shb.From(w)[id1b], b2_2 = shb.From(w)[1-id1b];    }
          ( xpass == 1 ? Z1_1 : Z1_2 ) = shb.Cat( px1, pr, py1 );
          ( xpass == 1 ? Z2_1 : Z2_2 ) = shb.Cat( px2, pr, py2 );
          ( xpass == 1 ? v_1 : v_2 ) = v;
          ( xpass == 1 ? w_1 : w_2 ) = w;     }

     // Check for an illegal condition.

     if ( eq && ReverseComplement(Z1_1) != Z2_1 ) 
     {    if ( logc.PULL_APART_DEBUG ) std::cout << "illegal, bailing\n";
          return False;    }

     for ( int xpass = 1; xpass <= npasses; xpass++ )
     {    if ( xpass == 2 && rx1 < 0 ) continue;
          int px1, px2, py1, py2, pr;
          if ( xpass == 1 )
          {    px1 = x1, px2 = x2, py1 = y1, py2 = y2, pr = r;    }
          else
          {    px1 = rx1, px2 = rx2, py1 = ry1, py2 = ry2, pr = rr;    }
          vec<int> p11, p22;
          p11.push_back(px1,pr,py1), p22.push_back(px2,pr,py2);
          int v = ( xpass == 1 ? v_1 : v_2 );
          int w = ( xpass == 1 ? w_1 : w_2 );

          // Announce.

          if ( logc.verb[ "PULL_APART" ] >= 1 || logc.PULL_APART_DEBUG )
          {    std::cout << "pulling part " << px1 << "," << pr << "," << py1 
                    << " from " << px2 << "," << pr << "," << py2 << std::endl;    }

          // Edit.

          basevector Z1 = ( xpass == 1 ? Z1_1 : Z1_2 );
          basevector Z2 = ( xpass == 1 ? Z2_1 : Z2_2 );
          int a1 = ( xpass == 1 ? a1_1 : a1_2 );
          int a2 = ( xpass == 1 ? a2_1 : a2_2 );
          int b1 = ( xpass == 1 ? b1_1 : b1_2 );
          int b2 = ( xpass == 1 ? b2_1 : b2_2 );

          if ( logc.PULL_APART_DEBUG ) 
               std::cout << "deleting edges at " << v << " and " << w << std::endl;
          shb.DeleteEdgesAtVertex(v), shb.DeleteEdgesAtVertex(w);
          int z1 = shb.EdgeObjectCount( ), z2 = shb.EdgeObjectCount( ) + 1;
          if ( logc.PULL_APART_DEBUG ) 
          {    std::cout << "adding edge from " << a1 << " to " << b1 << std::endl;
               std::cout << "adding edge from " << a2 << " to " << b2 << std::endl;
               PRINT2( Z1.size( ), Z2.size( ) );    }
          shb.AddEdge(a1,b1,Z1), shb.AddEdge(a2,b2,Z2);
          if (eq) shb.InvMutable( ).push_back( z2, z1 );
          else if ( rx1 < 0 ) shb.InvMutable( ).push_back( -1, -1 );
          else if ( xpass == 1 ) 
          {    shb.InvMutable( ).push_back( z1 + 2, z2 + 2 );    }
          else if ( xpass == 2 ) 
          {    shb.InvMutable( ).push_back( z1 - 2, z2 - 2 );    }
          to_left.push_back(a1,a2), to_right.push_back(b1,b2);

          e.push_back( p11, p22 );
          f.push_back( z1, z2 );     }
     shb.TransformPaths( e, f, paths_index );
     return True;     }

