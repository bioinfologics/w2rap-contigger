///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "CoreTools.h"
//#include "ParallelVecUtilities.h"
#include "Qualvector.h"
#include "graph/FindCells.h"
#include "kmers/KmerRecord.h"
#include "kmers/MakeLookup.h"
#include "math/Functions.h"
#include "paths/HyperBasevector.h"
#include "paths/RemodelGapTools.h"
#include "paths/long/KmerCount.h"
#include "paths/long/ReadPath.h"
#include "paths/long/large/GapToyTools.h"

// AnalyzeBranches: note not adjusting to_right.  This is wrong.

void AnalyzeBranches( HyperBasevector& hb, vec<int>& to_right, const vec<int>& inv2, 
     ReadPathVec& paths2, const Bool ANALYZE_BRANCHES_REV, 
     const int min_ratio2, const Bool ANALYZE_BRANCHES_VERBOSE )

{    double clock0 = WallClockTime( );
     vec<int> to_left;
     hb.ToLeft(to_left);
     for ( int64_t i = 0; i < (int64_t) paths2.size( ); i++ )
     {    ReadPath& p = paths2[i];
          for ( int64_t j = 0; j < (int64_t) p.size( ); j++ )
          {    if ( p[j] >= hb.EdgeObjectCount( ) ) p[j] = -1;
               if ( j > 0 && p[j-1] >= 0 && p[j] >= 0 
                    && to_right[ p[j-1] ] != to_left[ p[j] ] )
               {    p[j] = -1;    }    }    }

     // Heuristics.

     const int max_dist = 4;
     const int min_ratio = 5;
     const int max_kill = 2;

     vec< std::pair<int,int> > breaks;
     vec< vec<int> > froms( hb.EdgeObjectCount( ) ), tos( hb.EdgeObjectCount( ) );
     LogTime( clock0, "analyzing branches 0" );
     double clock1 = WallClockTime( );
     for ( int pass = 1; pass <= 2; pass++ )
     {    const int batch = 10000;
          int64_t npids = paths2.size( )/2;
          #pragma omp parallel for
          for ( int64_t bi = 0; bi < npids; bi += batch )
          {    vec< std::pair<int,int> > PP;
               for ( int64_t pid = bi; pid < Min( bi + batch, npids ); pid++ )
               {    vec<int> x, y;
                    for ( int64_t j = 0; j < (int64_t) paths2[2*pid].size( ); j++ )
                         x.push_back( paths2[2*pid][j] );
                    for ( int64_t j = 0; j < (int64_t) paths2[2*pid+1].size( ); j++ )
                         y.push_back( paths2[2*pid+1][j] );
                    y.ReverseMe( );
                    for ( int j = 0; j < y.isize( ); j++ )
                         if ( y[j] >= 0 ) y[j] = inv2[ y[j] ];
                    if ( pass == 2 )
                    {    swap( x, y );
                         x.ReverseMe( ), y.ReverseMe( );
                         for ( int j = 0; j < x.isize( ); j++ )
                              if ( x[j] >= 0 ) x[j] = inv2[ x[j] ];
                         for ( int j = 0; j < y.isize( ); j++ )
                              if ( y[j] >= 0 ) y[j] = inv2[ y[j] ];    }
                    std::pair< vec<int>, vec<int> > p = std::make_pair( x, y );
                    vec< std::pair<int,int> > P;
                    for ( int j1 = 0; j1 < p.first.isize( ) - 1; j1++ )
                    {    if ( p.first[j1] >= 0 && p.first[j1+1] >= 0 )
                              P.push( p.first[j1], p.first[j1+1] );    }
                    for ( int j1 = 0; j1 < p.second.isize( ) - 1; j1++ )
                    {    if ( p.second[j1] >= 0 && p.second[j1+1] >= 0 )
                              P.push( p.second[j1], p.second[j1+1] );    }
                    for ( int j1 = 0; j1 < p.first.isize( ); j1++ )
                    {    int x1 = p.first[j1];
                         if ( x1 >= 0 )
                         {    int m = Position( p.second, x1 );
                              if ( m < 0 && p.second.nonempty( ) 
                                   && p.second[0] >= 0 ) 
                              {    P.push( x1, p.second[0] );    }    }    }
                    UniqueSort(P);
                    PP.append(P);    }
               #pragma omp critical
               {    for ( int j = 0; j < PP.isize( ); j++ )
                    {    froms[ PP[j].first ].
                              push_back( PP[j].second );
                         tos[ PP[j].second ].
                              push_back( PP[j].first );    }    }    }    }
     LogTime( clock1, "analyzing branches 1" );
     double clock1b = WallClockTime( );
     #pragma omp parallel for
     for ( int e = 0; e < hb.EdgeObjectCount( ); e++ )
     {    Sort( froms[e] );
          Sort( tos[e] );    }
     if (ANALYZE_BRANCHES_VERBOSE) std::cout << "\nforward reach:\n";
     LogTime( clock1b, "analyzing branches 1b" );
     double clock2 = WallClockTime( );



     for ( int e = 0; e < hb.EdgeObjectCount( ); e++ )
     {    int v = to_right[e];
          if ( hb.From(v).size( ) <= 1 ) continue;
          if ( hb.To(v).size( ) > 1 ) continue;
          vec< vec<int> > follow( hb.From(v).size( ) );
          vec<int> branches;
          for ( int j = 0; j < hb.From(v).isize( ); j++ )
          {    int f = hb.EdgeObjectIndexByIndexFrom( v, j );
               branches.push_back(f);    }
          int nbranches = branches.size( );

          for ( int j = 0; j < hb.From(v).isize( ); j++ )
          {    int f = hb.EdgeObjectIndexByIndexFrom( v, j );
               int w = to_right[f];
               for ( int l = 0; l < hb.From(w).isize( ); l++ )
                    follow[j].push_back( hb.EdgeObjectIndexByIndexFrom(w, l) );    }

          for ( int dpass = 1; dpass < max_dist; dpass++ )
          {    for ( int i = 0; i < nbranches; i++ )
               {    int n = follow[i].size( );
                    for ( int j = 0; j < n; j++ )
                    {    int w = to_right[ follow[i][j] ];
                         follow[i].append( hb.FromEdgeObj(w) );    }
                    UniqueSort( follow[i] );    }    }

          vec<int> fr, count;
          for ( int i = 0; i < froms[e].isize( ); i++ )
          {    int j = froms[e].NextDiff(i);
               int c = j - i, f = froms[e][i];
               fr.push_back(f), count.push_back(c);
               i = j - 1;    }
          vec<Bool> to_delete( fr.size( ), False );
          for ( int i = 0; i < fr.isize( ); i++ )
          {    vec<int> homes;
               for ( int j = 0; j < follow.isize( ); j++ )
                    if ( Member( follow[j], fr[i] ) ) homes.push_back(j);
               if ( homes.size( ) == follow.size( ) ) count[i] = 0;
               if ( homes.solo( ) )
               {    for ( int j = 0; j < fr.isize( ); j++ )
                    {    if ( fr[j] == hb.EdgeObjectIndexByIndexFrom( v, homes[0] ) )
                         {    count[j] += count[i];
                              count[i] = 0;    }    }    }    }
          for ( int i = 0; i < fr.isize( ); i++ )
               if ( count[i] == 0 ) to_delete[i] = True;
          EraseIf( fr, to_delete ), EraseIf( count, to_delete );
          vec<int> s1 = fr, s2 = branches;
          Sort(s1), Sort(s2);
          if ( s1 == s2 && s1.size( ) == 2 )
          {    if ( count[0] < min_ratio * count[1] 
                    && count[1] < min_ratio * count[0] )
               {    continue;    }    }
          ReverseSortSync( count, fr );
          if (ANALYZE_BRANCHES_VERBOSE)
          {    std::cout << e << " -->";
               for ( int i = 0; i < fr.isize( ); i++ )
                    std::cout << " " << fr[i] << "[" << count[i] << "]";    }
          if ( count.size( ) >= 2 && count[0] >= min_ratio2 * Max( 1, count[1] )
               && count[1] <= max_kill && Member( branches, fr[0] ) )
          {    if (ANALYZE_BRANCHES_VERBOSE) std::cout << " -- RECOMMEND PRUNING";
               for ( int j = 0; j < branches.isize( ); j++ )
                    if ( branches[j] != fr[0] ) breaks.push( e, branches[j] );    }
          if (ANALYZE_BRANCHES_VERBOSE) std::cout << "\n";    }
     UniqueSort(breaks);
     for ( int i = 0; i < breaks.isize( ); i++ )
     {    int e = breaks[i].first, f = breaks[i].second;
          int n = hb.N( );
          hb.AddVertices(2);
          hb.GiveEdgeNewFromVx( f, to_right[e], n );
          to_left[f] = n;
          int re = inv2[e], rf = inv2[f];
          if ( re >= 0 && rf >= 0 ) 
          {    hb.GiveEdgeNewToVx( rf, to_right[rf], n+1 );
               to_right[rf] = n+1;    }    }

     if (ANALYZE_BRANCHES_REV)
     {
     vec< std::pair<int,int> > breaksr;
     if (ANALYZE_BRANCHES_VERBOSE) std::cout << "\nbackward reach:\n";

     for ( int e = 0; e < hb.EdgeObjectCount( ); e++ )
     {    int v = to_left[e];
          if ( hb.To(v).size( ) <= 1 ) continue;
          if ( hb.From(v).size( ) > 1 ) continue;
          vec< vec<int> > preceed( hb.To(v).size( ) );
          vec<int> branches;
          for ( int j = 0; j < hb.To(v).isize( ); j++ )
          {    int f = hb.EdgeObjectIndexByIndexTo( v, j );
               branches.push_back(f);    }
          int nbranches = branches.size( );

          for ( int j = 0; j < hb.To(v).isize( ); j++ )
          {    int f = hb.EdgeObjectIndexByIndexTo( v, j );
               int w = to_left[f];
               for ( int l = 0; l < hb.To(w).isize( ); l++ )
                    preceed[j].push_back( hb.EdgeObjectIndexByIndexTo(w, l) );    }

          for ( int dpass = 1; dpass < max_dist; dpass++ )
          {    for ( int i = 0; i < nbranches; i++ )
               {    int n = preceed[i].size( );
                    for ( int j = 0; j < n; j++ )
                    {    int w = to_left[ preceed[i][j] ];
                         preceed[i].append( hb.ToEdgeObj(w) );    }
                    UniqueSort( preceed[i] );    }    }

          vec<int> fr, count;
          for ( int i = 0; i < tos[e].isize( ); i++ )
          {    int j = tos[e].NextDiff(i);
               int c = j - i, f = tos[e][i];
               if ( to_right[f] == to_left[e] )
               {   fr.push_back(f), count.push_back(c);    }
               i = j - 1;    }
          vec<Bool> to_delete( fr.size( ), False );
          for ( int i = 0; i < fr.isize( ); i++ )
          {    vec<int> homes;
               for ( int j = 0; j < preceed.isize( ); j++ )
                    if ( Member( preceed[j], fr[i] ) ) homes.push_back(j);
               if ( homes.size( ) == preceed.size( ) ) count[i] = 0;
               if ( homes.solo( ) )
               {    for ( int j = 0; j < fr.isize( ); j++ )
                    {    if ( fr[j] == hb.EdgeObjectIndexByIndexTo( v, homes[0] ) )
                         {    count[j] += count[i];
                              count[i] = 0;    }    }    }    }
          for ( int i = 0; i < fr.isize( ); i++ )
               if ( count[i] == 0 ) to_delete[i] = True;
          EraseIf( fr, to_delete ), EraseIf( count, to_delete );
          vec<int> s1 = fr, s2 = branches;
          Sort(s1), Sort(s2);
          if ( s1 == s2 && s1.size( ) == 2 )
          {    if ( count[0] < min_ratio * count[1] 
                    && count[1] < min_ratio * count[0] )
               {    continue;    }    }
          ReverseSortSync( count, fr );
          if (ANALYZE_BRANCHES_VERBOSE)
          {    std::cout << e << " <--";
               for ( int i = 0; i < fr.isize( ); i++ )
                    std::cout << " " << fr[i] << "[" << count[i] << "]";    }
          if ( count.size( ) >= 2 && count[0] >= min_ratio2 * Max( 1, count[1] )
               && count[1] <= max_kill && Member( branches, fr[0] ) )
          {    if (ANALYZE_BRANCHES_VERBOSE) std::cout << " -- RECOMMEND PRUNING";
               for ( int j = 0; j < branches.isize( ); j++ )
                    if ( branches[j] != fr[0] ) breaksr.push( branches[j], e );    }
          if (ANALYZE_BRANCHES_VERBOSE) std::cout << "\n";    }
     UniqueSort(breaksr);
     for ( int i = 0; i < breaksr.isize( ); i++ )
     {    int e = breaksr[i].first, f = breaksr[i].second;
          int n = hb.N( );
          hb.AddVertices(2);
          hb.GiveEdgeNewToVx( e, to_left[f], n );
          to_right[e] = n;
          
          int re = inv2[e], rf = inv2[f];
          if ( re >= 0 && rf >= 0 ) 
          {    hb.GiveEdgeNewFromVx( re, to_left[re], n+1 );
               to_left[re] = n+1;    }    }
     breaks.append(breaksr);
     }

     int nb = breaks.size( );
     for ( int i = 0; i < nb; i++ )
          breaks.push( inv2[ breaks[i].second ], inv2[ breaks[i].first ] );
     UniqueSort(breaks);
     #pragma omp parallel for
     for ( int64_t i = 0; i < (int64_t) paths2.size( ); i++ )
     {    ReadPath& p = paths2[i];
          Bool bad = False;
          for ( int j = 0; j < ( (int) p.size( ) ) - 1; j++ )
          {    std::pair<int,int> x = std::make_pair( p[j], p[j+1] );
               if ( BinMember( breaks, x ) ) bad = True;    }
          if (bad) p.resize(0);    }
     if (ANALYZE_BRANCHES_VERBOSE) std::cout << "\n";
     LogTime( clock2, "analyzing branches 2" );    }


void LayoutReads(const HyperBasevector &hb, const vec<int> &inv,
                 const vecbasevector &bases, const ReadPathVec &paths,
                 std::vector<std::vector<int>> &layout_pos, std::vector<std::vector<int64_t>> &layout_id,
                 std::vector<std::vector<bool>> &layout_or) {
     int nedges = hb.EdgeObjectCount();
     layout_pos.resize(nedges), layout_id.resize(nedges), layout_or.resize(nedges);
     for (int64_t i = 0; i < (int64_t) paths.size(); ++i) {
          vec<int> x;
          for (int64_t j = 0; j < (int64_t) paths[i].size(); ++j)
               x.push_back(paths[i][j]);
          if (x.empty()) continue;
          int pos = paths[i].getOffset();
          for (int j = 0; j < x.isize(); j++) {
               if (j > 0 && j < x.isize() - 1) continue;
               layout_pos[x[j]].push_back(pos);
               layout_id[x[j]].push_back(i);
               layout_or[x[j]].push_back(true);
               pos -= hb.EdgeLengthKmers(x[j]);
          }
          x.ReverseMe();
          for (int j = 0; j < x.isize(); j++)
               x[j] = inv[x[j]];
          pos = paths[i].getOffset() + bases[i].isize();
          int len = hb.EdgeLength(x[0]);
          for (int j = 1; j < x.isize(); j++)
               len += hb.EdgeLengthKmers(x[j]);
          pos = len - pos;
          for (int j = 0; j < x.isize(); j++) {
               if (j > 0 && j < x.isize() - 1) continue;
               layout_pos[x[j]].push_back(pos);
               layout_id[x[j]].push_back(i);
               layout_or[x[j]].push_back(false);
               pos -= hb.EdgeLengthKmers(x[j]);
          }
     }
     #pragma omp parallel for
     for (int e = 0; e < nedges; e++)
          SortSync(layout_pos[e], layout_or[e], layout_id[e]);
}



void RemoveHangs( HyperBasevector& hb, vec<int>& inv, ReadPathVec& paths, 
     const int max_del )
{
     const double junk_ratio = 10.0;
     vec<kmer_count> kc( hb.EdgeObjectCount( ) );
     for ( int e = 0; e < hb.EdgeObjectCount( ); e++ )
          kc[e].n = hb.EdgeObject(e).isize( ) - hb.K( ) + 1;
     digraphE<kmer_count> shb_kc( hb, kc );
     RemoveHangingEnds3( shb_kc, &kmer_count::N, max_del, junk_ratio, 100 );//XXX: 100 is an arbitrary parameter.
     vec<int> e_to_delete;
     vec<Bool> used;
     shb_kc.Used(used);
     for ( int e = 0; e < hb.EdgeObjectCount( ); e++ )
          if ( !used[e] ) e_to_delete.push_back(e);    
     hb.DeleteEdges(e_to_delete);
     }

void Patch(HyperBasevector &hb, vec<HyperBasevector> &mhbp, vecbvec &new_stuff) {
     double clock = WallClockTime();
     new_stuff.clear();
     new_stuff.reserve(std::accumulate(mhbp.begin(), mhbp.end(), 0ul,
                                       [](size_t nnn, HyperBasevector const &hbv) {
                                           nnn += hbv.E();
                                           auto end = hbv.To().end();
                                           auto itr2 = hbv.From().begin();
                                           for (auto itr = hbv.To().begin(); itr != end; ++itr, ++itr2)
                                                nnn += itr->size() * itr2->size();
                                           return nnn;
                                       }));
     int K = hb.K();
     for (int bl = 0; bl < mhbp.isize(); bl++) {

          if (mhbp[bl].N() > 0) {
               HyperBasevector const &hbp = mhbp[bl];
               // Insert patch.
               for (int e = 0; e < hbp.EdgeObjectCount(); e++)
                    new_stuff.push_back(hbp.EdgeObject(e));
               for (int v = 0; v < hbp.N(); v++)
                    for (int i1 = 0; i1 < hbp.To(v).isize(); i1++)
                         for (int i2 = 0; i2 < hbp.From(v).isize(); i2++) {
                              basevector const &e1 = hbp.EdgeObjectByIndexTo(v, i1);
                              basevector const &e2 = hbp.EdgeObjectByIndexFrom(v, i2);
                              new_stuff.push_back(TrimCat(K, e1, e2));
                         }
          }
     }

     std::cout << Date() << ": "<< TimeSince(clock) << " used patching" << std::endl;
}

template<class H> void DegloopCore( const int mode, H& hb, vec<int>& inv, 
     ReadPathVec& paths, const vecbasevector& bases, const VecPQVec& quals,
     const VecULongVec& paths_index, const int v, const int pass,
     const double min_dist, vec<int>& EDELS, const int verbosity,
     const vec<int>* ids )
{
     int K = hb.K( );
     int n = ( pass == 1 ? hb.From(v).size( ) : hb.To(v).size( ) );
     if ( n >= 2 )
     {    vec<vec<int>> qs(n);

          // Don't mess with homopolymers.

          Bool homop = False;
          for ( int i = 0; i < n; i++ )
          {    int e = ( pass == 1 ? hb.IFrom( v, i ) : hb.ITo( v, i ) );
               if ( hb.Bases(e) == 0 ) continue;
               int ne = hb.Bases(e);
               vec<char> b;
               const int hcount = 10;
               if ( pass == 1 )
               {    for ( int j = 0; j < hcount; j++ )
                         b.push_back( hb.EdgeObject(e)[ K - j - 1 ] );    }
               else
               {    for ( int j = 0; j < hcount; j++ )
                         b.push_back( hb.EdgeObject(e)[ ne - K + j ] );    }
               UniqueSort(b);
               if ( b.solo( ) ) homop = True;    }
          if (homop) return;

          // Proceed.

          int min_edge = 1000000000;
          for ( int i = 0; i < n; i++ )
          {    int e = ( pass == 1 ? hb.IFrom( v, i ) : hb.ITo( v, i ) );
               if ( hb.Bases(e) == 0 ) continue;
               min_edge = Min( min_edge, hb.Bases(e) );    }
          auto qvItr = quals.begin();
          for ( int i = 0; i < n; i++ )
          {    int e = ( pass == 1 ? hb.IFrom( v, i ) : hb.ITo( v, i ) );
               if ( hb.Bases(e) == 0 ) continue;
               int ne = hb.Bases(e), re = inv[e];
               for ( int xpass = 1; xpass <= 2; xpass++ )
               {    int x = ( xpass == 1 ? e : re );
                    for (int64_t j = 0; j < (int64_t) paths_index[x].size( ); j++)
                    {    int64_t id = paths_index[x][j];
                         const ReadPath& p = paths[id];
                         const basevector& b = bases[id];
                         const qualvector& q = qvItr[id];

                         // Set homopolymer base quality to min across it.

                         /*
                         vec<int> q( bases[id].size( ) );
                         for ( int j = 0; j < b.isize( ); j++ )
                              q[j] = quals[id][j];
                         for ( int j = 0; j < b.isize( ); j++ )
                         {    int k;
                              for ( k = j + 1; k < b.isize( ); k++ )
                                   if ( b[k] != b[j] ) break;
                              int m = 1000000000;
                              for ( int l = j; l < k; l++ )
                                   m = Min( m, (int) q[l] );
                              for ( int l = j; l < k; l++ )
                                   q[l] = m;
                              j = k - 1;    }
                         */

                         for ( int l = 0; l < (int) p.size( ); l++ )
                         {    if ( p[l] != x ) continue;
                              int estart = p.getOffset( ); // start of read on edge
                              for ( int m = 0; m < l; m++ )
                                   estart -= hb.Kmers( p[m] );
                              int estop = estart + b.size( ); // stop of read on edge
                              int rpos = ( xpass == 1 ^ pass == 1 ?
                                   -estart + ne - K : -estart + K - 1 );
                              if ( rpos < 0 || rpos >= b.isize( ) ) continue;

                              // Very stringent condition, probably too stringent
                              // in some cases.

                              if ( !( xpass == 1 ^ pass == 1 ) )
                              {    if ( IntervalOverlap(
                                        0, min_edge, estart, estop ) < K )
                                   {    continue;   }    }
                              else
                              {    if ( IntervalOverlap(
                                        ne - min_edge, ne, estart, estop ) < K )
                                   {    continue;   }    }

                              /*
                              if ( b[rpos] != hb.EdgeObject(x)[
                                   xpass == 1 ^ pass == 1 ? ne-K : K-1 ] )
                              {    continue;    }
                              */

                              if ( verbosity >= 3 )
                              {    ForceAssert( ids != NULL );
                                   std::cout << "read " << (*ids)[id]
                                        << " supports edge " << e << " with quality "
                                        << int(q[rpos]) << std::endl;    }

                              qs[i].push_back( q[rpos] );    }    }    }
               ReverseSort( qs[i] );    }

          // Assess quality score distribution difference.

          vec<double> m( n, -1 );
          vec<int> k(n);
          for ( int i = 0; i < n; i++ )
          {    k[i] = qs[i].size( );
               if ( qs[i].nonempty( ) ) m[i] = Mean( qs[i] );    }
          vec<int> dels;
          vec<double> dists;
          for ( int i1 = 0; i1 < n; i1++ )
          for ( int i2 = 0; i2 < n; i2++ )
          {    if ( i1 == i2 ) continue;
               int good1 = 0;
               for ( int j = 0; j < qs[i1].isize( ); j++ )
                    if ( qs[i1][j] >= 30 ) good1++;
               int good2 = 0;
               for ( int j = 0; j < qs[i2].isize( ); j++ )
                    if ( qs[i2][j] >= 30 ) good2++;
               int e2 = ( pass == 1 ? hb.IFrom( v, i2 ) : hb.ITo( v, i2 ) );
               int ne2 = hb.Kmers(e2);

               if ( mode >= 2 && k[i2] == 0 && good1 >= 10 && ne2 <= 200 )
                    dels.push_back(i2);

               if ( k[i1] == 0 || k[i2] == 0 ) continue;
               double dist = (m[i1]-m[i2]) 
                    / sqrt( m[i1]*m[i1]/k[i1] + m[i2]*m[i2]/k[i2] );
               Bool kill2 = ( dist >= min_dist && good2 <= 1 && ne2 <= 200 );
               if (kill2) dels.push_back(i2);
               if ( kill2 || ( dist >= 0 && verbosity >= 2 ) )
                    dists.push_back(dist);    }
          UniqueSort(dels);

          if ( verbosity >= 2 )
          {    std::cout << "\n" << ( pass == 1 ? "from" : "to" ) << "\n";
               for ( int i = 0; i < n; i++ )
               {    int e = ( pass == 1 ? hb.IFrom( v, i ) : hb.ITo( v, i ) );
                    std::cout << e << ": " << printSeq( qs[i] ) << std::endl;    }
               std::cout << "dist = " << printSeq(dists) << std::endl;    }

          // Record edges for deletion.

          if ( dels.nonempty( ) )
          {
               #pragma omp critical
               {    vec<int> edels;
                    if ( verbosity == 1 )
                         std::cout << "\n" << ( pass == 1 ? "from" : "to" ) << "\n";
                    for ( int i = 0; i < n; i++ )
                    {    int e = ( pass == 1 ? hb.IFrom( v, i ) : hb.ITo( v, i ) );
                         if ( Member( dels, i ) ) edels.push_back(e);
                         if ( verbosity == 1 )
                              std::cout << e << ": " << printSeq( qs[i] ) << std::endl;    }
                    Sort(edels);
                    EDELS.append(edels);
                    if ( verbosity >= 1 ) 
                         std::cout << "delete edges: " << printSeq(edels) << std::endl;
                    if ( verbosity == 1 )
                    {    std::cout << "dist = " << printSeq(dists)
                              << std::endl;    }    }    }    }    }

template void DegloopCore( const int mode, HyperBasevector& hb, vec<int>& inv, 
     ReadPathVec& paths, const vecbasevector& bases, const VecPQVec& quals,
     const VecULongVec& paths_index, const int v, const int pass,
     const double min_dist, vec<int>& EDELS, const int verbosity,
     const vec<int>* ids );

// Go through branch points.
// Score branches by computing quality score at Kth base.
// Uses version of quality score distribution test from DivineBubbles, but
// less sophisticated.

void Degloop( const int mode, HyperBasevector& hb, vec<int>& inv, ReadPathVec& paths,
     const vecbasevector& bases, const VecPQVec& quals, const double min_dist,
     const int verbosity )
{
     VecULongVec paths_index;
     invert( paths, paths_index );

     // Main loop.

     int K = hb.K( );
     vec<int> EDELS;
     #pragma omp parallel for
     for ( int v = 0; v < hb.N( ); v++ )
     {    for ( int pass = 1; pass <= 2; pass++ )
          {    DegloopCore( mode, hb, inv, paths, bases, quals, paths_index, 
                    v, pass, min_dist, EDELS, verbosity );    }    }
     int ed = EDELS.size( );
     for ( int i = 0; i < ed; i++ )
          EDELS.push_back( inv[ EDELS[i] ] );
     UniqueSort(EDELS);
     hb.DeleteEdges(EDELS);
     }
