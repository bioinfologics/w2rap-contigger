// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <util/OutputLog.h>
#include "CoreTools.h"
#include "Qualvector.h"
#include "graph/FindCells.h"
#include "kmers/KmerRecord.h"
#include "paths/HyperBasevector.h"
#include "paths/long/KmerCount.h"
#include "paths/long/ReadPath.h"
#include "paths/long/large/GapToyTools.h"



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

     OutputLog(2) << TimeSince(clock) << " creating patches" << std::endl;
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
                         const QualVec& q = qvItr[id];

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
