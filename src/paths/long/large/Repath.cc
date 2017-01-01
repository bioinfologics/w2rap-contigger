///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <util/OutputLog.h>
#include "CoreTools.h"
#include "ParallelVecUtilities.h"
#include "VecUtilities.h"
#include "paths/HyperBasevector.h"
#include "paths/HyperKmerPath.h"
#include "paths/KmerBaseBroker.h"
#include "paths/long/LongReadsToPaths.h"
#include "paths/long/LongProtoTools.h"
#include "paths/long/ReadPath.h"
#include "system/SortInPlace.h"

void RepathInMemory( const HyperBasevector& hb, const vecbasevector& edges,
             const vec<int>& inv, ReadPathVec& paths, const int K, const int K2,
             HyperBasevector& hb2, ReadPathVec& paths2 , const Bool REPATH_TRANSLATE, bool INVERT_PATHS,
             const Bool EXTEND_PATHS )
{
     // Build and unique sort places.  These are the paths, with the following
     // modifications:
     // (a) paths implying base sequence < K2 bases are discarded;
     // (b) if the inverse of a path is smaller we use it instead;
     // (c) places are unique sorted.

     OutputLog(2) << "beginning repathing "<<edges.size()<<" edges from K="<<K<<" to K2="<<K2<< std::endl;
     OutputLog(2) << "constructing places from "<<paths.size()<<" paths" << std::endl;
     uint64_t pathed=0,multipathed=0;
     for (auto &p:paths) {
          if (p.size()>0 ) pathed++;
          if (p.size()>2 ) multipathed++;
     }
     OutputLog(2) <<pathed<<" / "<<paths.size()<<" reads pathed, "<< multipathed << " spanning junctions"<< std::endl;
     std::vector< std::vector<int> > places;
     places.reserve( paths.size( ) );
     //TODO: same path -> same place, why don't we sort and unique paths then? before complicating all of this code!
     const int batch = 10000;
     #pragma omp parallel for
     for (int64_t m = 0; m < (int64_t) paths.size(); m += batch) {
          std::vector<std::vector<int> > placesm;
          placesm.reserve(batch);
          std::vector<int> x, y;
          int64_t n = Min(m + batch, (int64_t) paths.size());
          for (int64_t i = m; i < n; i++) {
               x.clear(), y.clear();
               for (int64_t j = 0; j < (int64_t) paths[i].size(); j++)
                    x.push_back(paths[i][j]);
               int nkmers = 0;
               for (int j = 0; j < x.size(); j++)
                    nkmers += edges[x[j]].size() - ((int) K - 1);
               if (nkmers + ((int) K - 1) < K2) continue;
               for (int j = x.size() - 1; j >= 0; j--)
                    y.push_back(inv[x[j]]);
               placesm.push_back(x < y ? x : y);
          }
          #pragma omp critical
          { places.insert(places.end(),placesm.begin(),placesm.end()); }
     }
     OutputLog(4) << "sorting "<<places.size()<<" path-places" << std::endl;
     sortInPlaceParallel(places.begin(), places.end());
     places.erase(std::unique(places.begin(),places.end()),places.end());
     places.shrink_to_fit();
     OutputLog(3) << places.size()<<" unique places" << std::endl;
     // Add extended places.

     if (EXTEND_PATHS)
     {    OutputLog(2) << "begin extending paths" << std::endl;
          vec<int> to_left, to_right;
          hb.ToLeft(to_left), hb.ToRight(to_right);
          std::vector<std::vector<int>> eplaces;
          for ( auto i = 0; i < places.size( ); i++ )
          {    vec<int> p = places[i];
               int v = to_left[ p.front( ) ], w = to_right[ p.back( ) ];
               while( hb.To(v).solo( ) )
               {    int e = hb.EdgeObjectIndexByIndexTo( v, 0 );
                    if ( !Member( p, e ) ) p.push_front(e);
                    else break;    }
               while( hb.From(w).solo( ) )
               {    int e = hb.EdgeObjectIndexByIndexFrom( w, 0 );
                    if ( !Member( p, e ) ) p.push_back(e);
                    else break;    }
               if ( p.size( ) > places[i].size( ) ) eplaces.push_back(p);    }
          places.insert(places.end(),eplaces.begin(),eplaces.end());
          OutputLog(2) << "resorting" << std::endl;
          sortInPlaceParallel(places.begin(),places.end());
          places.erase(std::unique(places.begin(),places.end()),places.end());
          places.shrink_to_fit();
          OutputLog(2) << "done extending paths" << std::endl;    }

     // Convert places to bases.  For paths of length > 1, we truncate at the
     // beginning and end so that they each contribute at most K2 bases.

     OutputLog(2) << "building all" << std::endl;
     vecbasevector all( places.size( ) );
     vec<int> left_trunc( places.size( ), 0 ), right_trunc( places.size( ), 0 );
#pragma omp parallel for
     for ( int64_t i = 0; i < (int64_t) places.size( ); i++ )
     {    vec<int> e;
          for ( int j = 0; j < (int) places[i].size( ); j++ )
               e.push_back( places[i][j] );
          basevector b = edges[ e[0] ];
          for ( int l = 1; l < e.isize( ); l++ )
          {    b.resize( b.isize( ) - ( K - 1 ) );
               b = Cat( b, edges[ e[l] ] );    }
          if ( e.size( ) > 1 )
          {    int x = e.back( );
               if ( edges[x].isize( ) > K2 )
               {    b.resize( b.isize( ) - ( edges[x].size( ) - K2 ) );
                    right_trunc[i] = edges[x].isize( ) - K2;    }
               x = e.front( );
               if ( edges[x].isize( ) > K2 )
               {    b.SetToSubOf( b, edges[x].size( ) - K2,
                                  b.isize( ) - ( edges[x].size( ) - K2 ) );
                    left_trunc[i] = edges[x].isize( ) - K2;    }    }
          all[i] = b;    }

     // Build HyperBasevector.

     //HyperBasevector hb2;
     vecKmerPath xpaths;
     HyperKmerPath h2;
     OutputLog(2) << "building new graph from places" << std::endl;
     unsigned const COVERAGE = 2u;
     LongReadsToPaths( all, K2, COVERAGE, &hb2, &h2, &xpaths );
     Destroy(all);

     // Write files.

     vec<int> inv2;
     hb2.Involution(inv2);

     // Translate paths to the K=200 graph.  Translation method is very ugly.
     if (REPATH_TRANSLATE)
     {    OutputLog(2) << "translating paths" << std::endl;
          vecKmerPath hpaths;
          vec<big_tagged_rpint> hpathsdb;
          for ( int64_t e = 0; e < h2.EdgeObjectCount( ); e++ )
               hpaths.push_back_reserve( h2.EdgeObject(e) );
          CreateDatabase( hpaths, hpathsdb );
          vec<int> sources, sinks, to_left, to_right;
          h2.Sources(sources), h2.Sinks(sinks);
          h2.ToLeft(to_left), h2.ToRight(to_right);
          vec< vec<int> > ipaths2( xpaths.size( ) );
          vec<int> starts( xpaths.size( ) ), stops( xpaths.size( ) );
          for ( int64_t id = 0; id < (int64_t) xpaths.size( ); id++ )
          {    vec<int> u;
               const KmerPath& p = xpaths[id];
               vec< triple<ho_interval,int,ho_interval> > M, M2;
               int rpos = 0;
               for ( int j = 0; j < p.NSegments( ); j++ )
               {    const KmerPathInterval& I = p.Segment(j);
                    vec<longlong> locs;
                    Contains( hpathsdb, I, locs );
                    for ( int l = 0; l < locs.isize( ); l++ )
                    {    const big_tagged_rpint& t = hpathsdb[ locs[l] ];
                         int hid = t.PathId( );
                         if ( hid < 0 ) continue;
                         longlong hpos = I.Start( ) - t.Start( );
                         longlong start = Max( I.Start( ), t.Start( ) );
                         longlong stop = Min( I.Stop( ), t.Stop( ) );
                         longlong hstart = start - t.Start( );
                         for ( int r = 0; r < t.PathPos( ); r++ )
                              hstart += hpaths[hid].Segment(r).Length( );
                         longlong hstop = hstart + stop - start;
                         longlong rstart = rpos + start - I.Start( );
                         longlong rstop = rstart + stop - start;
                         M.push( ho_interval( rstart, rstop ), hid,
                                 ho_interval( hstart, hstop ) );    }
                    rpos += I.Length( );    }
               Bool bad = False;
               for ( int i = 0; i < M.isize( ); i++ )
               {    int j;
                    for ( j = i + 1; j < M.isize( ); j++ )
                    {    if ( M[j].first.Start( ) != M[j-1].first.Stop( ) + 1 )
                              break;
                         if ( M[j].second != M[j-1].second ) break;
                         if ( M[j].third.Start( ) != M[j-1].third.Stop( ) + 1 )
                              break;    }
                    u.push_back( M[i].second );
                    Bool incomplete = False;
                    if ( i > 0 && M[i].third.Start( ) > 0 ) incomplete = True;
                    if ( j < M.isize( ) && M[j-1].third.Stop( )
                                           != hpaths[ M[i].second ].KmerCount( ) - 1 )
                    {    incomplete = True;
                         bad = True;    }
                    if ( i == 0 && j == M.isize( ) && !incomplete )
                    {    i = j - 1;
                         continue;    }
                    int last = ( i == 0 ? -1 : M2.back( ).first.Stop( ) );
                    if ( M[i].first.Start( ) > last + 1 ) bad = True;
                    M2.push( ho_interval( M[i].first.Start( ),
                                          M[j-1].first.Stop( ) ), M[i].second,
                             ho_interval( M[i].third.Start( ), M[j-1].third.Stop( ) ) );
                    if ( j == M.isize( ) && M[j-1].first.Stop( ) < p.KmerCount( ) - 1 )
                         bad = True;
                    i = j - 1;    }
               if ( !bad && u.nonempty( ) )
               {    starts[id] = M.front( ).third.Start( );
                    stops[id] = hb2.EdgeObject( u.back( ) ).isize( )
                                - ( M.back( ).third.Stop( ) + K2 );     }
               if ( !bad && u.nonempty( ) ) ipaths2[id] = u;    }
          // Parallelizing this loop does not speed it up.  Perhaps to speed it up
          // we have to do something smarter, so as to eliminate the binary search
          // inside the loop.
          for ( int64_t id = 0; id < (int64_t) paths.size( ); id++ )
          {    if ( paths[id].empty( ) ) continue;

               // Note that we have more info here: paths[id].getOffset( )
               // is the start position of the read on the original path.

               vec<int> x, y;
               for ( int64_t j = 0; j < (int64_t) paths[id].size( ); j++ )
                    x.push_back( paths[id][j] );
               int nkmers = 0;
               for ( int j = 0; j < x.isize( ); j++ )
                    nkmers += edges[x[j]].isize( ) - ( (int) K - 1 );
               if ( nkmers + ( (int) K - 1 ) < K2 ) continue;
               for ( int j = x.isize( ) - 1; j >= 0; j-- )
                    y.push_back( inv[ x[j] ] );
               Bool rc = ( y < x );
               x = Min( x, y );
               long pos = BinPosition( places, x );
               long n = ipaths2[pos].size( );

               paths2[id].resize(n);

               int offset;
               if ( !rc )
                    offset = paths[id].getOffset( ) + starts[pos] - left_trunc[pos];
               else offset = paths[id].getOffset( ) + stops[pos] - right_trunc[pos];
               paths2[id].setOffset(offset);

               if ( !rc )
               {    for ( int j = 0; j < n; j++ )
                         paths2[id][j] = ipaths2[pos][j];    }
               else
               {    for ( int j = 0; j < n; j++ )
                         paths2[id][j] = inv2[ ipaths2[pos][n-j-1] ];    }    }
          OutputLog(2) << "paths translation done" << std::endl;
     }
}
