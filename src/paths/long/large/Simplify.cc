///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <paths/PathFinder.h>
#include "CoreTools.h"
#include "Qualvector.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "paths/long/large/GapToyTools.h"
#include "paths/long/large/ImprovePath.h"
#include "paths/long/large/PullAparter.h"
#include "paths/long/large/Simplify.h"
#include "paths/long/large/tools/NhoodInfoCore.h"

void Trace( const String& TRACE_SEQ, const HyperBasevector& hb,
     const String& fin_dir, const int tid )
{    vec<int> hits;
     #pragma omp parallel for
     for ( int e = 0; e < (int) hb.E( ); e++ )
     {    String s = hb.EdgeObject(e).ToString( );
          if ( s.Contains(TRACE_SEQ) )
          {
               #pragma omp critical
               {    hits.push_back(e);    }    }    }
     Sort(hits);
     Mkdir777( fin_dir + "/trace" );
     HyperBasevectorX hbx(hb);
     BinaryWriter::writeFile( fin_dir + "/trace/a.hbx", hbx );
     String command = "DIR=" + fin_dir + "/trace OX=" + fin_dir + "/trace/"
          + ToString(tid) + " NH=True DEPTH=5 SEEDS=";
     for ( int j = 0; j < hits.isize( ); j++ )
     {    if ( j > 0 ) command += ",";
          command += ToString( hits[j] );    }
     NhoodInfoCore(command);    }

void Simplify( const String& fin_dir, HyperBasevector& hb, vec<int>& inv, 
     ReadPathVec& paths, const vecbasevector& bases, const VecPQVec& quals,
     const int MAX_SUPP_DEL, const Bool TAMP_EARLY, const int MIN_RATIO2, 
     const int MAX_DEL2, const Bool PLACE_PARTNERS, 
     const Bool ANALYZE_BRANCHES_VERBOSE2, const String& TRACE_SEQ, 
     const Bool DEGLOOP, const Bool EXT_FINAL, const int EXT_FINAL_MODE, 
     const Bool PULL_APART_VERBOSE, const vec<int>& PULL_APART_TRACE,
     const int DEGLOOP_MODE, const double DEGLOOP_MIN_DIST,
     const Bool IMPROVE_PATHS, const Bool IMPROVE_PATHS_LARGE, 
     const Bool FINAL_TINY, const Bool UNWIND3 )
{
     // Improve read placements and delete funky pairs.

     TestInvolution( hb, inv );
     if (PLACE_PARTNERS) PlacePartners( hb, inv, paths, bases, quals );
     ReroutePaths( hb, inv, paths, bases, quals );
     DeleteFunkyPathPairs( hb, inv, bases, paths, False );

     // Remove unsupported edges in certain situations.

     if ( TRACE_SEQ != "" ) Trace( TRACE_SEQ, hb, fin_dir, 1 );
     {    const int min_mult = 10;
          vec<int> dels;
          {
          vec<int> support( hb.EdgeObjectCount( ), 0 );
          for ( int64_t id = 0; id < (int64_t) paths.size( ); id++ )
          {    for ( int64_t j = 0; j < (int64_t) paths[id].size( ); j++ )
               {    int e = paths[id][j];
                    if ( j >= 1 ) support[e]++;
                    if ( inv[e] >= 0 && j < (int64_t) paths[id].size( ) - 1 )
                         support[ inv[e] ]++;    }    }
          #pragma omp parallel for
          for ( int v = 0; v < hb.N( ); v++ )
          {    if ( hb.From(v).size( ) == 2 )
               {    int e1 = hb.EdgeObjectIndexByIndexFrom( v, 0 );
                    int e2 = hb.EdgeObjectIndexByIndexFrom( v, 1 );
                    if ( support[e1] > support[e2] ) std::swap( e1, e2 );
                    int s1 = support[e1], s2 = support[e2];
                    if ( s1 <= MAX_SUPP_DEL && s2 >= min_mult * Max( 1, s1 ) )
                    {
                         #pragma omp critical
                         {    dels.push_back(e1);    }    }    }    }
          }
          {
          vec<int> support( hb.EdgeObjectCount( ), 0 );
          for ( int64_t id = 0; id < (int64_t) paths.size( ); id++ )
          {    for ( int64_t j = 0; j < (int64_t) paths[id].size( ); j++ )
               {    int e = paths[id][j];
                    if ( j < (int64_t) paths[id].size( ) - 1 ) support[e]++;
                    if ( inv[e] >= 0 && j >= 1 ) support[ inv[e] ]++;    }    }
          #pragma omp parallel for
          for ( int v = 0; v < hb.N( ); v++ )
          {    if ( hb.To(v).size( ) == 2 )
               {    int e1 = hb.EdgeObjectIndexByIndexTo( v, 0 );
                    int e2 = hb.EdgeObjectIndexByIndexTo( v, 1 );
                    if ( support[e1] > support[e2] ) std::swap( e1, e2 );
                    int s1 = support[e1], s2 = support[e2];
                    if ( s1 <= MAX_SUPP_DEL && s2 >= min_mult * Max( 1, s1 ) )
                    {
                         #pragma omp critical
                         {    dels.push_back(e1);    }    }    }    }
          }
          hb.DeleteEdges(dels);
          Cleanup( hb, inv, paths );    }
     if ( TRACE_SEQ != "" ) Trace( TRACE_SEQ, hb, fin_dir, 2 );

     // Clean up assembly.

     RemoveSmallComponents3(hb);
     Cleanup( hb, inv, paths );
     TestInvolution( hb, inv );
     Validate( hb, inv, paths );
     if ( TRACE_SEQ != "" ) Trace( TRACE_SEQ, hb, fin_dir, 3 );
     if (TAMP_EARLY) 
     {    Tamp( hb, inv, paths, 0 );
          TestInvolution( hb, inv );    }
     RemoveHangs( hb, inv, paths, 100 );
     Cleanup( hb, inv, paths );
     TestInvolution( hb, inv );
     Validate( hb, inv, paths );
     vec<int> to_right;
     hb.ToRight(to_right);
     AnalyzeBranches( hb, to_right, inv, paths, True, MIN_RATIO2,
          ANALYZE_BRANCHES_VERBOSE2 );
     Cleanup( hb, inv, paths );
     TestInvolution( hb, inv );
     RemoveHangs( hb, inv, paths, MAX_DEL2 );
     Cleanup( hb, inv, paths );
     TestInvolution( hb, inv );
     RemoveSmallComponents3(hb);
     Cleanup( hb, inv, paths );
     if ( TRACE_SEQ != "" ) Trace( TRACE_SEQ, hb, fin_dir, 4 );
     TestInvolution( hb, inv );
     PopBubbles( hb, inv, bases, quals, paths );
     Cleanup( hb, inv, paths );
     TestInvolution( hb, inv );
     DeleteFunkyPathPairs( hb, inv, bases, paths, False );
     Tamp( hb, inv, paths, 10 );
     TestInvolution( hb, inv );
     RemoveHangs( hb, inv, paths, 700 );
     Cleanup( hb, inv, paths );
     RemoveSmallComponents3(hb);
     Cleanup( hb, inv, paths );
     if ( TRACE_SEQ != "" ) Trace( TRACE_SEQ, hb, fin_dir, 5 );

     // Pull apart.

     {    std::cout << Date() << ": making paths index for pull apart" << std::endl;
          VecULongVec invPaths;
          invert( paths, invPaths, hb.EdgeObjectCount( ) );
          std::cout << Date() << ": pulling apart repeats" << std::endl;
          PullAparter pa(hb,inv,paths,invPaths,PULL_APART_TRACE,
                  PULL_APART_VERBOSE,5,5.0, true, true, true, true);
          size_t count = pa.SeparateAll();
          std::cout << Date() << ": there were " << count << " repeats pulled apart."
               << std::endl;
          std::cout << Date() << ": there were " << pa.getRemovedReadPaths() <<
                    " read paths removed during separation." << std::endl;
          Validate( hb, inv, paths );    }
     {
          std::cout << Date() << ": making paths index for PathFinder" << std::endl;
          VecULongVec invPaths;
          invert( paths, invPaths, hb.EdgeObjectCount( ) );
          std::cout << Date() << ": PathFinder: untangling simple choices" << std::endl;
          PathFinder(hb,inv,paths,invPaths).untangle_single_choices();
          std::cout<<"refreshing all structures as precaution"<<std::endl;
          inv.clear();
          hb.Involution(inv);
          std::cout<<"all structures refreshed"<<std::endl;
     }
     // Improve paths.

     if (IMPROVE_PATHS) 
     {    path_improver pimp;
          vec<int64_t> ids;
          ImprovePaths( paths, hb, inv, bases, quals, ids, pimp,
               IMPROVE_PATHS_LARGE, False );    }

     // Extend paths.

     if (EXT_FINAL)
     {    vec<int> to_left;
          hb.ToLeft(to_left), hb.ToRight(to_right);
          int ext = 0;
          auto qvItr = quals.begin();
          for ( int64_t id = 0; id < (int64_t) paths.size( ); id++,++qvItr )
          {    Bool verbose = False;
               const int min_gain = 20;
               ReadPath p = paths[id];
               ExtendPath2( paths[id], id, hb, to_left, to_right, bases[id], *qvItr,
                    min_gain, verbose, EXT_FINAL_MODE );    
               if ( p != paths[id] ) ext++;    }    
          std::cout << ext << " paths extended" << std::endl;    }

     // Degloop.

     if (DEGLOOP) 
     {    Degloop( DEGLOOP_MODE, hb, inv, paths, bases, quals, DEGLOOP_MIN_DIST );
          RemoveHangs( hb, inv, paths, 700 );
          Cleanup( hb, inv, paths );    }

     // Unwind three-edge plasmids.

     if (UNWIND3) UnwindThreeEdgePlasmids( hb, inv, paths );

     // Remove tiny stuff.

     if (FINAL_TINY)
     {    std::cout << Date( ) << ": removing small components" << std::endl;
          RemoveSmallComponents3( hb, True );
          Cleanup( hb, inv, paths );    
          CleanupLoops( hb, inv, paths );
          RemoveUnneededVerticesGeneralizedLoops( hb, inv, paths );    }    }
