///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "Basevector.h"
#include "CoreTools.h"
#include "ParallelVecUtilities.h"
#include "Qualvector.h"
#include "kmers/BigKPather.h"
#include "paths/HyperBasevector.h"
#include "paths/long/LargeKDispatcher.h"
#include "paths/long/MakeKmerStuff.h"
#include "paths/long/ReadPath.h"
#include "paths/long/large/AssembleGaps.h"
#include "paths/long/large/GapToyTools.h"
#include "paths/long/large/Unsat.h"
#include "system/SortInPlace.h"

template<int M> void MakeStartStop( const vecbasevector& bell, 
     const HyperBasevector& hb, const HyperBasevector& shb, const vec<int>& lefts, 
     const vec<int>& rights, vec<int>& starts, vec<int>& stops )
{    vec< triple<kmer<M>,int,int> > kmers_plus;
     MakeKmerLookup3( bell, kmers_plus );
     for ( int64_t i = 0; i < kmers_plus.jsize( ); i++ )
     {    int64_t j;
          for ( j = i + 1; j < kmers_plus.jsize( ); j++ )
               if ( kmers_plus[j].first != kmers_plus[i].first ) break;
          vec<int> es;
          Bool have_l = False, have_r = False;
          for ( int64_t k = i; k < j; k++ )
          {    int id = kmers_plus[k].second;
               if ( id < shb.EdgeObjectCount( ) ) es.push_back(id);
               else
               {    id -= shb.EdgeObjectCount( );
                    if ( id < lefts.isize( ) ) have_l = True;
                    else have_r = True;    }    }
          if (have_l) starts.append(es);
          if (have_r) stops.append(es);
          i = j - 1;    }    }

template <int M> struct MakeStartStopFunctor
{   void operator( )( const vecbasevector& bell, const HyperBasevector& hb, 
          const HyperBasevector& shb, const vec<int>& lefts, 
          const vec<int>& rights, vec<int>& starts, vec<int>& stops )
    { MakeStartStop<M>( bell, hb, shb, lefts, rights, starts, stops ); }
};

void AssembleGaps2( HyperBasevector& hb, vec<int>& inv2, ReadPathVec& paths2, 
     VecULongVec& paths2_index, vecbasevector& bases, VecPQVec const& quals,
     const String& work_dir, const Bool EXTEND, const Bool ANNOUNCE, 
     const Bool KEEP_ALL_LOCAL, const Bool CONSERVATIVE_KEEP, const Bool INJECT, 
     const Bool LOCAL_LAYOUT, const String DUMP_LOCAL, int K2_FLOOR, 
     const int DUMP_LOCAL_LROOT, const int DUMP_LOCAL_RROOT, 
     vecbvec& new_stuff, const Bool CYCLIC_SAVE,
     const int A2V, const int GAP_CAP, const int MAX_PROX_LEFT, 
     const int MAX_PROX_RIGHT, const int MAX_BPATHS )
{
     // Find clusters of unsatisfied links.

     vec< vec< std::pair<int,int> > > xs;
     Unsat( hb, inv2, paths2, xs, work_dir, A2V );

     // Condense to lists of lefts and rights.

     vec< std::pair<vec<int>,vec<int>> > LR( xs.size( ) );
     #pragma omp parallel for
     for ( int i = 0; i < xs.isize( ); i++ )
     {    vec<int> lefts, rights;
          for ( int j = 0; j < xs[i].isize( ); j++ )
          {    lefts.push_back( xs[i][j].first );
               rights.push_back( xs[i][j].second );    }
          UniqueSort(lefts), UniqueSort(rights);
          LR[i] = make_pair( lefts, rights );    }
     sortInPlaceParallel(LR.begin(),LR.end());

     // Remove inverted copies.  Should force symmetry first.

     {    vec<Bool> lrd( LR.size( ), False );
          PRINT( LR.size( ) );
          #pragma omp parallel for
          for ( int i = 0; i < LR.isize( ); i++ )
          {    vec<int> lefts, rights;
               for ( int j = 0; j < LR[i].first.isize( ); j++ )
                    rights.push_back( inv2[ LR[i].first[j] ] );
               for ( int j = 0; j < LR[i].second.isize( ); j++ )
                    lefts.push_back( inv2[ LR[i].second[j] ] );
               Sort(lefts), Sort(rights);
               if ( make_pair( lefts, rights ) <= LR[i] ) continue;
               if ( !BinMember( LR, make_pair( lefts, rights ) ) ) continue;
               lrd[i] = True;    }
          EraseIf( LR, lrd );
          PRINT( LR.size( ) );    }

     // Some setup stuff.

     int nedges = hb.EdgeObjectCount( );
     int K = hb.K( );

     // Layout reads.  Expensive, temporary (?).

     vec< vec<int> > layout_pos(nedges);
     vec< vec<int64_t> > layout_id(nedges);
     vec< vec<Bool> > layout_or(nedges);
     LayoutReads( hb, inv2, bases, paths2, layout_pos, layout_id, layout_or );

     // Extend terminal edges.

     if (EXTEND) 
          ExtendTerminalEdges( hb, layout_pos, layout_id, layout_or, bases, quals );

     // Make gap assemblies.

     vec<int> mgc = {2};
     int tmpdir_serial=0;
     int min_gap_count = mgc[0], nobj = hb.EdgeObjectCount( );
     vec<int> to_left, to_right;
     hb.ToLeft(to_left), hb.ToRight(to_right);
     vec< vec<basevector> > extras( LR.size( ) );
     vec<String> mreport( LR.size( ) );
     vec<HyperBasevector> mhbp( LR.size( ) );
     std::cout << Date( ) << ": now processing " << LR.size( ) << " blobs" << std::endl;
     std::cout << Date( ) << ": memory in use = " << MemUsageGBString( )
#ifdef __linux
          << ", peak = " << PeakMemUsageGBString( )
#endif
          << std::endl;
     double clockp1 = WallClockTime( );
     int nblobs = LR.size( ), dots_printed = 0, nprocessed = 0;
     int lrc = LR.size( );
     if ( GAP_CAP >= 0 ) lrc = GAP_CAP;
     #pragma omp parallel for schedule(dynamic, 1)
     for ( int bl = 0; bl < lrc; bl++ )
     {    
          // Get ready.

          double aclock1 = WallClockTime( );
          const vec<int> &lefts = LR[bl].first, &rights = LR[bl].second;
          //std::ostringstream mout;
          //mout << "\nbl = " << bl << ", lefts = " << printSeq(lefts)
               //<< ", rights = " << printSeq(rights) << std::endl;
          int K2_FLOOR_LOCAL = K2_FLOOR;

          /*if (ANNOUNCE)
          {   
               #pragma omp critical
               {    std::cout << "\n" << Date( ) << ": START " << bl << ", lefts = "
                         << printSeq(lefts) << ", rights = " 
                         << printSeq(rights) << std::endl;    }    }*/

          // Find relevant pids.
     
          vec<int64_t> pids;

          /*
          if ( A2V == 1 )
          {    const int max_delta = 120;
               for ( int l = 0; l < lefts.isize( ); l++ )
               for ( int k = layout_pos[ lefts[l] ].isize( ) - 1; k >= 0; k-- )
               {    int pos = layout_pos[ lefts[l] ][k]; 
                    int64_t id = layout_id[ lefts[l] ][k];
                    if ( pos + bases[id].isize( ) 
                         < hb.EdgeLengthBases( lefts[l] ) - max_delta ) 
                    {    break;    }
                    pids.push_back(id/2);    }
               for ( int l = 0; l < rights.isize( ); l++ )
               for ( int k = 0; k < layout_pos[ rights[l] ].isize( ); k++ )
               {    int pos = layout_pos[ rights[l] ][k]; 
                    int64_t id = layout_id[ rights[l] ][k];
                    if ( pos > max_delta ) break;
                    pids.push_back(id/2);    }    }
          */

          // if ( A2V >= 2 )
          {
               // Heuristics.

               const int pair_sample = 200;

               // First find the pairs that bridge from left to right, and mark
               // their endpoints.  Inefficient.

               vec<int64_t> pids1;
               vec<vec<int>> lstarts( lefts.size( ) ), rstarts( rights.size( ) );
               {    vec< quad<int64_t,Bool,int,int> > marks;
                    for ( int l = 0; l < lefts.isize( ); l++ )
                    for ( int k = 0; k < layout_pos[ lefts[l] ].isize( ); k++ )
                    {    if ( !layout_or[ lefts[l] ][k] ) continue;
                         int pos = layout_pos[ lefts[l] ][k]; 
                         int64_t id = layout_id[ lefts[l] ][k];
                         marks.push( id/2, False, pos, l );    }
                    for ( int l = 0; l < rights.isize( ); l++ )
                    for ( int k = 0; k < layout_pos[ rights[l] ].isize( ); k++ )
                    {    if ( layout_or[ rights[l] ][k] ) continue;
                         int pos = layout_pos[ rights[l] ][k]; 
                         int64_t id = layout_id[ rights[l] ][k];
                         marks.push( id/2, True, pos, l );    }
                    Sort(marks);
                    for ( int l = 0; l < marks.isize( ); l++ )
                    {    int m;
                         for ( m = l + 1; m < marks.isize( ); m++ )
                              if ( marks[m].first != marks[l].first ) break;
                         Bool left = False, right = False;
                         for ( int j = l; j < m; j++ )
                         {    if ( !marks[j].second ) left = True;
                              else right = True;    }
                         if ( left && right )
                         {    pids1.push_back( marks[l].first );
                              for ( int k = l; k < m; k++ )
                              {    if ( !marks[k].second )
                                   {    lstarts[ marks[k].fourth ].push_back( 
                                             marks[k].third );    }
                                   else
                                   {    rstarts[ marks[k].fourth ].push_back( 
                                             marks[k].third );    }    }    }
                         l = m - 1;    }    }
               UniqueSort(pids1);
               for ( int l = 0; l < lefts.isize( ); l++ )
                    Sort( lstarts[l] );
               for ( int l = 0; l < rights.isize( ); l++ )
                    Sort( rstarts[l] );

               // Now find the pairs that start close to one of the bridge pairs.

               vec<int64_t> pids2;
               for ( int l = 0; l < lefts.isize( ); l++ )
               for ( int k = 0; k < layout_pos[ lefts[l] ].isize( ); k++ )
               {    int pos = layout_pos[ lefts[l] ][k]; 
                    int64_t id = layout_id[ lefts[l] ][k];
                    Bool fw = layout_or[ lefts[l] ][k];
                    if ( BinMember( pids1, id/2 ) ) continue;
                    /*
                    if ( lstarts[l].empty( ) ) // XXXXXXXXXXXXXXXXXXXXXXXXXXXXX
                    {    std::cout << "No lstarts." << std::endl; // XXXXXXXXXXXXXXXXXXXX
                         PRINT( lefts[l] ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
                         Scram(0);    } // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
                    */
                    int low = lstarts[l].front( ), high = lstarts[l].back( );
                    Bool close = False;
                    if ( low <= pos && pos <= high ) close = True;
                    else
                    {    if (fw)
                         {    if ( low > pos && low - pos <= MAX_PROX_LEFT ) 
                              {    close = True;    }
                              else if ( high < pos && pos - high <= MAX_PROX_RIGHT ) 
                              {    close = True;    }    }
                         else
                         {    if ( low > pos && low - pos <= MAX_PROX_RIGHT ) 
                              {    close = True;    }
                              else if ( high < pos && pos - high <= MAX_PROX_LEFT ) 
                              {    close = True;    }    }    }
                    if (close) pids2.push_back(id/2);    }
               for ( int l = 0; l < rights.isize( ); l++ )
               for ( int k = 0; k < layout_pos[ rights[l] ].isize( ); k++ )
               {    int pos = layout_pos[ rights[l] ][k]; 
                    int64_t id = layout_id[ rights[l] ][k];
                    Bool fw = layout_or[ rights[l] ][k];
                    if ( BinMember( pids1, id/2 ) ) continue;
                    /*
                    if ( rstarts[l].empty( ) ) // XXXXXXXXXXXXXXXXXXXXXXXXXXXXX
                    {    std::cout << "No rstarts." << std::endl; // XXXXXXXXXXXXXXXXXXXX
                         PRINT( rights[l] ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
                         Scram(0);    } // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
                    */
                    int low = rstarts[l].front( ), high = rstarts[l].back( );
                    Bool close = False;
                    if ( low <= pos && pos <= high ) close = True;
                    else
                    {    if (fw)
                         {    if ( low > pos && low - pos <= MAX_PROX_LEFT ) 
                              {    close = True;    }
                              else if ( high < pos && pos - high <= MAX_PROX_RIGHT ) 
                              {    close = True;    }    }
                         else
                         {    if ( low > pos && low - pos <= MAX_PROX_RIGHT ) 
                              {    close = True;    }
                              else if ( high < pos && pos - high <= MAX_PROX_LEFT ) 
                              {    close = True;    }    }    }
                    if (close) pids2.push_back(id/2);    }
               UniqueSort(pids2);

               // Now subsample if needed.

               int keep = pair_sample/2;
               if ( pids1.isize( ) + pids2.isize( ) <= pair_sample )
                    pids.append(pids1);
               else if ( pids1.isize( ) <= keep ) pids1.append(pids1);
               else
               {    for ( int l = 0; l < keep; l++ )
                    {    int m = ( l * pids1.isize( ) ) / keep;
                         pids.push_back( pids1[m] );    }    }
               if ( pids.isize( ) + pids2.isize( ) <= pair_sample )
                    pids.append(pids2);
               else
               {    keep = pair_sample - pids.isize( );
                    for ( int l = 0; l < keep; l++ )
                    {    int m = ( l * pids2.isize( ) ) / keep;
                         pids.push_back( pids2[m] );    }    }    }

          UniqueSort(pids);
          //mout << "using " << pids.size( ) << " pairs" << std::endl;
          //mout << "pids = " << printSeq(pids) << std::endl;

          if (ANNOUNCE)
          {   
               #pragma omp critical
               {    std::cout << "\n" << Date( ) << ": MID1 " << bl
                         << ", time used so far = " << TimeSince(aclock1)
                         << ", npids = " << pids.size( ) << std::endl;    }    }

          // Attempt local assembly.

          String TMP;
          if ( KEEP_ALL_LOCAL ) {
               #pragma omp atomic
               tmpdir_serial++;
               TMP = work_dir + "/local/" + ToString( tmpdir_serial );
                     } else
               TMP = work_dir + "/local/" + ToString( omp_get_thread_num() );

          //mout << "assembling in " << TMP << "\n";
          //mout << "total setup time = " << TimeSince(aclock1) << std::endl;
          VecEFasta corrected;
          vecbasevector creads;
          vec<pairing_info> cpartner;
          vec<int> cid;
          LongProtoTmpDirManager tmp_mgr(TMP);
          SupportedHyperBasevector shb;

          int lroot = lefts[0], rroot = rights[0];

          MakeLocalAssembly1( lroot, rroot, hb, bases, quals, pids, TMP,
               LOCAL_LAYOUT, K2_FLOOR_LOCAL, work_dir, corrected, creads, cpartner, 
               cid, tmp_mgr );

          if (ANNOUNCE)
          {   
               #pragma omp critical
               {    std::cout << "\n" << Date( ) << ": MID2 " << bl
                         << ", time used so far = " << TimeSince(aclock1)
                         << std::endl;    }    }

          retry:
          MakeLocalAssembly2( corrected, hb, lefts, rights, shb, INJECT,
               K2_FLOOR_LOCAL, creads, tmp_mgr, cid, cpartner );

          if (ANNOUNCE)
          {   
               #pragma omp critical
               {    std::cout << "\n" << Date( ) << ": MID3 " << bl << ", K = "
                         << shb.K( ) << ", time used so far = " 
                         << TimeSince(aclock1) << std::endl;    }    }

          if ( shb.K( ) == 0 )
          {    // TODO: no more dots in advances...
               /*mreport[bl] = mout.str( );
               Dot( nblobs, nprocessed, dots_printed, ANNOUNCE, bl );*/
               continue;    }

          if ( DUMP_LOCAL != "" )
          {    static Bool special(False);
               if ( !special && lroot == DUMP_LOCAL_LROOT
                         && rroot == DUMP_LOCAL_RROOT )
               {    
                    #pragma omp critical
                    {    std::cout << "dumping special1" << std::endl;
                         Ofstream( out, work_dir + "/special/1.dot" );
                         shb.PrintSummaryDOT0w( out, True, False, True );
                         Ofstream( fout, work_dir + "/special/1.fasta" );
                         for ( int e = 0; e < shb.EdgeObjectCount( ); e++ )
                              shb.EdgeObject(e).Print( fout, e );
                         special = True;    }    }    }

          // Find edges "starts" and "stops" overlapping root edges.

          double sclock = WallClockTime( );
          vec<int> starts, stops;
          vecbasevector bell;
          for ( int e = 0; e < shb.EdgeObjectCount( ); e++ )
               bell.push_back( shb.EdgeObject(e) );
          for ( int l = 0; l < lefts.isize( ); l++ )
               bell.push_back( hb.EdgeObject( lefts[l] ) );
          for ( int r = 0; r < rights.isize( ); r++ )
               bell.push_back( hb.EdgeObject( rights[r] ) );
          BigK::dispatch<MakeStartStopFunctor>( 
               shb.K( ), bell, hb, shb, lefts, rights, starts, stops );
          UniqueSort(starts), UniqueSort(stops);

          if (ANNOUNCE)
          {   
               #pragma omp critical
               {    std::cout << "\n" << Date( ) << ": MID4 " << bl
                         << ", time used so far = " << TimeSince(aclock1)
                         << std::endl;    }    }

          // Reduce shb to those edges between starts and stops.

          vec<int> yto_left, yto_right;
          shb.ToLeft(yto_left), shb.ToRight(yto_right);
          vec<int> keep = Intersection( starts, stops );
          if ( !CONSERVATIVE_KEEP )
          {    keep.append(starts);
               keep.append(stops);    }
          for ( int j1 = 0; j1 < starts.isize( ); j1++ )
          for ( int j2 = 0; j2 < stops.isize( ); j2++ )
          {    int v = yto_right[ starts[j1] ], w = yto_left[ stops[j2] ];
               vec<int> b = shb.EdgesSomewhereBetween( v, w );
               keep.append(b);
               if ( CONSERVATIVE_KEEP && ( b.nonempty( ) || v == w ) )
               {     keep.push_back( starts[j1], stops[j2] );    }    }
          UniqueSort(keep);
          vec<int> ydels;
          for ( int e = 0; e < shb.EdgeObjectCount( ); e++ )
               if ( !BinMember( keep, e ) ) ydels.push_back(e);
          HyperBasevector xshb(shb);
          xshb.DeleteEdges(ydels);
          xshb.RemoveUnneededVertices( );
          xshb.RemoveDeadEdgeObjects( );
          //mout << TimeSince(sclock) << " used contracting" << std::endl;

          if ( DUMP_LOCAL != "" )
          {    static Bool special(False);
               if ( !special && lroot == DUMP_LOCAL_LROOT
                         && rroot == DUMP_LOCAL_RROOT )
               {    
                    #pragma omp critical
                    {    std::cout << "dumping special2" << std::endl;
                         Ofstream( out, work_dir + "/special/2.dot" );
                         xshb.PrintSummaryDOT0w( out, True, False, True );
                         Ofstream( fout, work_dir + "/special/2.fasta" );
                         for ( int e = 0; e < xshb.EdgeObjectCount( ); e++ )
                              xshb.EdgeObject(e).Print( fout, e );
                         special = True;    }    }    }

          // Attempt to recover assemblies with cycles by raising K.

          if ( CYCLIC_SAVE && !xshb.Acyclic( ) )
          {    vec<int> K2s = {100,128,144,172,200};
               int j;
               for ( j = 0; j < K2s.isize( ); j++ )
                    if ( xshb.K( ) < K2s[j] ) break;
               if ( j < K2s.isize( ) )
               {    K2_FLOOR_LOCAL = K2s[j];
                    //PRINT_TO( mout, K2_FLOOR_LOCAL );
                    goto retry;    }    }

          if (ANNOUNCE)
          {   
               #pragma omp critical
               {    std::cout << "\n" << Date( ) << ": MID5 " << bl
                         << ", time used so far = " << TimeSince(aclock1)
                         << std::endl;    }    }

          if ( !xshb.Acyclic( ) || xshb.N( ) == 0 )
          {    //if ( !xshb.Acyclic( ) ) mout << "has cycle, not using" << std::endl;
               //if ( xshb.N( ) == 0 ) mout << "local assembly empty" << std::endl;
               //mreport[bl] += mout.str( );
               //Dot( nblobs, nprocessed, dots_printed, ANNOUNCE, bl );
               continue;    }
          //mout << "local assembly has " << xshb.NComponents( )
          //     << " components" << "\n";

          // Make bpaths.  These are all source-sink paths through the
          // local graph.

          double aclock2 = WallClockTime( );
          vec<basevector> bpaths;
          vec<int> sources, sinks;
          xshb.Sources(sources), xshb.Sinks(sinks);
          vec<int> zto_left, zto_right;
          xshb.ToLeft(zto_left), xshb.ToRight(zto_right);
          for ( int i1 = 0; i1 < sources.isize( ); i1++ )
          for ( int i2 = 0; i2 < sinks.isize( ); i2++ )
          {    vec<vec<int>> p;
               xshb.EdgePaths( zto_left, zto_right, sources[i1], sinks[i2], p );
               for ( int l = 0; l < p.isize( ); l++ )
               {    basevector b = xshb.EdgeObject( p[l][0] );
                    for ( int m = 1; m < p[l].isize( ); m++ )
                    {    b.resize( b.isize( ) - ( xshb.K( ) - 1 ) );
                         b = Cat( b, xshb.EdgeObject( p[l][m] ) );    }
                    bpaths.push_back(b);    
                    if ( bpaths.isize( ) > MAX_BPATHS ) break;    }    }
          //PRINT_TO( mout, bpaths.size( ) );
          if ( bpaths.isize( ) > MAX_BPATHS )
          {    //mout << "Too many bpaths." << std::endl;
               //mreport[bl] += mout.str( );
               //Dot( nblobs, nprocessed, dots_printed, ANNOUNCE, bl );
               continue;    }

          // Make more bpaths.  

          for ( int l = 0; l < lefts.isize( ); l++ )
          {    Bool ext = False;
               for ( int m = 0; m < lefts.isize( ); m++ )
               {    if ( to_right[ lefts[m] ] == to_left[ lefts[l] ] )
                    {    basevector b = hb.EdgeObject( lefts[m] );
                         b.resize( b.isize( ) - (K-1) );
                         b = Cat( b, hb.EdgeObject( lefts[l] ) );
                         bpaths.push_back(b);
                         ext = True;    }    }
               if ( !ext ) bpaths.push_back( hb.EdgeObject( lefts[l] ) );    }
          for ( int r = 0; r < rights.isize( ); r++ )
          {    Bool ext = False;
               for ( int m = 0; m < rights.isize( ); m++ )
               {    if ( to_left[ rights[m] ] == to_right[ rights[r] ] )
                    {    basevector b = hb.EdgeObject( rights[r] );
                         b.resize( b.size( ) - (K-1) );
                         b = Cat( b, hb.EdgeObject( rights[m] ) );
                         bpaths.push_back(b);
                         ext = True;    }    }
               if ( !ext ) bpaths.push_back( hb.EdgeObject( rights[r] ) );    }
          //PRINT_TO( mout, bpaths.size( ) );

          if ( lroot == DUMP_LOCAL_LROOT && rroot == DUMP_LOCAL_RROOT )
               DumpBPaths(bpaths,lroot,rroot,work_dir+"/loc/bpaths");

          // Make the bpaths into a HyperBasevector.

          //mout << "initial patch creation time = " << TimeSince(aclock2) << std::endl;
          vecbasevector bpathsx;
          for ( int l = 0; l < bpaths.isize( ); l++ )
               bpathsx.push_back( bpaths[l] );
          BasesToGraph( bpathsx, K, mhbp[bl] );
          //mout << "patch creation time = " << TimeSince(aclock2) << std::endl;

          if (ANNOUNCE)
          {   
               #pragma omp critical
               {    std::cout << "\n" << Date( ) << ": SAVING " << bl << ", lefts = "
                         << printSeq(lefts) << ", rights = " 
                         << printSeq(rights) << ", npids = " << pids.size( )
                         << ", time used so far = " << TimeSince(aclock1)
                         << std::endl;    }    }

          // Save.
     
          //mreport[bl] += mout.str( );
          //Dot( nblobs, nprocessed, dots_printed, ANNOUNCE, bl );
            }

     std::cout << TimeSince(clockp1) << " spent in local assemblies, "
          << "memory in use = " << MemUsageGBString( )
#ifdef __linux
          << ", peak = " << PeakMemUsageGBString( )
#endif
          << std::endl;
     // Do the patching.

     const vec< std::pair<int,int> > blobs( LR.size( ) );
     Patch( hb, blobs, mhbp, work_dir, mreport, new_stuff );    }
