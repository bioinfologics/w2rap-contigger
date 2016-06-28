///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Experimental code to close gaps at the K=60 stage.  This works according to
// the following schema:
// (1) Somehow find pairs of edges (e1,e2), the "left root" and "right root" that
//     it seems to make sense to try to close the gap between.
// (2) Somehow associate a list of pids to each such pair.
// (3) Somehow form a local assembly from the associated reads and extract the
//     part of it bridging from e1 to e2.
// (4) Formally combine the the local assemblies with the preexisting assembly and
//     rebuild.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "Basevector.h"
#include "MainTools.h"
#include "PairsManager.h"
#include "Qualvector.h"
#include "kmers/KmerSpectrumCore.h"
#include "paths/FindErrorsCore.h"
#include "paths/HyperBasevector.h"
#include "paths/KmerBaseBroker.h"
#include "paths/long/BubbleFreeN50.h"
#include "paths/long/BuildReadQGraph.h"
#include "paths/long/Correct1Pre.h"
#include "paths/long/FillPairs.h"
#include "paths/long/LoadCorrectCore.h"
#include "paths/long/LongProtoTools.h"
#include "paths/long/LongReadsToPaths.h"
#include "paths/long/ReadPath.h"
#include "paths/long/large/GapToyTools.h"

void Preclose( vecbvec const& bases, VecPQVec const& quals,
               const String& work_dir, HyperBasevector& hb, vec<int>& inv,
               ReadPathVec& paths, const Bool preclose_verbose ) {
    double clock = WallClockTime( );
    cout << "\n" << Date( ) << ": begin preclose" << std::endl;
    BubbleFreeN50 bub( hb, 1000 );
    PRINT3( hb.N( ), hb.EdgeObjectCount( ), bub.PostN50( ) );

    // Heuristics.

    const int min_links = 3;

    // Set up assembly.

    vec<int> to_left, to_right;
    hb.ToLeft(to_left), hb.ToRight(to_right);
    double yclock = WallClockTime( );
    vec< std::pair< vec<int>, vec<int> > > pairs( paths.size( ) );
    for ( int64_t i = 0; i < (int64_t) paths.size( ); i += 2 ) {
        int64_t id1 = i, id2 = i + 1;
        {
            vec<int> p1, p2;
            for ( int64_t j = 0; j < (int64_t) paths[id1].size( ); j++ )
                p1.push_back( paths[id1][j] );
            for ( int64_t j = (int64_t) paths[id2].size( ) - 1; j >= 0; j-- )
                p2.push_back( inv[ paths[id2][j] ] );
            pairs[i] = make_pair( p1, p2 );
        }
        {
            vec<int> p1, p2;
            for ( int64_t j = 0; j < (int64_t) paths[id2].size( ); j++ )
                p1.push_back( paths[id2][j] );
            for ( int64_t j = (int64_t) paths[id1].size( ) - 1; j >= 0; j-- )
                p2.push_back( inv[ paths[id1][j] ] );
            pairs[i+1] = make_pair( p1, p2 );
        }
    }
    cout << TimeSince(yclock) << " used defining pairs" << std::endl;

    // Phase 1.  Define edge pairs.  These are pairs of assembly edges (e1, e2)
    // which we think it would be a good idea to fill the gap between.  These are
    // called root edges.  The particular definition used here is as follows:
    // (a) Let e1 range over all terminal edges.
    // (b) Identify the path pairs (p1,p2) for which e1 is in p1 but not p2, and p2
    //     is nonempty.
    // (c) For a given e1 there must be at least min_links such pairs.
    // (d) Take all the edges that occur in the p2s and find one that occurs the
    //     largest number of times.  If that's at least min_links, call that
    //     edge e2.
    // (e) If both an edge pair and its inverse occur, delete one.

    double iclock = WallClockTime( );
    vec< std::pair<int,int> > edge_pairs;
    {
        vec<int> lends;
        for ( int e = 0; e < hb.EdgeObjectCount( ); e++ )
            if ( hb.From( to_right[e] ).empty( ) ) lends.push_back(e);
        vec< vec< vec<int> > > lfollow( lends.size( ) );
        for ( int i = 0; i < pairs.isize( ); i++ ) {
            const vec<int> &p1 = pairs[i].first, &p2 = pairs[i].second;
            for ( int j = 0; j < p1.isize( ); j++ ) {
                int x = BinPosition( lends, p1[j] );
                if ( x < 0 ) continue;
                if ( p2.empty( ) || Member( p2, p1[j] ) ) continue;
                lfollow[x].push_back(p2);
            }
        }
        // int count = 0;
        for ( int i = 0; i < lends.isize( ); i++ ) {
            if ( lfollow[i].isize( ) < min_links ) continue;
            Sort( lfollow[i] );
            vec<int> all;
            for ( int m = 0; m < lfollow[i].isize( ); m++ )
                all.append( lfollow[i][m] );
            Sort(all);
            int M = 0;
            int best = -1;
            for ( int r = 0; r < all.isize( ); r++ ) {
                int s = all.NextDiff(r);
                if ( s - r > M ) {
                    best = all[r];
                    M = s - r;
                }
                r = s - 1;
            }
            if ( M < min_links ) continue;
            edge_pairs.push( lends[i], best );
            /*
            cout << "[" << count++ << "] " << lends[i] << " --> ";
            for ( int m = 0; m < lfollow[i].isize( ); m++ )
            {    int n = lfollow[i].NextDiff(m);
                 if ( m > 0 ) std::cout << "; ";
                 cout << printSeq( lfollow[i][m] ) << "<" << n-m << ">";
                 m = n - 1;    }
            cout << "\n";
            */
        }
    }
    Sort(edge_pairs);
    vec<Bool> edel( edge_pairs.size( ), False );
    for ( int i = 0; i < edge_pairs.isize( ); i++ ) {
        pair<int,int> p = make_pair(
                              inv[ edge_pairs[i].second ], inv[ edge_pairs[i].first ] );
        if ( p < edge_pairs[i] && BinMember( edge_pairs, p ) ) edel[i] = True;
    }
    EraseIf( edge_pairs, edel );
    cout << TimeSince(iclock) << " used defining edge pairs" << std::endl;

    // Phase 2.  Define read pairs associated to each edge pair.  The particular
    // definition used here is as follows: we associate (p1,p2) to (e1,e2) if
    // e1 is in p1 but not p2, or e2 is in p2 but not p1.
    // NO: MORE STUFF ADDED, SEE BELOW.

    double xclock = WallClockTime( );
    vec< vec<int> > pids( edge_pairs.size( ) );
    {
        vec< vec<int> > pun1( pairs.size( ) ), pun2( pairs.size( ) );
        #pragma omp parallel for
        for ( int i = 0; i < pairs.isize( ); i++ ) {
            vec<int> p1 = pairs[i].first, p2 = pairs[i].second;
            Sort(p1), Sort(p2);
            for ( int j = 0; j < p1.isize( ); j++ )
                if ( !BinMember( p2, p1[j] ) ) pun1[i].push_back( p1[j] );
            for ( int j = 0; j < p2.isize( ); j++ )
                if ( !BinMember( p1, p2[j] ) ) pun2[i].push_back( p2[j] );
        }
        vec<vec<int>>  px1( hb.EdgeObjectCount( ) ), px2( hb.EdgeObjectCount( ) );
        for ( int i = 0; i < pairs.isize( ); i++ ) {
            for ( int j = 0; j < pun1[i].isize( ); j++ )
                px1[ pun1[i][j] ].push_back(i/2);
            for ( int j = 0; j < pun2[i].isize( ); j++ )
                px2[ pun2[i][j] ].push_back(i/2);
        }
        #pragma omp parallel for
        for ( int e = 0; e < hb.EdgeObjectCount( ); e++ ) {
            UniqueSort( px1[e] );
            UniqueSort( px2[e] );
        }
        #pragma omp parallel for
        for ( int j = 0; j < edge_pairs.isize( ); j++ ) {
            pids[j].append( px1[ edge_pairs[j].first ] );
            pids[j].append( px2[ edge_pairs[j].second ] );
            UniqueSort( pids[j] );
        }
    }
    cout << TimeSince(xclock) << " used selecting read pairs" << std::endl;

    // Add to pids.  This does not necessarily help.

    int nedges = hb.EdgeObjectCount( );
    vec< vec<int> > layout_pos(nedges);
    vec< vec<int64_t> > layout_id(nedges);
    vec< vec<Bool> > layout_or(nedges);
    LayoutReads( hb, inv, bases, paths, layout_pos, layout_id, layout_or );
    {
        vec< std::pair<vec<int>,vec<int>> > pairs;
        vec<int64_t> pairs_pid;
        DefinePairs( paths, inv, pairs, pairs_pid, work_dir + "/a.60" );
        // vec<vec<int>> left_empty, right_empty;
        // Empty( hb, pairs, pairs_pid, left_empty, right_empty, EMPTY2 );
        vec< std::pair<int,int> > pairsx;
        for ( int i = 0; i < pairs.isize( ); i++ ) {
            int ii = pairs.NextDiff(i);
            pairsx.push( i, ii );
            i = ii - 1;
        }
        vec<int> pairs_simple_start( nedges, -1 ), pairs_simple_stop( nedges, -1 );
        for ( int i = 0; i < pairsx.isize( ); i++ ) {
            int j1 = pairsx[i].first, j2 = pairsx[i].second;
            if ( pairs[j1].first.solo( ) && pairs[j1].first == pairs[j1].second ) {
                pairs_simple_start[ pairs[j1].first[0] ] = j1;
                pairs_simple_stop[ pairs[j1].first[0] ] = j2;
            }
        }
        for ( int ei = 0; ei < edge_pairs.isize( ); ei++ ) {
            const int lroot = edge_pairs[ei].first, rroot = edge_pairs[ei].second;
            const int max_delta = 120;
            vec<int> pids1, pids2;
            for ( int k = layout_pos[lroot].isize( ) - 1; k >= 0; k-- ) {
                int pos = layout_pos[lroot][k];
                int id = layout_id[lroot][k];
                if ( pos + bases[id].isize( )
                        < hb.EdgeLengthBases(lroot) - max_delta ) {
                    break;
                }
                pids1.push_back(id/2);
            }
            for ( int k = 0; k < layout_pos[rroot].isize( ); k++ ) {
                int pos = layout_pos[rroot][k];
                int id = layout_id[rroot][k];
                if ( pos > max_delta ) break;
                pids1.push_back(id/2);
            }
            UniqueSort(pids1);
            if ( pairs_simple_start[lroot] >= 0 ) {
                for ( int k = pairs_simple_start[lroot];
                        k < pairs_simple_stop[lroot]; k++ ) {
                    int64_t pid = pairs_pid[k];
                    for ( int l = 0; l < 2; l++ ) {
                        int id = 2*pid + l;
                        if ( paths[id].getOffset( ) + bases[id].isize( )
                                >= hb.EdgeLengthBases(lroot) - max_delta ) {
                            pids2.push_back(pid);
                        }
                    }
                }
            }
            if ( pairs_simple_start[rroot] >= 0 ) {
                for ( int k = pairs_simple_start[rroot];
                        k < pairs_simple_stop[rroot]; k++ ) {
                    int64_t pid = pairs_pid[k];
                    for ( int l = 0; l < 2; l++ ) {
                        int id = 2*pid + l;
                        if ( paths[id].getOffset( ) <= max_delta )
                            pids2.push_back(pid);
                    }
                }
            }
            UniqueSort(pids2);
            const int max_pids = 100;
            if ( pids1.isize( ) <= max_pids ) pids[ei].append(pids1);
            else pids[ei].append(pids2);
            UniqueSort( pids[ei] );
        }
    }

    // Phase 3.  Build K=60 assembly across each edge pair.  The particular method
    // used here is as follows:
    // (a) Form precorrected reads as in the generation of the mod0 reads in
    //     LongProto.
    // (b) Form the K=60 unipath graph associated to these reads and the root
    //     edges e1, e2.
    // (c) Find the edges containing the rightmost kmer of e1 and the leftmost
    //     kmer of e2, then find all edges between these.  Reduce the graph to
    //     this set of edges.

    int K = hb.K( );
    long_logging logc( "", "" );
    logc.STATUS_LOGGING = False;
    logc.MIN_LOGGING = False;
    long_heuristics heur( "" );
    vec<String> reports( edge_pairs.size( ) );
    vec<vecbasevector> extra( edge_pairs.size( ) );
    double pclock = WallClockTime( );
    #pragma omp parallel for
    for ( int ej = 0; ej < edge_pairs.isize( ); ej++ ) {
        ostringstream mout;
        mout << "\n[" << ej << "] " << edge_pairs[ej].first << " --> "
             << edge_pairs[ej].second << std::endl;
        PRINT_TO( mout, pids[ej].size( ) );

        vecbasevector lbases;
        vecqualvector lquals;
        for ( int i = 0; i < pids[ej].isize( ); i++ ) {
            int id1 = 2*pids[ej][i], id2 = 2*pids[ej][i] + 1;
            lbases.push_back( bases[id1] ), lbases.push_back( bases[id2] );
            lquals.push_back( quals.begin()[id1] );
            lquals.push_back( quals.begin()[id2] );
        }
        double clock = WallClockTime( );
        const int SEP = 0;
        const int STDEV = 100;
        const String LIB = "woof";
        const size_t nreads = bases.size( );
        PairsManager lpairs(nreads);
        lpairs.addLibrary( SEP, STDEV, LIB );
        size_t npairs = nreads / 2;
        for ( size_t pi = 0; pi < npairs; pi++ )
            lpairs.addPairToLib( 2 * pi, 2 * pi + 1, 0 );

        // Run Correct1.

        vec<int> trace_ids, precorrect_seq;
        ParseIntSet( "{" + heur.PRECORRECT_SEQ + "}", precorrect_seq, false );
        const int max_freq = heur.FF_MAX_FREQ;
        size_t nReads = lbases.size();

        {
            vecbasevector lbases0(lbases);
            PC_Params pcp;
            const int K_PC = 25;
            KmerSpectrum kspec(K_PC);
            pre_correct_parallel( pcp, K_PC, &lbases, &lquals, &kspec, -1, 1 );
            for ( int id = 0; id < (int) lbases.size( ); id++ ) {
                for ( int j = 0; j < lbases[id].isize( ); j++ ) {
                    if ( lbases[id][j] != lbases0[id][j] )
                        lquals[id][j] = 0;
                }
            }
        }

        // Carry out initial pair filling.

        vecbasevector lbases_done;
        vec<Bool> to_edit( nReads, True );
        vec<Bool> done( nReads, False );
        {
            lbases_done = lbases;
            vecbasevector filled;
            const int MIN_FREQ = 5;
            FillPairs( lbases, lpairs, MIN_FREQ, filled, heur.FILL_PAIRS_ALT,
                       False );
            int64_t fill_count = 0;
            for ( int64_t id = 0; id < (int64_t) filled.size( ); id++ ) {
                if ( filled[id].size( ) == 0 ) continue;
                fill_count++;
                int n = lbases[id].size( );
                lbases_done[id] = filled[id];
                lquals[id].resize(0);
                lquals[id].resize( filled[id].size( ), 40 );
                lbases[id] = lbases_done[id];
                if ( n < lbases[id].isize( ) ) {
                    lquals[id].resize(n);
                    if ( lpairs.getPartnerID(id) >= id ) lbases[id].resize(n);
                    else {
                        lbases[id].SetToSubOf( lbases[id],
                                               lbases[id].isize( ) - n, n );
                    }
                }
                done[id] = True;
                if ( lpairs.getPartnerID(id) < id ) lbases_done[id].resize(0);
                to_edit[id] = False;
            }
        }

        // Cap quality scores.

        CapQualityScores( lquals, done );

        // Do precorrection.

        vec<int> trim_to;
        for ( int j = 0; j < precorrect_seq.isize( ); j++ ) {
            String TMP = "does_not_exist";
            Correct1Pre( TMP, precorrect_seq[j], max_freq, lbases, lquals, lpairs,
                         to_edit, trim_to, trace_ids, logc, heur );
        }

        vecbasevector edges;
        edges.push_back( hb.EdgeObject( edge_pairs[ej].first ) );
        edges.push_back( hb.EdgeObject( edge_pairs[ej].second ) );
        vecbasevector all(lbases);
        all.Append(edges);
        HyperBasevector hbx;
        unsigned const COVERAGE = 50u;
        LongReadsToPaths(all, K, COVERAGE, False, False, &hbx);

        int start = -1, stop = -1;
        for ( int i = 0; i < hbx.EdgeObjectCount( ); i++ ) {
            String e = hbx.EdgeObject(i).ToString( );

            /*
            String l = edges[0].ToString( );
            l.resize(K);
            String r
                 = basevector( edges[1], edges[1].isize( ) - K, K ).ToString( );
            */

            String l
                = basevector( edges[0], edges[0].isize( ) - K, K ).ToString( );
            String r = edges[1].ToString( );
            r.resize(K);

            if ( e.Contains(l) ) start = i;
            if ( e.Contains(r) ) stop = i;
        }

        vec<int> to_left, to_right;
        hbx.ToLeft(to_left), hbx.ToRight(to_right);

        vec<int> bet = hbx.EdgesSomewhereBetween( to_right[start], to_left[stop] );
        if ( bet.empty( ) ) {
            reports[ej] = mout.str( );
            continue;
        }
        bet.push_back( start, stop );
        UniqueSort(bet);

        vec<int> dels;
        for ( int e = 0; e < hbx.EdgeObjectCount( ); e++ )
            if ( !BinMember( bet, e ) ) dels.push_back(e);
        hbx.DeleteEdges(dels);
        hbx.RemoveUnneededVertices( );
        hbx.RemoveDeadEdgeObjects( );

        PRINT_TO( mout, hbx.EdgeObjectCount( ) );

        for ( int i = 0; i < hbx.EdgeObjectCount( ); i++ )
            extra[ej].push_back( hbx.EdgeObject(i) );
        for ( int v = 0; v < hbx.N( ); v++ )
            for ( int l1 = 0; l1 < hbx.To(v).isize( ); l1++ )
                for ( int l2 = 0; l2 < hbx.From(v).isize( ); l2++ ) {
                    basevector b = TrimCat( K, hbx.EdgeObjectByIndexTo( v, l1 ),
                                            hbx.EdgeObjectByIndexFrom( v, l2 ) );
                    extra[ej].push_back(b);
                }
        reports[ej] = mout.str( );

        /*
        Ofstream( out, work_dir + "/glunk.dot" );
        hbx.PrintSummaryDOT0w( out, True, False, True );
        Ofstream( eout, work_dir + "/glunk.fasta" );
        for ( int e = 0; e < hbx.EdgeObjectCount( ); e++ )
        {   hbx.EdgeObject(e).Print( eout, e );    }
        Scram(0);
        */

    }
    cout << TimeSince(pclock) << " used in local assembly loop" << std::endl;

    if (preclose_verbose) {
        for ( int i = 0; i < reports.isize( ); i++ )
            cout << reports[i];
        cout << "\n";
    }

    // Phase 4. Rebuild global assembly.

    HyperBasevector hbp;
    {
        double aclock = WallClockTime( );
        vecbasevector all;
        for ( int pass = 1; pass <= 2; pass++ ) {
            int edges = 0;
            for ( int e = 0; e < hb.EdgeObjectCount( ); e++ ) {
                if ( inv[e] >= e ) {
                    if ( pass == 1 ) edges++;
                    else all.push_back( hb.EdgeObject(e) );
                }
            }
            for ( int v = 0; v < hb.N( ); v++ )
                for ( int l1 = 0; l1 < hb.To(v).isize( ); l1++ ) {
                    int x1 = hb.EdgeObjectIndexByIndexTo( v, l1 );
                    if ( pass == 1 ) {
                        for ( int l2 = 0; l2 < hb.From(v).isize( ); l2++ ) {
                            int x2 = hb.EdgeObjectIndexByIndexFrom( v, l2 );
                            if ( make_pair( x1, x2 )
                                    <= make_pair( inv[x2], inv[x1] ) ) {
                                edges++;
                            }
                        }
                    } else {
                        basevector e1 = hb.EdgeObject(x1);
                        e1 = basevector( e1, e1.isize( ) - K, K );
                        for ( int l2 = 0; l2 < hb.From(v).isize( ); l2++ ) {
                            int x2 = hb.EdgeObjectIndexByIndexFrom( v, l2 );
                            basevector e2 = hb.EdgeObject(x2);
                            if ( make_pair( x1, x2 )
                                    <= make_pair( inv[x2], inv[x1] ) ) {
                                e2.resize(K);
                                all.push_back(
                                    TrimCat( K, e1, e2 ) );
                            }
                        }
                    }
                }
            if ( pass == 1 ) {
                for ( int j = 0; j < edge_pairs.isize( ); j++ )
                    edges += extra[j].size( );
                all.reserve(edges);
            } else {
                for ( int j = 0; j < edge_pairs.isize( ); j++ )
                    all.Append( extra[j] );
            }
        }
        cout << TimeSince(aclock) << " used building all" << std::endl;
        double uclock = WallClockTime( );
        unsigned const COVERAGE = 50u;
        LongReadsToPaths(all, K, COVERAGE, False, False, &hbp);
        cout << TimeSince(uclock) << " used in unipath graph reconstruction"
             << std::endl;
    }
    BubbleFreeN50 bubp( hbp, 1000 );
    PRINT3( hbp.N( ), hbp.EdgeObjectCount( ), bubp.PostN50( ) );

    // Reconstruct paths.

    double rclock = WallClockTime( );
    HyperBasevector hbn;
    ReadPathVec pathsn;
    rePath( hbp, bases, quals, true, false, &hbn, &pathsn );
    BubbleFreeN50 bubn( hbn, 1000 );
    PRINT3( hbn.N( ), hbn.EdgeObjectCount( ), bubn.PostN50( ) );
    cout << TimeSince(rclock) << " used reconstructing paths" << std::endl;
    double tclock = WallClockTime( );

    vec<int> invn;
    hbn.Involution(invn);

    hb = hbn;
    inv = invn;
    paths = pathsn;
    cout << TimeSince(tclock) << " used in tail" << std::endl;
    double wclock = WallClockTime( );
    cout << Date( ) << ": done with preclose" << std::endl;
    cout << TimeSince(clock) << " used in total in preclose\n" << std::endl;
}
