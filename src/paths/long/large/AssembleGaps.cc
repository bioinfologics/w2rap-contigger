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
//#include "ParallelVecUtilities.h"
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
#include <util/w2rap_timers.h>
#include <paths/long/LoadCorrectCore.h>
#include <paths/long/ReadStack.h>


template<int M>
void MakeStartStop(const vecbasevector &bell,
                   const HyperBasevector &hb, const HyperBasevector &shb, const vec<int> &lefts,
                   const vec<int> &rights, vec<int> &starts, vec<int> &stops) {
    vec<triple<kmer<M>, int, int> > kmers_plus;
    MakeKmerLookup3(bell, kmers_plus);
    for (int64_t i = 0; i < kmers_plus.jsize(); i++) {
        int64_t j;
        for (j = i + 1; j < kmers_plus.jsize(); j++)
            if (kmers_plus[j].first != kmers_plus[i].first) break;
        vec<int> es;
        Bool have_l = False, have_r = False;
        for (int64_t k = i; k < j; k++) {
            int id = kmers_plus[k].second;
            if (id < shb.EdgeObjectCount()) es.push_back(id);
            else {
                id -= shb.EdgeObjectCount();
                if (id < lefts.isize()) have_l = True;
                else have_r = True;
            }
        }
        if (have_l) starts.append(es);
        if (have_r) stops.append(es);
        i = j - 1;
    }
}

template<int M>
struct MakeStartStopFunctor {
    void operator( )(const vecbasevector &bell, const HyperBasevector &hb,
                     const HyperBasevector &shb, const vec<int> &lefts,
                     const vec<int> &rights, vec<int> &starts, vec<int> &stops) {
        MakeStartStop<M>(bell, hb, shb, lefts, rights, starts, stops);
    }
};

void FindPidsST(vec<int64_t> & pids, const vec<int> &lefts, const vec<int> &rights,
                const vec<vec<int>> &layout_pos, const vec<vec<int64_t>> &layout_id, const vec<vec<Bool>> & layout_or,
                const int MAX_PROX_LEFT, const int MAX_PROX_RIGHT, const int pair_sample){

    // Get ready.



    {
        // Heuristics.


        // First find the pairs that bridge from left to right, and mark
        // their endpoints.  Inefficient.

        vec<int64_t> pids1;
        vec<vec<int>> lstarts(lefts.size()), rstarts(rights.size());
        {
            vec<quad<int64_t, Bool, int, int> > marks;
            for (int l = 0; l < lefts.isize(); l++)
                for (int k = 0; k < layout_pos[lefts[l]].isize(); k++) {
                    if (!layout_or[lefts[l]][k]) continue;
                    int pos = layout_pos[lefts[l]][k];
                    int64_t id = layout_id[lefts[l]][k];
                    marks.push(id / 2, False, pos, l);
                }
            for (int l = 0; l < rights.isize(); l++)
                for (int k = 0; k < layout_pos[rights[l]].isize(); k++) {
                    if (layout_or[rights[l]][k]) continue;
                    int pos = layout_pos[rights[l]][k];
                    int64_t id = layout_id[rights[l]][k];
                    marks.push(id / 2, True, pos, l);
                }
            Sort(marks);
            for (int l = 0; l < marks.isize(); l++) {
                int m;
                for (m = l + 1; m < marks.isize(); m++)
                    if (marks[m].first != marks[l].first) break;
                Bool left = False, right = False;
                for (int j = l; j < m; j++) {
                    if (!marks[j].second) left = True;
                    else right = True;
                }
                if (left && right) {
                    pids1.push_back(marks[l].first);
                    for (int k = l; k < m; k++) {
                        if (!marks[k].second) {
                            lstarts[marks[k].fourth].push_back(
                                    marks[k].third);
                        } else {
                            rstarts[marks[k].fourth].push_back(
                                    marks[k].third);
                        }
                    }
                }
                l = m - 1;
            }
        }
        UniqueSort(pids1);
        for (int l = 0; l < lefts.isize(); l++)
            Sort(lstarts[l]);
        for (int l = 0; l < rights.isize(); l++)
            Sort(rstarts[l]);

        // Now find the pairs that start close to one of the bridge pairs.

        vec<int64_t> pids2;
        for (int l = 0; l < lefts.isize(); l++)
            for (int k = 0; k < layout_pos[lefts[l]].isize(); k++) {
                int pos = layout_pos[lefts[l]][k];
                int64_t id = layout_id[lefts[l]][k];
                Bool fw = layout_or[lefts[l]][k];
                if (BinMember(pids1, id / 2)) continue;

                int low = lstarts[l].front(), high = lstarts[l].back();
                Bool close = False;
                if (low <= pos && pos <= high) close = True;
                else {
                    if (fw) {
                        if (low > pos && low - pos <= MAX_PROX_LEFT) { close = True; }
                        else if (high < pos && pos - high <= MAX_PROX_RIGHT) { close = True; }
                    } else {
                        if (low > pos && low - pos <= MAX_PROX_RIGHT) { close = True; }
                        else if (high < pos && pos - high <= MAX_PROX_LEFT) { close = True; }
                    }
                }
                if (close) pids2.push_back(id / 2);
            }

        for (int l = 0; l < rights.isize(); l++)
            for (int k = 0; k < layout_pos[rights[l]].isize(); k++) {
                int pos = layout_pos[rights[l]][k];
                int64_t id = layout_id[rights[l]][k];
                Bool fw = layout_or[rights[l]][k];
                if (BinMember(pids1, id / 2)) continue;

                int low = rstarts[l].front(), high = rstarts[l].back();
                Bool close = False;
                if (low <= pos && pos <= high) close = True;
                else {
                    if (fw) {
                        if (low > pos && low - pos <= MAX_PROX_LEFT) { close = True; }
                        else if (high < pos && pos - high <= MAX_PROX_RIGHT) { close = True; }
                    } else {
                        if (low > pos && low - pos <= MAX_PROX_RIGHT) { close = True; }
                        else if (high < pos && pos - high <= MAX_PROX_LEFT) { close = True; }
                    }
                }
                if (close) pids2.push_back(id / 2);
            }

        UniqueSort(pids2);

        // Now subsample if needed.

        int keep = pair_sample / 2;
        if (pids1.isize() + pids2.isize() <= pair_sample)
            pids.append(pids1);
        else if (pids1.isize() <= keep) pids1.append(pids1);
        else {
            for (int l = 0; l < keep; l++) {
                int m = (l * pids1.isize()) / keep;
                pids.push_back(pids1[m]);
            }
        }
        if (pids.isize() + pids2.isize() <= pair_sample)
            pids.append(pids2);
        else {
            keep = pair_sample - pids.isize();
            for (int l = 0; l < keep; l++) {
                int m = (l * pids2.isize()) / keep;
                pids.push_back(pids2[m]);
            }
        }
    }
    UniqueSort(pids);
}


void AssembleGaps2(HyperBasevector &hb, vec<int> &inv2, ReadPathVec &paths2,
                   VecULongVec &paths2_index, const vecbasevector &bases, VecPQVec const &quals,
                   const String &work_dir, std::vector<int> k2floor_sequence,
                   vecbvec &new_stuff, const Bool CYCLIC_SAVE,
                   const int A2V, const int GAP_CAP, const int MAX_PROX_LEFT,
                   const int MAX_PROX_RIGHT, const int MAX_BPATHS, const int pair_sample) {
    // Find clusters of unsatisfied links.

    vec<vec<std::pair<int, int> > > xs;
    Unsat(hb, inv2, paths2, xs, work_dir, A2V);

    //Should print some stats about Unsats here.

    // Condense to lists of lefts and rights.

    vec<std::pair<vec<int>, vec<int>>> LR(xs.size());
    #pragma omp parallel for
    for (auto i = 0; i < xs.size(); i++) {
        vec<int> lefts, rights;
        for (int j = 0; j < xs[i].size(); j++) {
            lefts.push_back(xs[i][j].first);
            rights.push_back(xs[i][j].second);
        }
        UniqueSort(lefts), UniqueSort(rights);
        LR[i] = make_pair(lefts, rights);
    }

    //sortInPlaceParallel(LR.begin(), LR.end());
    __gnu_parallel::sort(LR.begin(),LR.end());

    // Remove inverted copies.  Should force symmetry first.

    {
        vec<Bool> lrd(LR.size(), False);
        PRINT(LR.size());
        #pragma omp parallel for
        for (int i = 0; i < LR.isize(); i++) {
            vec<int> lefts, rights;
            for (int j = 0; j < LR[i].first.isize(); j++)
                rights.push_back(inv2[LR[i].first[j]]);
            for (int j = 0; j < LR[i].second.isize(); j++)
                lefts.push_back(inv2[LR[i].second[j]]);
            Sort(lefts), Sort(rights);
            if (make_pair(lefts, rights) <= LR[i]) continue;
            if (!BinMember(LR, make_pair(lefts, rights))) continue;
            lrd[i] = True;
        }
        EraseIf(LR, lrd);
        PRINT(LR.size());
    }

    // Some setup stuff.

    int nedges = hb.EdgeObjectCount();
    int K = hb.K();

    // Layout reads.  Expensive, temporary (?).

    vec<vec<int> > layout_pos(nedges);
    vec<vec<int64_t> > layout_id(nedges);
    vec<vec<Bool> > layout_or(nedges);
    LayoutReads(hb, inv2, bases, paths2, layout_pos, layout_id, layout_or);

    // Make gap assemblies.

    vec<int> mgc = {2};
    int tmpdir_serial = 0;
    int min_gap_count = mgc[0], nobj = hb.EdgeObjectCount();
    vec<int> to_left, to_right;
    hb.ToLeft(to_left), hb.ToRight(to_right);
    vec<vec<basevector> > extras(LR.size());//this is accumulation, generates memory blocks
    vec<HyperBasevector> mhbp(LR.size());//this is accumulation, generates memory blocks
    std::cout << Date() << ": now processing " << LR.size() << " blobs" << std::endl;
    std::cout << Date() << ": memory in use = " << MemUsageGBString()
              #ifdef __linux
              << ", peak = " << PeakMemUsageGBString( )
              #endif
              << std::endl;
    double clockp1 = WallClockTime();
    int nblobs = LR.size(), dots_printed = 0, nprocessed = 0;
    int lrc = LR.size();
    if (GAP_CAP >= 0) lrc = GAP_CAP;



    std::cout << "And now for the really slow part..." << std::endl;
    //TODO: check local variable usage, should be made minimal!!!
    TIMELOG_DECLARE_ATOMIC(AG2_FindPids);
    TIMELOG_DECLARE_ATOMIC(AG2_ReadSetCreation);
    TIMELOG_DECLARE_ATOMIC(AG2_CorrectionSuite);
    TIMELOG_DECLARE_ATOMIC(AG2_LocalAssembly2);
    TIMELOG_DECLARE_ATOMIC(AG2_LocalAssemblyEval);
    TIMELOG_DECLARE_ATOMIC(AG2_CreateBpaths);
    TIMELOG_DECLARE_ATOMIC(AG2_PushBpathsToGraph);
    //Init readstacks, we'll need them!
    readstack::init_LUTs();

    #pragma omp parallel for schedule(dynamic, 5)
    for (int bl = lrc-1; bl >= 0; --bl) {

        TIMELOG_DECLARE_LOCAL(AG2_FindPids,Loop);
        TIMELOG_DECLARE_LOCAL(AG2_ReadSetCreation,Loop);
        TIMELOG_DECLARE_LOCAL(AG2_CorrectionSuite,Loop);
        TIMELOG_DECLARE_LOCAL(AG2_LocalAssembly2,Loop);
        TIMELOG_DECLARE_LOCAL(AG2_LocalAssemblyEval,Loop);
        TIMELOG_DECLARE_LOCAL(AG2_CreateBpaths,Loop);
        TIMELOG_DECLARE_LOCAL(AG2_PushBpathsToGraph,Loop);


        vec<int64_t> pids;
        const vec<int> &lefts = LR[bl].first, &rights = LR[bl].second; //TODO: how big is this? can we copy it?


        //PART1-------------------------------
        TIMELOG_START_LOCAL(AG2_FindPids,Loop);

        FindPidsST(pids,lefts,rights,layout_pos,layout_id,layout_or,MAX_PROX_LEFT,MAX_PROX_RIGHT,pair_sample);

        TIMELOG_STOP_LOCAL(AG2_FindPids,Loop);

        //PART2-----------------------------------------------------
        TIMELOG_START_LOCAL(AG2_ReadSetCreation,Loop);


        //Assembly starts

        VecEFasta corrected;
        vecbasevector creads;
        vec<pairing_info> cpartner;
        vec<int> cid;
        SupportedHyperBasevector shb;

        //Local readset creation and error correction.
        vecbasevector gbases;
        vecqualvector gquals;
        PairsManager gpairs;

        qvec qv;
        gbases.reserve(2*pids.isize());
        gquals.reserve(2*pids.isize());

        for ( int l = 0; l < pids.isize( ); l++ )
        {    int64_t pid = pids[l];
            int64_t id1 = 2*pid, id2 = 2*pid + 1;
            gbases.push_back( bases[id1] );
            gbases.push_back( bases[id2] );
            quals[id1].unpack(&qv);
            gquals.push_back( qv );
            quals[id2].unpack(&qv);
            gquals.push_back( qv );    }

        const int SEP = 0;
        const int STDEV = 100;
        const String LIB = "woof";
        const size_t nreads = gbases.size( );
        // PairsManager gpairs(nreads);

        gpairs = PairsManager(nreads);
        gpairs.addLibrary( SEP, STDEV, LIB );
        size_t npairs = nreads / 2;
        for ( size_t pi = 0; pi < npairs; pi++ ) gpairs.addPairToLib( 2 * pi, 2 * pi + 1, 0 );
        TIMELOG_STOP_LOCAL(AG2_ReadSetCreation,Loop);

        TIMELOG_START_LOCAL(AG2_CorrectionSuite,Loop);

        uint NUM_THREADS = 1;
        long_heuristics heur( "" ); //TODO: this is allocated in stack and wastes both time and space!
        heur.K2_FLOOR = k2floor_sequence[0];
        CorrectionSuite( gbases, gquals, gpairs, heur, creads, corrected, cid, cpartner, NUM_THREADS, "", False);



        TIMELOG_STOP_LOCAL(AG2_CorrectionSuite,Loop);
        //PART3-----------------------------------------------------
        HyperBasevector xshb;

        for (auto K2_FLOOR_LOCAL: k2floor_sequence) {
            TIMELOG_START_LOCAL(AG2_LocalAssembly2, Loop);
            MakeLocalAssembly2(corrected, lefts, rights, shb, K2_FLOOR_LOCAL, creads, cid, cpartner);
            TIMELOG_STOP_LOCAL(AG2_LocalAssembly2, Loop);

            if (shb.K() == 0) continue;

            TIMELOG_START_LOCAL(AG2_LocalAssemblyEval, Loop);

            // Find edges "starts" and "stops" overlapping root edges.
            vec<int> starts, stops;
            vecbasevector bell;
            for (int e = 0; e < shb.EdgeObjectCount(); e++)
                bell.push_back(shb.EdgeObject(e));
            for (int l = 0; l < lefts.isize(); l++)
                bell.push_back(hb.EdgeObject(lefts[l]));
            for (int r = 0; r < rights.isize(); r++)
                bell.push_back(hb.EdgeObject(rights[r]));

            BigK::dispatch<MakeStartStopFunctor>( shb.K(), bell, hb, shb, lefts, rights, starts, stops);

            UniqueSort(starts), UniqueSort(stops);

            // Reduce shb to those edges between starts and stops.

            vec<int> yto_left, yto_right;
            shb.ToLeft(yto_left), shb.ToRight(yto_right);
            vec<int> keep = Intersection(starts, stops);
            keep.append(starts);
            keep.append(stops);
            for (int j1 = 0; j1 < starts.isize(); j1++)
                for (int j2 = 0; j2 < stops.isize(); j2++) {
                    int v = yto_right[starts[j1]], w = yto_left[stops[j2]];
                    vec<int> b = shb.EdgesSomewhereBetween(v, w);
                    keep.append(b);
                }
            UniqueSort(keep);
            vec<int> ydels;
            for (int e = 0; e < shb.EdgeObjectCount(); e++)
                if (!BinMember(keep, e)) ydels.push_back(e);
            xshb=shb;
            xshb.DeleteEdges(ydels);
            xshb.RemoveUnneededVertices();
            xshb.RemoveDeadEdgeObjects();
            //mout << TimeSince(sclock) << " used contracting" << std::endl;

            TIMELOG_STOP_LOCAL(AG2_LocalAssemblyEval,Loop);

            if (!CYCLIC_SAVE || xshb.Acyclic() ) break;

        };

        //PART4-----------------------------------------------------
        if (!xshb.Acyclic() || xshb.N() == 0) {
            continue;
        }

        //mout << "local assembly has " << xshb.NComponents( )
        //     << " components" << "\n";

        // Make bpaths.  These are all source-sink paths through the
        // local graph.
        TIMELOG_START_LOCAL(AG2_CreateBpaths,Loop);

        double aclock2 = WallClockTime();
        vec<basevector> bpaths;
        vec<int> sources, sinks;
        xshb.Sources(sources), xshb.Sinks(sinks);
        vec<int> zto_left, zto_right;
        xshb.ToLeft(zto_left), xshb.ToRight(zto_right);
        for (int i1 = 0; i1 < sources.isize(); i1++) {
            for (int i2 = 0; i2 < sinks.isize(); i2++) {
                vec<vec<int>> p;
                xshb.EdgePaths(zto_left, zto_right, sources[i1], sinks[i2], p);
                for (int l = 0; l < p.isize(); l++) {
                    basevector b = xshb.EdgeObject(p[l][0]);
                    for (int m = 1; m < p[l].isize(); m++) {
                        b.resize(b.isize() - (xshb.K() - 1));
                        b = Cat(b, xshb.EdgeObject(p[l][m]));
                    }
                    bpaths.push_back(b);
                    if (bpaths.isize() > MAX_BPATHS) break;
                }
            }
        }

        //PRINT_TO( mout, bpaths.size( ) );
        if (bpaths.isize() > MAX_BPATHS) {    //mout << "Too many bpaths." << std::endl;
            //mreport[bl] += mout.str( );
            //Dot( nblobs, nprocessed, dots_printed, ANNOUNCE, bl );
            TIMELOG_STOP_LOCAL(AG2_CreateBpaths,Loop);
            continue;
        }

        // Make more bpaths.

        for (int l = 0; l < lefts.isize(); l++) {
            Bool ext = False;
            for (int m = 0; m < lefts.isize(); m++) {
                if (to_right[lefts[m]] == to_left[lefts[l]]) {
                    basevector b = hb.EdgeObject(lefts[m]);
                    b.resize(b.isize() - (K - 1));
                    b = Cat(b, hb.EdgeObject(lefts[l]));
                    bpaths.push_back(b);
                    ext = True;
                }
            }
            if (!ext) bpaths.push_back(hb.EdgeObject(lefts[l]));
        }
        for (int r = 0; r < rights.isize(); r++) {
            Bool ext = False;
            for (int m = 0; m < rights.isize(); m++) {
                if (to_left[rights[m]] == to_right[rights[r]]) {
                    basevector b = hb.EdgeObject(rights[r]);
                    b.resize(b.size() - (K - 1));
                    b = Cat(b, hb.EdgeObject(rights[m]));
                    bpaths.push_back(b);
                    ext = True;
                }
            }
            if (!ext) bpaths.push_back(hb.EdgeObject(rights[r]));
        }

        TIMELOG_STOP_LOCAL(AG2_CreateBpaths,Loop);
        // Make the bpaths into a HyperBasevector.

        TIMELOG_START_LOCAL(AG2_PushBpathsToGraph,Loop);
        vecbasevector bpathsx;
        for (int l = 0; l < bpaths.isize(); l++)
            bpathsx.push_back(bpaths[l]);
        BasesToGraph(bpathsx, K, mhbp[bl]);
        TIMELOG_STOP_LOCAL(AG2_PushBpathsToGraph,Loop);

    }
    std::cout << " ... finally finished."<< std::endl;
    std::cout << TimeSince(clockp1) << " spent in local assemblies, "
              << "memory in use = " << MemUsageGBString()
              #ifdef __linux
              << ", peak = " << PeakMemUsageGBString( )
              #endif
              << std::endl;

    TIMELOG_REPORT(std::cout,AssembleGaps,AG2_FindPids,AG2_ReadSetCreation,AG2_CorrectionSuite,AG2_LocalAssembly2,AG2_LocalAssemblyEval,AG2_CreateBpaths,AG2_PushBpathsToGraph);
    TIMELOG_REPORT(std::cout,Correct1Pre,C1P_Align,C1P_InitBasesQuals,C1P_Correct,C1P_UpdateBasesQuals);
    TIMELOG_REPORT(std::cout,CorrectPairs1,CP1_Align,CP1_MakeStacks,CP1_Correct);
    // Do the patching.
    const vec<std::pair<int, int> > blobs(LR.size());
    Patch(hb, blobs, mhbp, work_dir, new_stuff);
}
