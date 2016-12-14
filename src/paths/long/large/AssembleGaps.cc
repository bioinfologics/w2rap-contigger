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
void MakeStartStop(const std::vector<basevector> &bell,
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
    void operator( )(const std::vector<basevector>  &bell, const HyperBasevector &hb,
                     const HyperBasevector &shb, const vec<int> &lefts,
                     const vec<int> &rights, vec<int> &starts, vec<int> &stops) {
        MakeStartStop<M>(bell, hb, shb, lefts, rights, starts, stops);
    }
};

void FindPidsST(std::vector<int64_t> & pids, const vec<int> &lefts, const vec<int> &rights,
                const std::vector<std::vector<int>> &layout_pos, const std::vector<std::vector<int64_t>> &layout_id, const std::vector<std::vector<bool>> & layout_or,
                const int MAX_PROX_LEFT, const int MAX_PROX_RIGHT, const int pair_sample){


    // First find the pairs that bridge from left to right, and mark
    // their endpoints.  Inefficient.

    std::vector<int64_t> pids1;
    vec<vec<int>> lstarts(lefts.size()), rstarts(rights.size());
    {
        vec<quad<int64_t, Bool, int, int> > marks;
        for (int l = 0; l < lefts.isize(); l++)
            for (int k = 0; k < layout_pos[lefts[l]].size(); k++) {
                if (!layout_or[lefts[l]][k]) continue;
                int pos = layout_pos[lefts[l]][k];
                int64_t id = layout_id[lefts[l]][k];
                marks.push(id / 2, False, pos, l);
            }
        for (int l = 0; l < rights.isize(); l++)
            for (int k = 0; k < layout_pos[rights[l]].size(); k++) {
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
    //UniqueSort(pids1);
    std::sort(pids1.begin(), pids1.end());
    pids1.erase(std::unique(pids1.begin(), pids1.end()), pids1.end());
    for (int l = 0; l < lefts.isize(); l++)
        Sort(lstarts[l]);
    for (int l = 0; l < rights.isize(); l++)
        Sort(rstarts[l]);

    // Now find the pairs that start close to one of the bridge pairs.

    std::vector<int64_t> pids2;
    for (int l = 0; l < lefts.isize(); l++)
        for (int k = 0; k < layout_pos[lefts[l]].size(); k++) {
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
        for (int k = 0; k < layout_pos[rights[l]].size(); k++) {
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

    //UniqueSort(pids2);
    std::sort(pids2.begin(), pids2.end());
    pids2.erase(std::unique(pids2.begin(), pids2.end()), pids2.end());

    // Now subsample if needed.

    int keep = pair_sample / 2;
    pids.reserve(std::min((int) (pids1.size() + pids2.size()),pair_sample));

    //Pids1
    if (pids1.size() + pids2.size() <= pair_sample or pids1.size() <= keep)
        pids.insert(pids.end(), pids1.begin(), pids1.end());
    else {
        for (int l = 0; l < keep; l++) {
            int m = (l * pids1.size()) / keep;
            pids.push_back(pids1[m]);
        }
    }

    //Pids2
    if (pids.size() + pids2.size() <= pair_sample or pids2.size() <= keep)
        pids.insert(pids.end(), pids2.begin(), pids2.end());
    else {
        for (int l = 0; l < keep; l++) {
            int m = (l * pids2.size()) / keep;
            pids.push_back(pids2[m]);
        }
    }

    //UniqueSort(pids);
    std::sort(pids.begin(),pids.end());
    pids.erase(std::unique(pids.begin(),pids.end()),pids.end());

}

void CreateLocalReadSet(vecbasevector &gbases,vecqualvector &gquals, PairsManager &gpairs, std::vector<int64_t> & pids,
                        const vecbasevector &bases, VecPQVec const &quals) {


        qvec qv;
        gbases.reserve(2 * pids.size());
        gquals.reserve(2 * pids.size());

        for (int l = 0; l < pids.size(); l++) {
            int64_t pid = pids[l];
            int64_t id1 = 2 * pid, id2 = 2 * pid + 1;
            gbases.push_back(bases[id1]);
            gbases.push_back(bases[id2]);
            quals[id1].unpack(&qv); //TODO it is probably better to move this to the task!
            gquals.push_back(qv);
            quals[id2].unpack(&qv);
            gquals.push_back(qv);
        }

        const int SEP = 0;
        const int STDEV = 100;
        const String LIB = "woof";
        const size_t nreads = gbases.size();
        // PairsManager gpairs(nreads);

        gpairs = PairsManager(nreads);
        gpairs.addLibrary(SEP, STDEV, LIB);
        size_t npairs = nreads / 2;
        for (size_t pi = 0; pi < npairs; pi++) gpairs.addPairToLib(2 * pi, 2 * pi + 1, 0);
}

void AssembleGaps2(HyperBasevector &hb, vec<int> &inv2, ReadPathVec &paths2,
                   VecULongVec &paths2_index, const vecbasevector &bases, VecPQVec const &quals,
                   const String &work_dir, std::vector<int> k2floor_sequence,
                   vecbvec &new_stuff, const Bool CYCLIC_SAVE,
                   const int A2V, const int MAX_PROX_LEFT,
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

    __gnu_parallel::sort(LR.begin(), LR.end());

    // Remove inverted copies.  Should force symmetry first.

    {
        vec<Bool> lrd(LR.size(), False);
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
    }
    std::cout << Date() << ": " << LR.size() << " non-inverted clusters" << std::endl;
    // Some setup stuff.

    int nedges = hb.EdgeObjectCount();
    int K = hb.K();

    // Layout reads.  Expensive, temporary (?).

    std::vector<std::vector<int> > layout_pos(nedges);
    std::vector<std::vector<int64_t> > layout_id(nedges);
    std::vector<std::vector<bool>> layout_or(nedges);
    LayoutReads(hb, inv2, bases, paths2, layout_pos, layout_id, layout_or);

    // Make gap assemblies.

    vec<int> mgc = {2};
    int tmpdir_serial = 0;
    int min_gap_count = mgc[0], nobj = hb.EdgeObjectCount();
    vec<int> to_left, to_right;
    hb.ToLeft(to_left), hb.ToRight(to_right);
    vec<vec<basevector> > extras(LR.size());//this is accumulation, generates memory blocks
    vec<HyperBasevector> mhbp(LR.size());//this is accumulation, generates memory blocks
    std::cout << Date() << ": processing " << LR.size() << " blobs" << std::endl;
    double clockp1 = WallClockTime();
    int nblobs = LR.size();
    std::atomic_uint_fast64_t solved(0);

    //TODO: check local variable usage, should be made minimal!!!
    //Init readstacks, we'll need them!
    readstack::init_LUTs();
#define BATCH_SIZE 5000

    for (uint64_t bstart = 0; bstart < nblobs; bstart += BATCH_SIZE) {
        #pragma omp parallel
        {
            #pragma omp for schedule(dynamic,1)
            for (uint64_t bl = bstart; bl < bstart + BATCH_SIZE; ++bl) {
                if (bl >= nblobs) continue;
                //First part: create the gbases and gquals. this is locked by memory accesses and very convoluted
                const vec<int> &lefts = LR[bl].first, &rights = LR[bl].second; //TODO: how big is this? can we copy it?

                //Local readset
                std::vector<int64_t> pids;
                vecbasevector gbases;
                vecqualvector gquals;
                PairsManager gpairs;

                //Corrected reads
                VecEFasta corrected;
                vecbasevector creads;
                vec<pairing_info> cpartner;
                vec<int> cid;

                //Local assembly graph
                HyperBasevector xshb;

                //PART1-------------------------------
                FindPidsST(pids, lefts, rights, layout_pos, layout_id, layout_or, MAX_PROX_LEFT, MAX_PROX_RIGHT,
                           pair_sample);


                CreateLocalReadSet(gbases, gquals, gpairs, pids, bases, quals);
                HyperBasevector *mhbp_t = &mhbp[bl];

                //#pragma omp task shared(lefts,rights)
                //{
                uint NUM_THREADS = 1;
                long_heuristics heur("");
                heur.K2_FLOOR = k2floor_sequence[0];
                CorrectionSuite(gbases, gquals, gpairs, heur, creads, corrected, cid, cpartner, NUM_THREADS, "",
                                False);

                for (auto K2_FLOOR_LOCAL: k2floor_sequence) {
                    SupportedHyperBasevector shb;

                    MakeLocalAssembly2(corrected, lefts, rights, shb, K2_FLOOR_LOCAL, creads, cid, cpartner);

                    if (shb.K() == 0) continue;

                    // Find edges "starts" and "stops" overlapping root edges.
                    vec<int> starts, stops;
                    std::vector<basevector> bell;
                    bell.reserve(shb.EdgeObjectCount() + lefts.isize() + rights.isize());
                    for (int e = 0; e < shb.EdgeObjectCount(); e++)
                        bell.push_back(shb.EdgeObject(e));
                    for (int l = 0; l < lefts.isize(); l++)
                        bell.push_back(hb.EdgeObject(lefts[l]));
                    for (int r = 0; r < rights.isize(); r++)
                        bell.push_back(hb.EdgeObject(rights[r]));

                    BigK::dispatch<MakeStartStopFunctor>(shb.K(), bell, hb, shb, lefts, rights, starts, stops);

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
                    xshb = shb;
                    xshb.DeleteEdges(ydels);
                    xshb.RemoveUnneededVertices();
                    xshb.RemoveDeadEdgeObjects();


                    if (!CYCLIC_SAVE || xshb.Acyclic()) break;

                }

                if (xshb.Acyclic() && xshb.N() > 0) {

                    // Make bpaths.  These are all source-sink paths through the
                    // local graph.

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

                    if (bpaths.isize() <= MAX_BPATHS) {
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
                        // Make the bpaths into a HyperBasevector.
                        vecbasevector bpathsx;
                        for (int l = 0; l < bpaths.isize(); l++)
                            bpathsx.push_back(bpaths[l]);
                        BasesToGraph(bpathsx, K, *mhbp_t);
                        ++solved;
                    }
                }
                //}//---OMP TASK END---
            }
        }

        std::cout << Date() << ": "<< std::min(bstart+BATCH_SIZE,(uint64_t)nblobs) <<" blobs processed, paths found for " << solved << std::endl;
    }
    std::cout << Date() << TimeSince(clockp1) << " spent in local assemblies." << std::endl;

    TIMELOG_REPORT(std::cout,AssembleGaps,AG2_FindPids,AG2_ReadSetCreation,AG2_CorrectionSuite,AG2_LocalAssembly2,AG2_LocalAssemblyEval,AG2_CreateBpaths,AG2_PushBpathsToGraph);
    TIMELOG_REPORT(std::cout,Correct1Pre,C1P_Align,C1P_InitBasesQuals,C1P_Correct,C1P_UpdateBasesQuals);
    TIMELOG_REPORT(std::cout,CorrectPairs1,CP1_Align,CP1_MakeStacks,CP1_Correct);
    // Do the patching.
    const vec<std::pair<int, int> > blobs(LR.size());
    Patch(hb, blobs, mhbp, work_dir, new_stuff);
}
