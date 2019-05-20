#include <iostream>
#include <Vec.h>
#include <util/OutputLog.h>
#include <CoreTools.h>
#include <feudal/BinaryStream.h>
#include <util/OutputLog.h>
#include "CoreTools.h"
#include "Equiv.h"
//#include "ParallelVecUtilities.h"
#include "ParseSet.h"
#include "VecUtilities.h"
#include "paths/HyperBasevector.h"
#include "paths/long/MakeKmerStuff.h"
#include "paths/long/ReadPath.h"
#include "paths/long/large/Unsat.h"
#include "system/SortInPlace.h"

void time_edge_paths(const HyperBasevector &xshb, const vec<int> &zto_left, const vec<int> &zto_right, const vec<int> &sources, const vec<int> &sinks, vec<vec<int>> &paths, int MAX_BPATHS = 100000) {
    auto start = std::chorno::high_resolution_clock::now();

    xshb.EdgePaths(zto_left, zto_right, sources[i1], sinks[i2], paths, -1, MAX_BPATHS);
    auto stop = std::chorno::high_resolution_clock::now();

    for (const auto &path:paths) {
        min_path_size = std::min(min_path_size, path.size());
        max_path_size = std::max(max_path_size, path.size());
        avg_path_size += path.size();
    }
    avg_path_size /= paths.size();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    std::cout << duration.count() << "," << min_path_size << "," << max_path_size << "," << "," << avg_path_size << "," << paths.size() << std::endl;
    return duration.count();
}

int main(int argc, char **argv) {
    OutputLogLevel = 2;

    HyperBasevector hb;
    vec<int> inv2;
    ReadPathVec paths2;
    int MAX_BPATHS = 100000;

    //Load hbv
    OutputLog(2) <<"Loading graph..." << std::endl;
    BinaryReader::readFile("step6.hbv", &hb);
    //Create inversion
    OutputLog(4) <<"Creating graph involution..." << std::endl;
    inv2.clear();
    hb.Involution(inv2);
    //load paths
    OutputLog(2) <<"Loading paths..." << std::endl;
    LoadReadPathVec(paths2,"step6.paths");
    //create path inversion
    OutputLog(2) << "Graph and paths loaded" << std::endl << std::endl;


    vec<vec<std::pair<int, int> > > xs;
    OutputLog(2) <<"Loading clusters..." << std::endl;
    BinaryReader::readFile("clusters.bin", &xs);
    //Create inversion
    OutputLog(4) <<"Clusters loaded..." << std::endl;


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
    OutputLog(2) << LR.size() << " unique clusters to be processed as blobs" << std::endl;
    OutputLog(2) << "laying out reads" << std::endl;
    // Some setup stuff.

    int nedges = hb.EdgeObjectCount();
    int K = hb.K();

    // Layout reads.  Expensive, temporary (?).

    std::vector<std::vector<int> > layout_pos(nedges);
    std::vector<std::vector<int64_t> > layout_id(nedges);
    std::vector<std::vector<bool>> layout_or(nedges);
    LayoutReads(hb, inv2, bases, paths2, layout_pos, layout_id, layout_or);

    // Make gap assemblies.

    vec<int> to_left, to_right;
    hb.ToLeft(to_left), hb.ToRight(to_right);
    vec<vec<basevector> > extras(LR.size());//this is accumulation, generates memory blocks
    vec<HyperBasevector> mhbp(LR.size());//this is accumulation, generates memory blocks
    //std::cout << Date() << ": processing " << LR.size() << " blobs" << std::endl;
    double clockp1 = WallClockTime();
    int nblobs = LR.size();
    uint64_t solved=0;
    //uint64_t solutionK[500];
    //for (auto &sk:solutionK) sk=0;

    //TODO: check local variable usage, should be made minimal!!!
    //Init readstacks, we'll need them!
    readstack::init_LUTs();
    OutputLog(2) << "processing blobs" << std::endl;

    auto bstart = std::max(0ul , bstart_arg);
    auto bend = std::min(xs.size(), bend_arg);

    for (uint64_t bl = bstart; bl < bend; ++bl) {
        //First part: create the gbases and gquals. this is locked by memory accesses and very convoluted
        const vec<int> &lefts = LR[bl].first, &rights = LR[bl].second; //TODO: how big is this? can we copy it?

        //Local readset
        std::vector<int64_t> pids;
        vecbasevector gbases;
        QualVecVec gquals;
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
            std::cout << "Found a solution for blob " << bl << std::endl;
            for (int i1 = 0; i1 < sources.isize(); i1++) {
                for (int i2 = 0; i2 < sinks.isize(); i2++) {
                    vec<vec<int>> p;
                    std::cout << "\tPATH," << sources[i1] << "," << sinks[i2] << ",";
                    time_edge_paths(xshb, zto_left, zto_right, sources[i1], sinks[i2], p, MAX_BPATHS);
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
                //++solutionK[xshb.K()];
            }
        }
    }

    return 0;
}