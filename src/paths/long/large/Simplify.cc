///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "CoreTools.h"
#include "Qualvector.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "paths/long/large/GapToyTools.h"
#include "paths/long/large/ImprovePath.h"
#include "paths/long/large/PullAparter.h"
#include "paths/long/large/Simplify.h"
// Standard bj pathfinder
#include "paths/PathFinder.h"
// Pabcio modified pathfinder pathfinder
#include "pacbio/PathFinder_pb.h"

void SimplifyWpb(const String &fin_dir, HyperBasevector &hb, vec<int> &inv,
              ReadPathVec &paths, const vecbasevector &bases, const VecPQVec &quals,
              const int MAX_SUPP_DEL, const Bool TAMP_EARLY, const int MIN_RATIO2,
              const int MAX_DEL2,
              const Bool ANALYZE_BRANCHES_VERBOSE2, const String &TRACE_SEQ,
              const Bool DEGLOOP, const Bool EXT_FINAL, const int EXT_FINAL_MODE,
              const Bool PULL_APART_VERBOSE, const vec<int> &PULL_APART_TRACE,
              const int DEGLOOP_MODE, const double DEGLOOP_MIN_DIST,
              const Bool IMPROVE_PATHS, const Bool IMPROVE_PATHS_LARGE,
              const Bool FINAL_TINY, const Bool UNWIND3, const bool RUN_PATHFINDER, const bool dump_pf_files, vecbvec& pb_bases) {
    // Improve read placements and delete funky pairs.
    std::cout << "Edge count: " << hb.EdgeObjectCount() << " Path count:" << paths.size() << std::endl;
    std::cout << "Simplify: rerouting paths" << std::endl;
    ReroutePaths(hb, inv, paths, bases, quals);
    DeleteFunkyPathPairs(hb, inv, bases, paths, False);

    // Remove unsupported edges in certain situations.
    std::cout << "Edge count: " << hb.EdgeObjectCount() << " Path count:" << paths.size() << std::endl;
    std::cout << "Simplify: removing unsupported edges" << std::endl;
    {
        const int min_mult = 10;
        vec<int> dels;
        {
            vec<int> support(hb.EdgeObjectCount(), 0);
            for (int64_t id = 0; id < (int64_t) paths.size(); id++) {
                for (int64_t j = 0; j < (int64_t) paths[id].size(); j++) {
                    int e = paths[id][j];
                    if (j >= 1) support[e]++;
                    if (inv[e] >= 0 && j < (int64_t) paths[id].size() - 1)
                        support[inv[e]]++;
                }
            }
#pragma omp parallel for
            for (int v = 0; v < hb.N(); v++) {
                if (hb.From(v).size() == 2) {
                    int e1 = hb.EdgeObjectIndexByIndexFrom(v, 0);
                    int e2 = hb.EdgeObjectIndexByIndexFrom(v, 1);
                    if (support[e1] > support[e2]) std::swap(e1, e2);
                    int s1 = support[e1], s2 = support[e2];
                    if (s1 <= MAX_SUPP_DEL && s2 >= min_mult * Max(1, s1)) {
#pragma omp critical
                        { dels.push_back(e1); }
                    }
                }
            }
        }
        {
            vec<int> support(hb.EdgeObjectCount(), 0);
            for (int64_t id = 0; id < (int64_t) paths.size(); id++) {
                for (int64_t j = 0; j < (int64_t) paths[id].size(); j++) {
                    int e = paths[id][j];
                    if (j < (int64_t) paths[id].size() - 1) support[e]++;
                    if (inv[e] >= 0 && j >= 1) support[inv[e]]++;
                }
            }
#pragma omp parallel for
            for (int v = 0; v < hb.N(); v++) {
                if (hb.To(v).size() == 2) {
                    int e1 = hb.EdgeObjectIndexByIndexTo(v, 0);
                    int e2 = hb.EdgeObjectIndexByIndexTo(v, 1);
                    if (support[e1] > support[e2]) std::swap(e1, e2);
                    int s1 = support[e1], s2 = support[e2];
                    if (s1 <= MAX_SUPP_DEL && s2 >= min_mult * Max(1, s1)) {
#pragma omp critical
                        { dels.push_back(e1); }
                    }
                }
            }
        }
        hb.DeleteEdges(dels);
        Cleanup(hb, inv, paths);
    }

    std::cout << "Edge count: " << hb.EdgeObjectCount() << " Path count:" << paths.size() << std::endl;
    std::cout << "Simplify: removing small Components" << std::endl;
    // Clean up assembly.

    RemoveSmallComponents3(hb);
    Cleanup(hb, inv, paths);

    if (TAMP_EARLY) {
        std::cout << "Edge count: " << hb.EdgeObjectCount() << " Path count:" << paths.size() << std::endl;
        std::cout << "Simplify: Tamping" << std::endl;
        Tamp(hb, inv, paths, 0);
    }

    RemoveHangs(hb, inv, paths, 100);
    Cleanup(hb, inv, paths);

    std::cout << "Edge count: " << hb.EdgeObjectCount() << " Path count:" << paths.size() << std::endl;
    std::cout << "Simplify: analysing branches" << std::endl;
    vec<int> to_right;
    hb.ToRight(to_right);

    AnalyzeBranches(hb, to_right, inv, paths, True, MIN_RATIO2, ANALYZE_BRANCHES_VERBOSE2);
    Cleanup(hb, inv, paths);
    RemoveHangs(hb, inv, paths, MAX_DEL2);
    Cleanup(hb, inv, paths);
    RemoveSmallComponents3(hb);
    Cleanup(hb, inv, paths);

    std::cout << "Edge count: " << hb.EdgeObjectCount() << " Path count:" << paths.size() << std::endl;
    std::cout << "Simplify: popping bubbles" << std::endl;
    PopBubbles(hb, inv, bases, quals, paths);
    Cleanup(hb, inv, paths);

    DeleteFunkyPathPairs(hb, inv, bases, paths, False);

    std::cout << "Edge count: " << hb.EdgeObjectCount() << " Path count:" << paths.size() << std::endl;
    std::cout << "Simplify: Tamping (700)" << std::endl;

    Tamp(hb, inv, paths, 10);
    RemoveHangs(hb, inv, paths, 700);
    Cleanup(hb, inv, paths);
    RemoveSmallComponents3(hb);
    Cleanup(hb, inv, paths);

    // Pull apart.

    {
        std::cout << Date() << ": making paths index for pull apart" << std::endl;
        VecULongVec invPaths;
        invert(paths, invPaths, hb.EdgeObjectCount());
        std::cout << Date() << ": pulling apart repeats" << std::endl;
        PullAparter pa(hb, inv, paths, invPaths, PULL_APART_TRACE, PULL_APART_VERBOSE, 5, 5.0);
        size_t count = pa.SeparateAll();
        std::cout << Date() << ": there were " << count << " repeats pulled apart." << std::endl;
        std::cout << Date() << ": there were " << pa.getRemovedReadPaths() << " read paths removed during separation."
                  << std::endl;
    }

    if (RUN_PATHFINDER) {
        std::cout << Date() << ": making paths index for PathFinder_pb" << std::endl;
        VecULongVec invPaths;
        invert(paths, invPaths, hb.EdgeObjectCount());
        if (dump_pf_files) {
            BinaryWriter::writeFile(fin_dir + "/pf_start.hbv", hb);
            //paths.WriteAll(fin_dir + "/pf_start.paths");
            WriteReadPathVec(paths,(fin_dir + "/pf_start.paths").c_str());
        }
        std::cout << Date() << ": PathFinder_pb: unrolling loops" << std::endl;
        PathFinder_pb(hb, inv, paths, invPaths).unroll_loops(800);
        std::cout << "Removing Unneded Vertices" << std::endl;
        RemoveUnneededVertices2(hb, inv, paths);
        Cleanup(hb, inv, paths);

        if (dump_pf_files) {
            BinaryWriter::writeFile(fin_dir + "/pf_unrolled_loops.hbv", hb);
            //paths.WriteAll(fin_dir + "/pf_unrolled_loops.paths");
            WriteReadPathVec(paths, (fin_dir + "/pf_unrolled_loops.paths").c_str());
        }
        
        // Create the paths just before the event
//        LongReadPather pbp(&pb_bases, &hb);
//        pbp.Hbv2Map(&hb);
//        auto pb_paths = pbp.mapReads();

        // Write inmediate before resolving paths to check
        std::string out_prefix = "graph";
        BinaryWriter::writeFile(fin_dir + "/" + out_prefix + ".before_pathfinder.hbv", hb);
        WriteReadPathVec(paths,(fin_dir + "/" + out_prefix + ".before_pathfinder.paths").c_str());
        
//        // dump paths to file to check the regions manually
//        std::ofstream pb_paths_file;
//        pb_paths_file.open(fin_dir + "/pacbiopaths.txt");
//        uint64_t paht_id = 0;

        auto edges = hb.Edges();
        std::cout << "edges loaded" << std::endl;
        vec<int> inversos;
        hb.Involution(inversos);

//        for (auto pbp: pb_paths){
//          pb_paths_file << paht_id;
//          for (auto eid: pbp){
//            pb_paths_file <<"," << eid << "("<<inversos[eid]<<")";
//          }
//          pb_paths_file << std::endl;
//          paht_id++;
//        }
//        pb_paths_file.close();

        for (auto i=3000;i < 3100 ; i+=500) {
            //std::cout << "-------------Parameter for pathfinder: " << i << "-----------" << std::endl
//            auto totalpaths=paths;
//            totalpaths.insert(totalpaths.end(),pb_paths.begin(),pb_paths.end());
            VecULongVec invtotalPaths;
            auto totalpaths=paths;
            invert(totalpaths, invtotalPaths, hb.EdgeObjectCount());

            std::cout << Date() << ": PathFinder_pb: resolving repeats of size " << i << std::endl;

            invtotalPaths.clear();
            invert(totalpaths, invtotalPaths, hb.EdgeObjectCount());
            std::cout << Date() << ": PathFinder_pb: analysing single-direction repeats" << std::endl;
            PathFinder_pb(hb, inv, totalpaths, invtotalPaths).untangle_complex_in_out_choices(i);
            std::cout << "Removing Unneded Vertices" << std::endl;
            RemoveUnneededVertices2(hb, inv, paths);
            Cleanup(hb, inv, paths);
        }

        if (dump_pf_files) {
            BinaryWriter::writeFile(fin_dir + "/pf_end.hbv", hb);
            //paths.WriteAll(fin_dir + "/pf_end.paths");
            WriteReadPathVec(paths,(fin_dir + "/pf_end.paths").c_str());
        }




    }
    // Improve paths.

    if (IMPROVE_PATHS) {
        path_improver pimp;
        vec<int64_t> ids;
        ImprovePaths(paths, hb, inv, bases, quals, ids, pimp,
                     IMPROVE_PATHS_LARGE, False);
    }

    // Extend paths.

    if (EXT_FINAL) {
        vec<int> to_left;
        hb.ToLeft(to_left), hb.ToRight(to_right);
        int ext = 0;
        auto qvItr = quals.begin();
        for (int64_t id = 0; id < (int64_t) paths.size(); id++, ++qvItr) {
            Bool verbose = False;
            const int min_gain = 20;
            ReadPath p = paths[id];
            ExtendPath2(paths[id], id, hb, to_left, to_right, bases[id], *qvItr,
                        min_gain, verbose, EXT_FINAL_MODE);
            if (p != paths[id]) ext++;
        }
        std::cout << ext << " paths extended" << std::endl;
    }

    // Degloop.

    if (DEGLOOP) {
        Degloop(DEGLOOP_MODE, hb, inv, paths, bases, quals, DEGLOOP_MIN_DIST);
        std::cout << Date() << ": removing Hangs" << std::endl;
        RemoveHangs(hb, inv, paths, 700);
        std::cout << Date() << ": cleanup" << std::endl;
        Cleanup(hb, inv, paths);
        std::cout << Date() << ": cleanup finished" << std::endl;

    }

    // Unwind three-edge plasmids.

    if (UNWIND3) UnwindThreeEdgePlasmids(hb, inv, paths);

    // Remove tiny stuff.

    if (FINAL_TINY) {
        std::cout << Date() << ": removing small components" << std::endl;
        RemoveSmallComponents3(hb, True);
        Cleanup(hb, inv, paths);
        CleanupLoops(hb, inv, paths);
        RemoveUnneededVerticesGeneralizedLoops(hb, inv, paths);
    }
}
