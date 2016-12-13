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
#include <kmers/kmatch/KMatch.h>
#include "CoreTools.h"
#include "Qualvector.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "paths/long/large/GapToyTools.h"
#include "paths/long/large/ImprovePath.h"
#include "paths/long/large/PullAparter.h"
#include "paths/long/large/Simplify.h"
#include "kmers/kmatch/KMatch.h"
#include "GFADump.h"
#include "OverlapValidator.h"


void graph_status(const HyperBasevector &hb) {
    uint64_t total=0;
    uint64_t null_sized=0;
    for (auto i=0; i<hb.EdgeObjectCount(); ++i) {
        if (hb.EdgeObject(i).size()>0)
        total+=hb.EdgeObject(i).size()-hb.K()+1;
        else ++null_sized;
    }
    std::cout << Date() << ": GRAPH contains " << total << " " <<hb.K()<<"-mers in " <<hb.EdgeObjectCount()<<" edges";
    if (null_sized>0) std::cout << " ("<<null_sized<<" gap edges)";
    std::cout << std::endl;
}
void path_status(const ReadPathVec &paths){
    uint64_t pe=0,ps=0,pm=0;
    for (auto &p:paths)
        if (p.size()==0) ++pe;
        else if (p.size()==1) ++ps;
        else ++pm;
    std::cout << Date() << ": PATHS -> empty: " << pe << "  single-edge: " << ps <<"  multi-edge: " << pm << std::endl;

}

void update_read_placements(HyperBasevector &hb, vec<int> &inv, ReadPathVec &paths, const vecbasevector &bases, const VecPQVec &quals){
    path_improver pimp;
    vec<int64_t> ids;
    ImprovePaths(paths, hb, inv, bases, quals, ids, pimp,
                 /*IMPROVE_PATHS_LARGE*/False, False);
    ReroutePaths(hb, inv, paths, bases, quals);
    DeleteFunkyPathPairs(hb, inv, bases, paths, False);

}

void update_read_placements_kmatch(HyperBasevector &hb, vec<int> &inv, ReadPathVec &paths, const vecbasevector &bases, const VecPQVec &quals) {
    std::cout<<Date()<<": creating kmatch object"<<std::endl;
    KMatch km(31);
    km.Hbv2Index(hb,31);
    //uint64_t mapped=0,c=0;
    /*km.Hbv2Map(hb);
    uint64_t TEST_READ_COUNT=100000;
    auto wc1=WallClockTime();
    std::cout<<Date()<<": Counting matches on "<<TEST_READ_COUNT<<" reads"<<std::endl;

    for (auto i=0;i<TEST_READ_COUNT;++i){
        if (km.countReadMatches(bases[i].ToString())>0) ++mapped;
    }
    std::cout<<Date()<<": "<<mapped<<" / "<<TEST_READ_COUNT<<" reads have hits to edges, "<<TimeSince(wc1)<<" spent in countReadMatches "<<std::endl;
    wc1=WallClockTime();
    std::cout<<Date()<<": Looking up "<<TEST_READ_COUNT<<" reads"<<std::endl;
    mapped=0;
    for (auto i=0;i<TEST_READ_COUNT;++i){
        if (km.lookupRead(bases[i].ToString()).size()>0) ++mapped;
    }
    std::cout<<Date()<<": "<<mapped<<" / "<<TEST_READ_COUNT<<" reads have hits to edges, "<<TimeSince(wc1)<<" spent in lookupRead "<<std::endl;
    wc1=WallClockTime();
    std::cout<<Date()<<": Looking up "<<TEST_READ_COUNT<<" reads in MAP"<<std::endl;
    mapped=0;
    for (auto i=0;i<TEST_READ_COUNT;++i){
        if (km.lookupReadInMap(bases[i].ToString()).size()>0) ++mapped;
    }
    std::cout<<Date()<<": "<<mapped<<" / "<<TEST_READ_COUNT<<" reads have hits to edges, "<<TimeSince(wc1)<<" spent in lookupReadInMap "<<std::endl;
    */
    std::atomic_uint_fast64_t mapped(0),mapped50(0),mappedUniquely(0),c(0);
    #pragma omp parallel for
    for (auto i=0;i<bases.size();++i){
        auto hits=km.lookupRead(bases[i].ToString());;
        if (hits.size()>0) ++mapped;
        if (hits.size()>(bases[i].size() - 31 + 1)*.5) {
            ++mapped50;
            bool unique=true;
            int64_t e=-1;
            for (auto & h:hits){
                if (e==-1) e=h.edge_id;
                if (e!=h.edge_id) {
                    unique =false;
                    break;
                }
            }
            if (unique) ++mappedUniquely;
        }

        ++c;
        if (c%100000==0) std::cout<<Date()<<": "<<mapped<<" / "<<c<<" reads have hits to edges"<<std::endl;
    }
    std::cout<<Date()<<": "<<mapped<<" / "<<bases.size()<<" reads have hits to edges"<<std::endl;
    std::cout<<Date()<<": "<<mapped50<<" / "<<bases.size()<<" reads have 50%+ kmers hitting to edges"<<std::endl;
    std::cout<<Date()<<": "<<mappedUniquely<<" / "<<bases.size()<<" reads have 50%+ kmers hitting same edge and no kmers hitting others"<<std::endl;
}

#include <iterator>

template<class ForwardIt, class T>
ForwardIt binary_find(ForwardIt first, ForwardIt last, const T& value)
{
    // Note: BOTH type T and the type after ForwardIt is dereferenced
    // must be implicitly convertible to BOTH Type1 and Type2, used in Compare.
    // This is stricter then lower_bound requirement (see above)

    first = std::lower_bound(first, last, value);
    return (first != last && value == *first) ? first : last;
}



void remove_unsupported_edges(HyperBasevector &hb, vec<int> &inv, ReadPathVec &paths, const vecbasevector &bases, const VecPQVec &quals, const int MAX_SUPP_DEL, const int min_mult){
    uint64_t delcount=1;
    uint64_t pass=0;
    while(delcount) {
        vec<int> toLeft,toRight;
        hb.ToLeft(toLeft);
        hb.ToRight(toRight);
        vec<int> dels;
        {
            std::vector<int> support_in(hb.EdgeObjectCount(), 0);
            std::vector<int> support_out(hb.EdgeObjectCount(), 0);
            //this is support of the edge as exit from its input vertex.

            for (int64_t id = 0; id < (int64_t) paths.size(); id++) {
                for (int64_t j = 0; j < (int64_t) paths[id].size(); j++) {
                    int e = paths[id][j];
                    if (j >= 1) {
                        support_in[e]++;
                        if (inv[e] >= 0) support_out[inv[e]]++;
                    }
                    if (j < (int64_t) paths[id].size() - 1) {
                        support_out[e]++;
                        if (inv[e] >= 0) support_in[inv[e]]++;
                    }
                }
            }

            //First (trivial) decisions: edges with no in and no out support
#pragma omp parallel for
            for (auto e=0;e<hb.EdgeObjectCount();e++){
                if (hb.EdgeObject(e).size()>1000) continue;
                if (support_in[e]<=MAX_SUPP_DEL && support_out[e]<=MAX_SUPP_DEL){
                    auto vfrom=toLeft[e];
                    auto vto=toRight[e];
                    bool alternative_from=hb.To(vfrom).size()==0;
                    bool alternative_to=hb.From(vto).size()==0;

                    for (auto fi=0;fi<hb.From(vfrom).isize();++fi) {
                        auto other = hb.EdgeObjectIndexByIndexFrom(vfrom, fi);
                        if (other != e and (inv[e] < 0 or other != inv[e]) and
                            support_in[other] >= min_mult * support_in[e] and
                            support_in[other] > min_mult and support_in[other] > MAX_SUPP_DEL) {
                            alternative_from = true;
                            break;
                        }
                    }

                    for (auto ti=0;ti<hb.To(vto).isize();++ti) {
                        auto other = hb.EdgeObjectIndexByIndexTo(vto, ti);
                        if (other != e and (inv[e] < 0 or other != inv[e]) and
                            support_out[other] >= min_mult * support_out[e] and
                            support_out[other] > min_mult and support_out[other] > MAX_SUPP_DEL) {
                            alternative_to = true;
                            break;
                        }
                    }
#pragma omp critical
                    if (alternative_from and alternative_to) {
                        dels.push_back(e);
                        if (inv[e]>0) dels.push_back(inv[e]);
                    }
                }
            }
        }
        auto before = hb.EdgeObjectCount();
        delcount = dels.size();
        std::sort(dels.begin(), dels.end());
        //Update paths first, just to be sure

        for (int64_t i = 0; i < (int64_t) paths.size(); i++) {
            for (int64_t j = 0; j < (int64_t) paths[i].size(); j++) {
                if (binary_find(dels.begin(), dels.end(), paths[i][j]) != dels.end()) {
                    paths[i].resize(j);
                    break;
                }
            }
        }

        GFADumpDetail("unsupported_paths_marked_detail"+std::to_string(pass),hb,inv,dels);
        hb.DeleteEdges(dels);
        Cleanup(hb, inv, paths);
        std::cout << Date() << ": " << delcount << " / " << before << " unsupported edges removed, "
                  << hb.EdgeObjectCount() << " edges after cleanup" << std::endl;
        // Improve read placements and delete funky pairs.
        std::cout << Date() << ": rerouting paths and cleaning pairs" << std::endl;
        ReroutePaths(hb, inv, paths, bases, quals);
        DeleteFunkyPathPairs(hb, inv, bases, paths, False);
        std::cout << Date() << ": improving paths" << std::endl;
        path_improver pimp;
        vec<int64_t> ids;
        ImprovePaths(paths, hb, inv, bases, quals, ids, pimp, false, False);
        pass++;
    }
}

void full_cleanup(HyperBasevector &hb, vec<int> &inv, ReadPathVec &paths, const int tampsize, const int hangssize){
    Tamp(hb, inv, paths, tampsize);
    Cleanup(hb, inv, paths);
    RemoveHangs(hb, inv, paths, 100);
    Cleanup(hb, inv, paths);
    RemoveSmallComponents3(hb);
    Cleanup(hb, inv, paths);
}



void Simplify(const String &fin_dir, HyperBasevector &hb, vec<int> &inv,
              ReadPathVec &paths, const vecbasevector &bases, const VecPQVec &quals,
              const int MAX_SUPP_DEL, const Bool TAMP_EARLY, const int MIN_RATIO2,
              const int MAX_DEL2,
              const Bool ANALYZE_BRANCHES_VERBOSE2, const String &TRACE_SEQ,
              const Bool DEGLOOP, const Bool EXT_FINAL, const int EXT_FINAL_MODE,
              const Bool PULL_APART_VERBOSE, const vec<int> &PULL_APART_TRACE,
              const int DEGLOOP_MODE, const double DEGLOOP_MIN_DIST,
              const Bool IMPROVE_PATHS, const Bool IMPROVE_PATHS_LARGE,
              const Bool FINAL_TINY, const Bool UNWIND3, const bool RUN_PATHFINDER, const bool dump_pf_files, const bool VERBOSE_PATHFINDER) {


    std::cout<<Date()<<": "<<(check_from_to(hb)? "graph adjacencies OK":"graph has incorrect vertex-edge adjacencies")<<std::endl;
    graph_status(hb);
    path_status(paths);

    // Improve read placements and delete funky pairs.
    std::cout << Date() << ": rerouting paths and cleaning pairs" << std::endl;
    ReroutePaths(hb, inv, paths, bases, quals);
    DeleteFunkyPathPairs(hb, inv, bases, paths, False);
    if (IMPROVE_PATHS) {
        std::cout << Date() << ": improving paths" << std::endl;

        path_improver pimp;
        vec<int64_t> ids;
        ImprovePaths(paths, hb, inv, bases, quals, ids, pimp, IMPROVE_PATHS_LARGE, False);
    }
    path_status(paths);
    OverlapValidator oval(hb,inv,paths);
    GFADumpDetail("before_ovlpval_detail",hb,inv);
    oval.compute_overlap_support();
    oval.analyse_complex_overlaps();
    std::vector<int> paint;
    paint.reserve(hb.EdgeObjectCount());
    for (int e:oval.find_perfect_tips(1000,5)) paint.push_back(e);
    std::cout<<Date()<<": "<<paint.size()<<" perfect tips found"<<std::endl;
    GFADumpDetail("ovlpval_perfect_tips_detail",hb,inv,paint);
    hb.DeleteEdges(paint);
    Cleanup(hb,inv,paths);
    graph_status(hb);
    path_status(paths);


    // Remove unsupported edges in certain situations.
    //const int min_mult=5;
    //std::cout << Date() << ": removing alternative edges with input support <="<<MAX_SUPP_DEL << std::endl;
    //remove_unsupported_edges(hb,inv,paths,bases,quals,MAX_SUPP_DEL,min_mult);




    // Clean up assembly.
    std::cout << Date() <<": removing small components"<<std::endl;
    RemoveSmallComponents3(hb);
    Cleanup(hb, inv, paths);


    /*
    if (TAMP_EARLY) {
        std::cout << Date() << ": Tamping" << std::endl;
        Tamp(hb, inv, paths, 0);
    }

    RemoveHangs(hb, inv, paths, 100);
    Cleanup(hb, inv, paths);*/
    graph_status(hb);
    path_status(paths);

    /*std::cout << Date() << ": analysing branches" << std::endl;


    AnalyzeBranches(hb, to_right, inv, paths, True, MIN_RATIO2, ANALYZE_BRANCHES_VERBOSE2);
    Cleanup(hb, inv, paths);
    RemoveHangs(hb, inv, paths, MAX_DEL2);
    Cleanup(hb, inv, paths);
    RemoveSmallComponents3(hb);
    Cleanup(hb, inv, paths);

    graph_status(hb);
    path_status(paths);*/

    /*std::cout << Date() << ": popping bubbles" << std::endl;
    PopBubbles(hb, inv, bases, quals, paths);

    Cleanup(hb, inv, paths);
    DeleteFunkyPathPairs(hb, inv, bases, paths, False);
    Tamp(hb, inv, paths, 10);
    RemoveHangs(hb, inv, paths, 700);
    Cleanup(hb, inv, paths);
    RemoveSmallComponents3(hb);
    Cleanup(hb, inv, paths);
    graph_status(hb);
    path_status(paths);*/

    // Pull apart.



    if (RUN_PATHFINDER) {
        //TODO: remove pull aparter once the pathfinder solves all repeats
        VecULongVec invPaths;
        std::cout << Date() << ": pulling apart canonical repeats" << std::endl;
        invert(paths, invPaths, hb.EdgeObjectCount());
        PullAparter pa(hb, inv, paths, invPaths, PULL_APART_TRACE, PULL_APART_VERBOSE, 5, 5.0);
        size_t count = pa.SeparateAll();
        std::cout << Date() << ": " << count << " repeats separated, " << pa.getRemovedReadPaths()
                  << " read paths removed" << std::endl;
        graph_status(hb);
        path_status(paths);



        std::cout << Date() << ": running pathfinder" << std::endl;
        invert(paths, invPaths, hb.EdgeObjectCount());
        if (dump_pf_files) {
            BinaryWriter::writeFile(fin_dir + "/pf_start.hbv", hb);
            WriteReadPathVec(paths,(fin_dir + "/pf_start.paths").c_str());
        }

        PathFinder pf(hb, inv, paths, invPaths,5,VERBOSE_PATHFINDER);
        pf.unroll_loops(800);
        RemoveUnneededVertices2(hb, inv, paths);
        Cleanup(hb, inv, paths);

        if (dump_pf_files) {
            BinaryWriter::writeFile(fin_dir + "/pf_unrolled_loops.hbv", hb);
            WriteReadPathVec(paths,(fin_dir + "/pf_unrolled_loops.paths").c_str());
        }
        invPaths.clear();
        invert( paths, invPaths, hb.EdgeObjectCount( ) );
        pf.untangle_complex_in_out_choices(700);
        RemoveUnneededVertices2(hb, inv, paths);
        Cleanup(hb, inv, paths);
        DeleteFunkyPathPairs(hb, inv, bases, paths, False);
        Tamp(hb, inv, paths, 10);
        RemoveHangs(hb, inv, paths, 700);
        Cleanup(hb, inv, paths);
        RemoveSmallComponents3(hb);
        graph_status(hb);
        path_status(paths);

        if (dump_pf_files) {
            BinaryWriter::writeFile(fin_dir + "/pf_end.hbv", hb);
            WriteReadPathVec(paths,(fin_dir + "/pf_end.paths").c_str());
        }
    } else {
        std::cout << Date() << ": pulling apart canonical repeats" << std::endl;
        VecULongVec invPaths;
        invert(paths, invPaths, hb.EdgeObjectCount());
        PullAparter pa(hb, inv, paths, invPaths, PULL_APART_TRACE, PULL_APART_VERBOSE, 5, 5.0);
        size_t count = pa.SeparateAll();
        std::cout << Date() << ": " << count << " repeats separated, " << pa.getRemovedReadPaths() << " read paths removed" << std::endl;
        graph_status(hb);
        path_status(paths);
    }
    // Improve paths.

    if (IMPROVE_PATHS) {
        path_improver pimp;
        vec<int64_t> ids;
        std::cout << Date() << ": improving paths" << std::endl;
        ImprovePaths(paths, hb, inv, bases, quals, ids, pimp,IMPROVE_PATHS_LARGE, False);
        path_status(paths);
    }

    // Extend paths.

    if (EXT_FINAL) {
        std::cout << Date() << ": extending paths" << std::endl;
        vec<int> to_left;
        vec<int> to_right;

        hb.ToLeft(to_left), hb.ToRight(to_right);
        int ext = 0;
        auto qvItr = quals.begin();
        for (int64_t id = 0; id < (int64_t) paths.size(); id++, ++qvItr) {
            Bool verbose = False;
            const int min_gain = 20;
            ReadPath p = paths[id];
            ExtendPath2(paths[id], id, hb, to_left, to_right, bases[id], *qvItr, min_gain, verbose, EXT_FINAL_MODE);
            if (p != paths[id]) ext++;
        }
        std::cout << Date() << ": "<< ext << " paths extended" << std::endl;
        path_status(paths);
    }

    // Degloop.

    if (DEGLOOP) {
        std::cout << Date() << ": deglooping" << std::endl;
        Degloop(DEGLOOP_MODE, hb, inv, paths, bases, quals, DEGLOOP_MIN_DIST);
        RemoveHangs(hb, inv, paths, 700);
        Cleanup(hb, inv, paths);
        graph_status(hb);
        path_status(paths);
    }

    // Unwind three-edge plasmids.

    if (UNWIND3) {
        std::cout << Date() << ": unwinding 3-edge plasmids" << std::endl;
        UnwindThreeEdgePlasmids(hb, inv, paths);
        graph_status(hb);
        path_status(paths);
    }

    // Remove tiny stuff.

    if (FINAL_TINY) {
        std::cout << Date() << ": removing small components" << std::endl;
        RemoveSmallComponents3(hb, True);
        Cleanup(hb, inv, paths);
        CleanupLoops(hb, inv, paths);
        RemoveUnneededVerticesGeneralizedLoops(hb, inv, paths);
        graph_status(hb);
        path_status(paths);
    }
}