//
// Created by Bernardo Clavijo (TGAC) on 11/07/2016.
//

#ifndef W2RAP_CONTIGGER_PATHFINDER_H
#define W2RAP_CONTIGGER_PATHFINDER_H
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "paths/long/large/GapToyTools.h"


class PathFinder {
public:
    PathFinder( HyperBasevector& hbv, vec<int>& inv, ReadPathVec& paths, VecULongVec& invPaths, int min_reads = 5, bool verbose=false ) :
    mHBV(hbv),
    mInv(inv),
    mPaths(paths),
    mEdgeToPathIds(invPaths),
    mMinReads(min_reads),
    mVerbose(verbose)
    {


    }
    //Simplification routines
    void clip_tips(int max_length, int max_reads);
    void extend_until_repeat();
    //Graph-related methods
    std::vector<std::vector<uint64_t>> AllPathsFromTo(std::vector<uint64_t> in_edges, std::vector<uint64_t> out_edges, uint64_t max_length);


    //ReadPath-related methods

    std::array<uint64_t,3> transition_votes(uint64_t left_e,uint64_t right_e);
    std::array<uint64_t,3> path_votes(std::vector<uint64_t> path);
    std::array<uint64_t,3> multi_path_votes(std::vector<std::vector<uint64_t>> path);
    void untangle_pins();
    void unroll_loops(uint64_t min_side_sizes);//untangles all single choices when support is uncontested
    void untangle_complex_in_out_choices(uint64_t large_frontier_size, bool verbose_separation=false);
    void update_prev_next();
    std::vector<std::vector<uint64_t>> is_unrollable_loop(uint64_t e,uint64_t min_side_sizes);//returns size of the unrolled loop
    uint64_t paths_per_kbp(uint64_t e);
    std::string edge_pstr(uint64_t e);
    std::string path_str(std::vector<uint64_t> e);
    std::map<uint64_t,std::vector<uint64_t>> separate_path(std::vector<uint64_t> p, bool verbose_separation=false);
    bool join_edges_in_path(std::vector<uint64_t> p);
    std::array<std::vector<uint64_t>,2>  get_all_long_frontiers(uint64_t e,uint64_t large_frontier_size);
    void migrate_readpaths(std::map<uint64_t,std::vector<uint64_t>> edgemap);
    void extend_bridging_paths();



private:
    HyperBasevector& mHBV;
    vec<int>& mInv;
    ReadPathVec& mPaths;
    VecULongVec& mEdgeToPathIds;
    vec<int> mToLeft;
    vec<int> mToRight;
    std::vector<std::vector<uint64_t>> next_edges,prev_edges;
    int mMinReads;
    bool mVerbose;


};

void simplifyWithPathFinder( HyperBasevector& hbv, vec<int>& inv, ReadPathVec& paths, VecULongVec& invPaths, int min_reads = 5, bool verbose=false, bool dump_intermediate_gfas = false );
#endif //W2RAP_CONTIGGER_PATHFINDER_H
