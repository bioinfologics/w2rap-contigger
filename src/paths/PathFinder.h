//
// Created by Bernardo Clavijo (TGAC) on 11/07/2016.
//

#ifndef W2RAP_CONTIGGER_PATHFINDER_H
#define W2RAP_CONTIGGER_PATHFINDER_H
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"

#include "paths/long/large/GapToyTools.h"

#include "Vec.h"


class PathFinder {
public:
    PathFinder( HyperBasevector& hbv, vec<int>& inv, ReadPathVec& paths, VecULongVec& invPaths, int min_reads = 5 ) :
    mHBV(hbv),
    mInv(inv),
    mPaths(paths),
    mEdgeToPathIds(invPaths),
    mMinReads(min_reads)
    {
        hbv.ToLeft(mToLeft);
        hbv.ToRight(mToRight);


    }
    void classify_forks();//how many forks of each type are there?
    std::vector<uint64_t> best_path_fw(uint64_t edge, int distance); //finds the best path forward for an edge
    std::array<uint64_t,3> transition_votes(uint64_t left_e,uint64_t right_e);
    std::array<uint64_t,3> path_votes(std::vector<uint64_t> path);
    std::array<uint64_t,3> multi_path_votes(std::vector<std::vector<uint64_t>> path);
    bool path_absolute_best(std::vector<uint64_t> path); //checks a path and its reverse, checks alternatives, true if shold be replaced
    void untangle_path(std::vector<uint64_t> path);
    void untangle_single_choices();//untangles all single choices when support is uncontested
    void untangle_complex_in_out_choices();
    void init_prev_next_vectors();
    uint64_t is_unrollable_loop(uint64_t e,uint64_t min_side_sizes);//returns size of the unrolled loop
    uint64_t paths_per_kbp(uint64_t e);
    std::string edge_pstr(uint64_t e);
    std::string path_str(std::vector<uint64_t> e);
    std::vector<uint64_t> get_all_long_frontiers(uint64_t e);



private:
    HyperBasevector& mHBV;
    vec<int>& mInv;
    ReadPathVec& mPaths;
    VecULongVec& mEdgeToPathIds;
    vec<int> mToLeft;
    vec<int> mToRight;
    std::vector<std::vector<uint64_t>> next_edges,prev_edges;
    int mMinReads;


};


#endif //W2RAP_CONTIGGER_PATHFINDER_H
