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
    PathFinder( HyperBasevector& hbv, vec<int>& inv, ReadPathVec& paths, VecULongVec& invPaths, const vecbasevector& bases, const VecPQVec& quals, int min_reads = 5, bool verbose=false ) :
    mHBV(hbv),
    mInv(inv),
    mPaths(paths),
    mEdgeToPathIds(invPaths),
    mMinReads(min_reads),
    mVerbose(verbose),
    mBases(bases),
    mQuals(quals)
    {


    }

    // Update and info-collecting functions
    void update_prev_next();
    void improve_paths();
    /**
     * Creates a single vector with the path from r1id and r1id+1, pe jump is specified with -1
     * @param r1id
     * @param collapse_overlaps WARNING: this can "false bridge" across small loops.
     * @return
     */
    std::vector<int64_t> get_full_paired_path(int r1id, bool collapse_overlaps=true);
    /**
     * Collects and counts all paths larger than 1 (including pairs, -1 are jumps) for this edge or its inversion
     * Needs to check read position to account for loop-stile linkage to itself!!!
     * @param e
     * @return vector of (count,path) pairs
     */
    std::vector<std::pair<uint64_t,std::vector<int64_t>>> collect_paths_for_edge(int e, bool collect_inverse=true, bool collapse_overlaps=true);

    // Current simplification heuristics
    void clip_tips(int max_length, int max_reads);
    void solve_perfect_repeats(int max_size);
    void pop_bubbles(int max_size);
    void unroll_loops(int repeat_size, int loop_size, int side_size);
    void extend_until_repeat(int max_collapsed_size);

    // Graph-related methods
    std::vector<std::vector<uint64_t>> AllPathsFromTo(std::vector<uint64_t> in_edges, std::vector<uint64_t> out_edges, uint64_t max_length);
    void separate_solutions(std::vector<std::vector<uint64_t >> sol_paths);

    //ReadPath-related methods



    std::array<uint64_t,3> transition_votes(uint64_t left_e,uint64_t right_e);
    std::array<uint64_t,3> path_votes(std::vector<uint64_t> path);
    std::array<uint64_t,3> multi_path_votes(std::vector<std::vector<uint64_t>> path);
    void untangle_pins();
    void unroll_loops(uint64_t min_side_sizes);//untangles all single choices when support is uncontested
    void untangle_complex_in_out_choices(uint64_t large_frontier_size, bool verbose_separation=false);

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
    const vecbasevector& mBases;
    const VecPQVec& mQuals;
    vec<int> mToLeft;
    vec<int> mToRight;
    std::vector<std::vector<uint64_t>> next_edges,prev_edges;
    int mMinReads;
    bool mVerbose;


};

void simplifyWithPathFinder( HyperBasevector& hbv, vec<int>& inv, ReadPathVec& paths, VecULongVec& invPaths, const vecbasevector& bases, const VecPQVec& quals, int min_reads = 5, bool verbose=false, bool dump_intermediate_gfas = false );
#endif //W2RAP_CONTIGGER_PATHFINDER_H
