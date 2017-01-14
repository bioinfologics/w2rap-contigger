//
// Created by Bernardo Clavijo (TGAC) on 11/07/2016.
//

#ifndef W2RAP_CONTIGGER_PATHFINDER_H_TX
#define W2RAP_CONTIGGER_PATHFINDER_H_TX

#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "paths/long/large/GapToyTools.h"
#include "Vec.h"
#include "TenX/TenX_pather.h"

class PathFinder_tx {
public:
    PathFinder_tx(TenXPather* txp, HyperBasevector* hbv, vec<int> inv, int min_reads = 5 );

    //Graph-related methods
    std::vector<std::vector<uint64_t>> AllPathsFromTo(std::vector<uint64_t> in_edges, std::vector<uint64_t> out_edges, uint64_t max_length);


    //ReadPath-related methods


    void classify_forks();//how many forks of each type are there?
    std::vector<uint64_t> best_path_fw(uint64_t edge, int distance); //finds the best path forward for an edge
    std::array<uint64_t,3> transition_votes(uint64_t left_e,uint64_t right_e);
    std::array<uint64_t,3> path_votes(std::vector<uint64_t> path);
    std::array<uint64_t,3> multi_path_votes(std::vector<std::vector<uint64_t>> path);
    bool path_absolute_best(std::vector<uint64_t> path); //checks a path and its reverse, checks alternatives, true if shold be replaced
    void untangle_path(std::vector<uint64_t> path);
    void untangle_pins();
    void unroll_loops(uint64_t min_side_sizes);//untangles all single choices when support is uncontested
    void untangle_complex_in_out_choices(uint64_t large_frontier_size, bool verbose_separation=false);
    void init_prev_next_vectors();
    std::vector<std::vector<uint64_t>> is_unrollable_loop(uint64_t e,uint64_t min_side_sizes);//returns size of the unrolled loop
    uint64_t paths_per_kbp(uint64_t e);
    std::string edge_pstr(uint64_t e);
    std::string path_str(std::vector<uint64_t> e);
    std::map<uint64_t,std::vector<uint64_t>> separate_path(std::vector<uint64_t> p, bool verbose_separation=false);
    bool join_edges_in_path(std::vector<uint64_t> p);
    std::array<std::vector<uint64_t>,2>  get_all_long_frontiers(uint64_t e,uint64_t large_frontier_size);
    void migrate_readpaths(std::map<uint64_t,std::vector<uint64_t>> edgemap);


private:
    TenXPather* mTxp;
    HyperBasevector* mHBV;
    vec<int> mInv;
    ReadPathVec& mPaths;
    VecULongVec& mEdgeToPathIds;
    vec<int> mToLeft;
    vec<int> mToRight;
    std::vector<std::vector<uint64_t>> next_edges,prev_edges;
    int mMinReads;

};

class LocalPaths {
public:
    LocalPaths::LocalPaths(HyperBasevector* hbv, std::vector<std::vector<int>> pair_solutions, vec<int>& to_right, TenXPather* txp, std::vector<BaseVec>& edges);

    int find_all_solution_paths();

    // Find all paths conecting 2 edges in the graph
    bool find_all_pair_conecting_paths(int edge_name , std::vector<int> path, int cont, int end_edge, int maxloop);

    std::vector<int> choose_best_path(std::vector<std::vector<int>>* alternative_paths);

    std::vector<std::vector<int>> frontier_solutions;

    // 1st dim is the pair order, second dim one of the paths for that pair and 3rd dim ins the path itself
    std::vector<std::vector<int>> all_paths;

private:
    TenXPather* mTxp;

    std::vector<BaseVec>* mEdges;

    HyperBasevector* mHBV;

//    std::vector<std::vector<int>> frontier_solutions;
    std::vector<int> ins;
    std::vector<int> outs;


    std::vector<std::vector<int>> pair_temp_paths;
    // 1st dim is the pair order, second dim one of the paths for that pair and 3rd dim ins the path itself
//    std::vector<std::vector<std::vector<int>>> all_paths;

    vec<int>* mToLeft;
    vec<int>* mToRight;

};
#endif //W2RAP_CONTIGGER_PATHFINDER_H_TX
