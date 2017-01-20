//
// Created by Gonzalo Garcia (TGAC) on 19/01/2017.
//

#ifndef W2RAP_CONTIGGER_LOCALPATHER_H1
#define W2RAP_CONTIGGER_LOCALPATHER_H1

#include "kmers/kmatch/KMatch.h"
#include "paths/long/ReadPath.h"
#include "paths/long/large/ExtractReads.h"
//#include "TenX/TenX_pather.h"

class LocalPaths {
public:
    LocalPaths::LocalPaths(HyperBasevector* hbv, std::vector<std::vector<uint64_t>> pair_solutions, vec<int>& to_right, std::vector<BaseVec>& edges);

    int find_all_solution_paths();

    // Find all paths conecting 2 edges in the graph
    bool find_all_pair_conecting_paths(uint64_t edge_name , std::vector<uint64_t > path, int cont, uint64_t end_edge, int maxloop);


    virtual std::vector<uint64_t> choose_best_path(std::vector<std::vector<uint64_t>>* alternative_paths) = 0;

    std::vector<std::vector<uint64_t>> frontier_solutions;

    // 1st dim is the pair order, second dim one of the paths for that pair and 3rd dim ins the path itself
    std::vector<std::vector<uint64_t>> all_paths;

    std::vector<BaseVec>* mEdges;

    HyperBasevector* mHBV;

    std::vector<uint64_t> ins;
    std::vector<uint64_t> outs;


    std::vector<std::vector<uint64_t>> pair_temp_paths;
    // 1st dim is the pair order, second dim one of the paths for that pair and 3rd dim ins the path itself
//    std::vector<std::vector<std::vector<int>>> all_paths;

    vec<int>* mToLeft;
    vec<int>* mToRight;

};
#endif //W2RAP_CONTIGGER_LOCALPATHER_H1
