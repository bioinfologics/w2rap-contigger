//
// Created by Bernardo Clavijo (TGAC) on 11/07/2016.
//

#include "TenX/PathFinder_tx.h"

PathFinder_tx::PathFinder_tx (TenXPather& txp, HyperBasevector& hbv, vec<int> inv, int min_reads, ReadPathVec& paths, VecULongVec& invPaths)
    : PathFinder(hbv, inv, paths, invPaths, 5), mTxp (txp) {
}

void PathFinder_tx::solve_region_using_TenX(uint64_t large_frontier_size, bool verbose_separation, float score_threshold) {
    //find a complex path
    // [GOnza] TODO: separate this function into 2 differets, one to create the map an the other to test the permutations
    auto edges = mHBV.Edges();

    uint64_t qsf=0,qsf_paths=0;
    uint64_t msf=0,msf_paths=0;
    init_prev_next_vectors();
    std::cout<<"vectors initialised"<<std::endl;
    std::set<std::array<std::vector<uint64_t>,2>> seen_frontiers,solved_frontiers;
    std::vector<std::vector<uint64_t>> paths_to_separate;
    int solved_regions = 0;
    int unsolved_regions = 0;
    for (int e = 0; e < mHBV.EdgeObjectCount(); ++e) {
        if (e < mInv[e] && mHBV.EdgeObject(e).size() < large_frontier_size) {
            // Get the frontiers fo the edge
            // [GONZA] TODO: check the return details to document
            auto f=get_all_long_frontiers(e, large_frontier_size);
            if (f[0].size()>1 and f[1].size()>1 and f[0].size() == f[1].size() and seen_frontiers.count(f)==0){
                seen_frontiers.insert(f);
                bool single_dir=true;
                std::map<std::string, float> shared_paths;
                for (auto in_e:f[0]) for (auto out_e:f[1]) if (in_e==out_e) {single_dir=false;break;}

                if (single_dir) {
                    // If there is a region to resolve within the frontiers
                    std::cout<<" Single direction frontiers for complex region on edge "<<e<<" IN:"<<path_str(f[0])<<" OUT: "<<path_str(f[1])<<std::endl;
                    std::vector<int> in_used(f[0].size(),0);
                    std::vector<int> out_used(f[1].size(),0);
                    std::vector<std::vector<uint64_t>> first_full_paths;
                    bool reversed=false;

                    // intersect all pairs of ins and outs to score, save the scores in a map to score the combinations later
                    for (auto in_i=0;in_i<f[0].size();++in_i) {
                        auto in_e=f[0][in_i];
                        for (auto out_i=0;out_i<f[1].size();++out_i) {
                            auto out_e=f[1][out_i];

                            int edges_in_path;
                            std::string pid;
                            auto in_e_seq = edges[in_e].ToString();
                            auto out_e_seq = edges[out_e].ToString();

                            // Intersect the tags for the edges
                            auto intersection_score = mTxp.edgeTagIntersection(in_e_seq, out_e_seq, 1500);
                            if (intersection_score>score_threshold){
                                // if the edges overlap in the tagspace thay are added to the map and the combination is markes in the used edges
                                pid = std::to_string(in_e) + "-" + std::to_string(out_e);
                                shared_paths[pid] += intersection_score; // This should score the link based in the number of tags that tha pair shares
                                out_used[out_i]++;
                                in_used[in_i]++;
                            }
                        }
                    } // When this is done i get a map of all combinations of ins and outs to score the permutations in the next step

                    // Here all combinations are counted, now i need to get the best configuration between nodes
                    auto in_frontiers = f[0];
                    auto out_frontiers = f[1];
                    float max_score = -9999.0;
                    int perm_number = 0;

                    std::vector<int> max_score_permutation;
                    do {
                        // Vectors to count seen edges (check that all edges are included in th solution)
                        std::vector<int> seen_in(in_frontiers.size(), 0);
                        std::vector<int> seen_out(in_frontiers.size(), 0);

                        //std::cout << "-----------------------------Testing permutaiton ------------------" << std::endl;
                        float current_score = 0;
                        for (auto pi=0; pi<in_frontiers.size(); ++pi){
                            std::string index = std::to_string(in_frontiers[pi])+"-"+std::to_string(out_frontiers[pi]);
                            if ( shared_paths.find(index) != shared_paths.end() ){
                                // Mark the pair as seen in this iteration
                                seen_in[pi]++;
                                seen_out[pi]++;
                                current_score += shared_paths[index];
                            }
                        }
                        // Check that all boundaries are used in the permutation
                        bool all_used = true;
                        for (auto a=0; a<seen_in.size(); ++a){
                            if (0==seen_in[a] or 0==seen_out[a]){
                                all_used = false;
                            }
                        }

                        if (current_score>max_score and all_used){
                            max_score = current_score;
                            max_score_permutation.clear();
                            for (auto aix: out_frontiers){
                                max_score_permutation.push_back(aix);
                            }
                        }
                        perm_number++;
                    } while (std::next_permutation(out_frontiers.begin(), out_frontiers.end()));

                    // Get the solution
                    if (max_score>score_threshold){
                        std::cout << " Found solution to region: " <<std::endl;
                        solved_regions++;

                        std::vector<std::vector<uint64_t>> wining_permutation;
                        for (auto ri=0; ri<max_score_permutation.size(); ++ri){
                            std::cout << in_frontiers[ri] << "(" << mInv[in_frontiers[ri]] << ") --> " << max_score_permutation[ri] << "("<< mInv[max_score_permutation[ri]] <<"), Score: "<<  max_score << std::endl;
                            std::vector<uint64_t> tp = {in_frontiers[ri], max_score_permutation[ri]};
                            wining_permutation.push_back(tp);
                        }
                        std::cout << "--------------------" << std::endl;

                        // Fill intermediate nodes
                        LocalPaths_TX lp (mHBV, wining_permutation, mToRight, mTxp, edges);
                        lp.find_all_solution_paths();
                        std::cout << "All paths done" << std::endl;
                        for (auto spi = 0; spi<lp.all_paths.size(); ++spi){
                            auto sv = lp.all_paths[spi];
                            paths_to_separate.push_back(sv);
                        }
                        std::cout << "--------------------" << std::endl;
                    } else {
//                        std::cout << "Region not resolved, not enough links or bad combinations" << max_score << std::endl;
                        unsolved_regions++;
                    }

                }
            }
        }
    }
    std::cout << "========================" << std::endl;
    std::cout << "Solved regions: " << solved_regions << ", Unsolved regions: " << unsolved_regions << " --> " << (float)solved_regions / (float)(solved_regions + unsolved_regions) * 100 << " %" << std::endl;
    std::cout << "========================" << std::endl;

    std::cout << "List of paths to separate" << std::endl;
    for (auto ss: paths_to_separate){
      for (auto sx: ss){
        std::cout << sx << ",";
      }
      std::cout << std::endl;
    }
    std::cout << "Paths to separate: " << paths_to_separate.size() << std::endl;

    uint64_t sep=0;
    std::map<uint64_t,std::vector<uint64_t>> old_edges_to_new;
    for (auto p: paths_to_separate){
        if (old_edges_to_new.count(p.front()) > 0 or old_edges_to_new.count(p.back()) > 0) {
            std::cout<<"WARNING: path starts or ends in an already modified edge, skipping"<<std::endl;
            continue;
        }

        auto oen=separate_path(p, verbose_separation);
        std::cout << "End separate paths, eon size: " << oen.size() << std::endl;
        if (oen.size()>0) {
            for (auto et:oen){
                if (old_edges_to_new.count(et.first)==0) old_edges_to_new[et.first]={};
                for (auto ne:et.second) old_edges_to_new[et.first].push_back(ne);
            }
            sep++;
            std::cout << "separated_paths: " << sep << std::endl;
        }
    }
//    if (old_edges_to_new.size()>0) {
//        migrate_readpaths(old_edges_to_new);
//    }
    std::cout<<" "<<sep<<" paths separated!"<<std::endl;
}