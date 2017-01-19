//
// Created by Gonzalo Garcia (TGAC) on 19/01/2017.
//

#include "LocalPather.h"
#include "paths/HyperBasevector.h"

LocalPaths::LocalPaths(HyperBasevector* hbv, std::vector<std::vector<uint64_t>> pair_solutions, vec<int> & to_right, TenXPather* txp, std::vector<BaseVec>& edges){
  mHBV = hbv;
  mTxp = txp;
  frontier_solutions = pair_solutions;
  mToRight = &to_right;
  mEdges = &edges;

  for (auto p: pair_solutions){
    ins.push_back(p[0]);
    outs.push_back(p[1]);
  }
}

bool LocalPaths::find_all_pair_conecting_paths(uint64_t edge_name , std::vector<uint64_t> path, int cont, uint64_t end_edge, int maxloop = 0) {
  // Find all paths between 2 given edges

  path.push_back(edge_name);
  cont ++;
  auto r_vertex = (*mToRight)[edge_name];

  // if the node is the end_edge or the run reaches a frontier limit the recursion completes !
  if (edge_name == end_edge){
    pair_temp_paths.push_back(path);

    return true;

  } else if (find(outs.begin(), outs.end(), edge_name) != outs.end()) {
    return false;
  }

  for (auto en=0; en<mHBV->FromSize(r_vertex); ++en){
    edge_name =mHBV->EdgeObjectIndexByIndexFrom(r_vertex, en);
    if (std::count(path.begin(), path.end(), edge_name) < maxloop+1){
      std::vector<uint64_t> subpath (path.begin(), path.begin()+cont);
      find_all_pair_conecting_paths(edge_name , subpath, cont, end_edge);
    }
  }
}

int LocalPaths::find_all_solution_paths(){
  // Get each one of the pairs for the frontier solutions finded and fill the edges and return every possible path
  // between the ends of each pair

  // For each of the pairs in the solved frontier pairs
  for (auto pn=0; pn < frontier_solutions.size(); ++pn){
    auto p_in = frontier_solutions[pn][0];
    auto p_out = frontier_solutions[pn][1];

    // Clear the temp pairs for this pair
    pair_temp_paths.clear();

    // Find all paths between the ends of the path
    std::vector<uint64_t> path;
    int cont = 0;
    auto status = find_all_pair_conecting_paths(p_in, path, cont, p_out);

    // Vote the best path from the all available
    auto best_path = choose_best_path(&pair_temp_paths);
    if (best_path.size()>0){
      all_paths.push_back(best_path);
    }
    else {
      std::cout << "Skipping an empty path between :" << p_in <<" and "<<  p_out << std::endl;
    }
  }
}

std::vector<uint64_t> LocalPaths::choose_best_path(std::vector<std::vector<uint64_t>>* alternative_paths){
  // Choose the best path for the combination from the pairs list

  // If there is only one posible path return that path and finish
  if (alternative_paths->size() == 0){
    std::cout <<"This path still CEROOOO" <<std::endl;
    std::vector<uint64_t> nopaths;
    return nopaths;
  }
  if (alternative_paths->size() == 1){
    std::cout << "Only one patha available: " << std::endl;
    return (*alternative_paths)[0];
  } else {
    // If there is more than one path vote for the best (the criteria here is most tag density (presentTags/totalKmers)
    float best_path = 0.0;
    float best_path_score = 0.0;
    for (auto path_index = 0; path_index < alternative_paths->size(); ++path_index) {
      float cpath_score = 0;
      for (auto ei = 0; ei < (*alternative_paths)[path_index].size() - 1; ++ei) {
        auto from_edge_string = (*mEdges)[(*alternative_paths)[path_index][ei]].ToString();
        auto to_edge_string = (*mEdges)[(*alternative_paths)[path_index][ei + 1]].ToString();
        cpath_score += mTxp->edgeTagIntersection(from_edge_string, to_edge_string, 1500);
      }
      if (cpath_score > best_path_score) {
        best_path = path_index;
        best_path_score = cpath_score;
      }
    }
    std::cout << "Best path selected: " << best_path << ", score: " << best_path_score << std::endl;
    for (auto p=0; p<(*alternative_paths)[best_path].size(); ++p){
      std::cout << (*alternative_paths)[best_path][p] << ",";
    }
    std::cout << std::endl;

    return (*alternative_paths)[best_path];
  }
}