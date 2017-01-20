//
// Created by Gonzalo Garcia (TGAC) on 19/01/2017.
//

#include "LocalPather.h"
#include "paths/HyperBasevector.h"

LocalPaths::LocalPaths(HyperBasevector* hbv, std::vector<std::vector<uint64_t>> pair_solutions, vec<int> & to_right, std::vector<BaseVec>& edges){
  mHBV = hbv;
//  mTxp = txp;
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

