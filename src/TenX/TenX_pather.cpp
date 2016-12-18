//
// Created by Gonzalo Garcia (TGAC) on 16/12/2016.
//

#include "TenX_pather.h"

TenXPather::TenXPather(std::vector<tenXRead>* aseqVector, HyperBasevector* ahbv): KMatch(31)  {
  //
  seqVector = aseqVector;
  hbv = ahbv;
}

int TenXPather::readsQc(){
  // Here calculate statistics regarding tags, indexes etcs
  std::map<std::string, int> index_map;
  return 0;
}

int TenXPather::makeIndexMap(bool to_disc=false){
  // Create a map with the index as key and a vector of indexes to the reads that have that index
  std::map<std::string, std::vector<int>> index_map;

  auto cont = 0;
  for (auto i=seqVector->begin(); i<seqVector->end(); ++i){
    auto idx = i->tag.ToString();
    index_map[idx].push_back(cont);
    cont++;
    }

  if (to_disc){
    for (auto const& i: index_map){
      std::cout << i.first << ":";
      for(auto const& v : i.second){
        std::cout << v << ",";
      }
      std::cout << std::endl;
    }
  }

  indexMap = index_map;
  return 0;
}

int TenXPather::makeTagMap(bool to_disc=false){
  // Create a map with the index as key and a vector of indexes to the reads that have that index
  std::map<std::string, std::vector<int>> tag_map;

  auto cont = 0;
  for (auto i=seqVector->begin(); i<seqVector->end(); ++i){
    auto idx = i->tag.ToString();
    tag_map[idx].push_back(cont);
    cont++;
  }

  if (to_disc){
    for (auto const& i: tag_map){
      std::cout << i.first << ":";
      for(auto const& v : i.second){
        std::cout << v << ",";
      }
      std::cout << std::endl;
    }
  }

  tagMap = tag_map;
  return 0;
}