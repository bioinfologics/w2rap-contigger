//
// Created by Gonzalo Garcia (TGAC) on 28/10/2016.
//

#include <iostream>
#include <kmers/kmatch/KMatch.h>
#include "testcode_hbv.h"

int main(){
  std::string fn = "/Users/ggarcia/Documents/test_dataset/test/testrun.large_K.final.hbv";
//  hbv_explorer hbvE(fn);
  vecbvec reads;
  reads.ReadAll("/Users/ggarcia/Documents/test_dataset/test/PE1_frag_reads_orig.fastb");

  std::cout << reads.size() << "Reads in the vector" << std::endl;

  HyperBasevector hbv;
  BinaryReader::readFile(fn, &hbv);

  KMatch kmt(31);
  kmt.Hbv2Map(&hbv);
  std::cout << kmt.edgeMap.size() << std::endl;

  auto ed = kmt.edgeMap;

//  std::vector<uint64_t> v;
//  for(std::map<uint64_t, std::vector<std::pair<int, int>>>::iterator it = ed.begin(); it != ed.end(); ++it) {
//    v.push_back(it->first);
//    std::cout << it->first;
//    for (auto a: it->second){
//      std::cout << "\t" << a.first <<"-"<<a.second;
//    }
//    std::cout << std::endl;
//  }
//  std::cout << "Numero de keys: " << v.size() << std::endl;

  auto paths = kmt.MapReads(reads);

}