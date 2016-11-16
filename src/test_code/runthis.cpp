//
// Created by Gonzalo Garcia (TGAC) on 28/10/2016.
//

#include <iostream>
#include <kmers/kmatch/KMatch.h>
#include "paths/long/large/ExtractReads.h"
#include "test_code/pacbio/pacbio_pather.h"

#include "testcode_hbv.h"


int main(){
  std::cout << "Loading reads..." << std::endl;
  InputDataMag dataMag("/Users/ggarcia/Documents/test_dataset/configuration_file.config", "/Users/ggarcia/Documents/test_dataset/test");
  vecbvec reads = dataMag.mag["PB1"]->bases;
  std::cout << "Reads already loaded..." << std::endl;
  std::cout << reads.size() << "Reads in the vector" << std::endl;

  std::cout << "Loading hbv file..." << std::endl;
  std::string fn = "/Users/ggarcia/Documents/ecoli_test/ecolitest/ecoli_k200.contig.hbv";
//  hbv_explorer hbvE(fn);
  HyperBasevector hbv;
  BinaryReader::readFile(fn, &hbv);

//  KMatch kmt(31);
//  kmt.Hbv2Map(&hbv);
//  std::cout << kmt.edgeMap.size() << std::endl;

//  auto ed = kmt.edgeMap;
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

//  auto paths = kmt.MapReads(reads, &hbv);
  PacbioPather pbp(&reads, &hbv);
  pbp.Hbv2Map(&hbv);
  auto pb_paths = pbp.mapReads();

}