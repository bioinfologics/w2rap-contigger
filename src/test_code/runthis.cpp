//
// Created by Gonzalo Garcia (TGAC) on 28/10/2016.
//

#include <iostream>
#include <kmers/kmatch/KMatch.h>
#include "paths/long/large/ExtractReads.h"
//#include "pacbio/pacbio_pather.h"
//#include "pacbio/PathFinder_pb.h"

#include "TenX/TenX_pather.h"
#include "TenX/PathFinder_tx.h"


int main(int argc, char *argv[]){
  std::cout << "Loading reads..." << std::endl;
  InputDataMag dataMag(argv[1], argv[2]);
//  InputDataMag dataMag("/Users/ggarcia/Documents/arabidopsis_test_ds/configuration_file.config", "/Users/ggarcia/Documents/arabidopsis_test_ds/k_200");

  auto reads = dataMag.mag["TEX"]->rReads;
  std::cout << "Reads already loaded..." << std::endl;
  std::cout << reads.size() << " Reads in the vector" << std::endl;

  std::cout << "Loading hbv file..." << std::endl;
//  std::string fn = "/Users/ggarcia/Documents/arabidopsis_test_ds/k_200/athal_k200.large_K.clean.hbv";
  std::string fn = argv[3];

  // Load hbv and create inversion
  HyperBasevector hbv;
  BinaryReader::readFile(fn, &hbv);
  vec<int> inv;
  inv.clear();
  hbv.Involution(inv);



  // Create the paths and invert them
  TenXPather txp (&reads, &hbv);
  std::cout<< Date() << " Map creation." << std::endl;
  txp.createEmptyMap(&hbv);
  std::cout<< Date() << " Map creation done..." << std::endl;
  std::cout<< Date() << " Map filling with reads..." << std::endl;
  txp.reads2kmerTagMap();
  std::cout<< Date() << " Map filling with reads done..." << std::endl;

//  TenXPather* txp2 = &txp;
  // Pathfinder
  std::cout<< Date() << " Starting pathfinder..." << std::endl;
  PathFinder_tx pf_tx (&txp, &hbv, inv, 5);
  std::cout<< Date() << " done pathfinder..." << std::endl;
  pf_tx.untangle_complex_in_out_choices(1000, true);

//  // Print the map content
//  for (auto &t: txp.kmerTagMap){
//    if (t.second.size()>0) {
//      std::cout << "Key: " << t.first << ", Count: " << t.second.size() << std::endl;
//      for (auto tt: t.second){
//        std::cout << "--->Tag: " << tt.first << "-->" <<tt.second <<std::endl;
//      }
//    }
////  }
//  std::cout << "Size of the dictionary: " << txp.kmerTagMap.size() << std::endl;
//
//  // Intersection all vs all
//  std::cout<<Date()<<" Intersecting:" << std::endl;
//  auto edges = hbv.Edges();
//#pragma omp parallel for
//  for (auto i=0; i<edges.size(); ++i){
//    for (auto j=0; j<edges.size(); ++j) {
//      if (edges[i].size()>5000 & edges[j].size()>5000) {
//        auto interseccion = txp.edgeTagIntersection(edges[i].ToString(), edges[j].ToString(), 500);
//        if (interseccion.size() > 0) {
//#pragma omp critical (printest)
//          std::cout << "Print intersection: " << i << "-" << j << "->" << interseccion.size() << std::endl;
//        }
////      for (auto elemento: interseccion){
////        std::cout << elemento << std::endl;
////    }
//      }
//    }
//  }



//  auto tagidx = txp.kmerize_tag("ATCCACCGTGGTGCAA");
//  std::cout << "Tagkmer: " << tagidx << std::endl;
//  txp.Hbv2Map(&hbv);

//  auto histogram = txp.readsTagQc();
//  for (auto g=0; g<histogram.size(); ++g){
//    std::cout << g << "," << histogram[g]<<std::endl;
//  }

//  auto g = txp.getTagLinks();


//  PacbioPather pbp(&reads, &hbv);
//  pbp.Hbv2Map(&hbv);
//
//  ReadPathVec pathsr = pbp.mapReads();
//  VecULongVec invPaths;
//
//
//  invert(pathsr, invPaths, hbv.EdgeObjectCount());

  // pathfinders
//  PathFinder_tx(hbv, inv, pathsr, invPaths).unroll_loops(800);
//  PathFinder_tx(hbv,inv,pathsr,invPaths).untangle_pins();
//  PathFinder_tx(hbv,inv,pathsr,invPaths).untangle_complex_in_out_choices(700, true);

//  RemoveUnneededVertices2(hbv, inv, pathsr);
//  Cleanup(hbv, inv, pathsr);
//
//  BinaryWriter::writeFile("/Users/ggarcia/Documents/test_dataset/test_ecoli_pb/pf_after_loops.hbv", hbv);
//  WriteReadPathVec(pathsr, "/Users/ggarcia/Documents/test_dataset/test_ecoli_pb/pf_after_loops.paths");
}