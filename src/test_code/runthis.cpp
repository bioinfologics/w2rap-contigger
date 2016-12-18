//
// Created by Gonzalo Garcia (TGAC) on 28/10/2016.
//

#include <iostream>
#include <kmers/kmatch/KMatch.h>
#include "paths/long/large/ExtractReads.h"
//#include "pacbio/pacbio_pather.h"
//#include "pacbio/PathFinder_pb.h"
#include "TenX/TenX_pather.h"

#include "testcode_hbv.h"


int main(){
  std::cout << "Loading reads..." << std::endl;
  InputDataMag dataMag("/Users/ggarcia/Documents/arabidopsis_test_ds/configuration_file.config", "/Users/ggarcia/Documents/arabidopsis_test_ds/k_200");
  auto reads = dataMag.mag["TEX"]->rReads;
  std::cout << "Reads already loaded..." << std::endl;
  std::cout << reads.size() << "Reads in the vector" << std::endl;

  std::cout << "Loading hbv file..." << std::endl;
  std::string fn = "/Users/ggarcia/Documents/arabidopsis_test_ds/k_200/athal_k200.large_K.clean.hbv";

  // Load hbv and create inversion
  HyperBasevector hbv;
  BinaryReader::readFile(fn, &hbv);
  vec<int> inv;
  inv.clear();
  hbv.Involution(inv);

  // Create the paths and invert them
  TenXPather txp(&reads, &hbv);
  txp.makeTagMap(true);
  txp.tagMap.size();

//  PacbioPather pbp(&reads, &hbv);
//  pbp.Hbv2Map(&hbv);
//
//  ReadPathVec pathsr = pbp.mapReads();
//  VecULongVec invPaths;
//
//
//  invert(pathsr, invPaths, hbv.EdgeObjectCount());

  // pathfinders
//  PathFinder_pb(hbv, inv, pathsr, invPaths).unroll_loops(800);
//  PathFinder_pb(hbv,inv,pathsr,invPaths).untangle_pins();
//  PathFinder_pb(hbv,inv,pathsr,invPaths).untangle_complex_in_out_choices(700, true);

//  RemoveUnneededVertices2(hbv, inv, pathsr);
//  Cleanup(hbv, inv, pathsr);
//
//  BinaryWriter::writeFile("/Users/ggarcia/Documents/test_dataset/test_ecoli_pb/pf_after_loops.hbv", hbv);
//  WriteReadPathVec(pathsr, "/Users/ggarcia/Documents/test_dataset/test_ecoli_pb/pf_after_loops.paths");
}