//
// Created by Gonzalo Garcia (TGAC) on 28/10/2016.
//

#include <iostream>
#include <kmers/kmatch/KMatch.h>
#include "paths/long/large/ExtractReads.h"
#include "test_code/pacbio/pacbio_pather.h"
#include "paths/PathFinder.h"

#include "testcode_hbv.h"


int main(){
  std::cout << "Loading reads..." << std::endl;
  InputDataMag dataMag("/Users/ggarcia/Documents/test_dataset/configuration_file.config", "/Users/ggarcia/Documents/test_dataset/test");
  vecbvec reads = dataMag.mag["PB1"]->bases;
  std::cout << "Reads already loaded..." << std::endl;
  std::cout << reads.size() << "Reads in the vector" << std::endl;

  std::cout << "Loading hbv file..." << std::endl;
  std::string fn = "/Users/ggarcia/Documents/ecoli_test/ecolitest/ecoli_k200.contig.hbv";

  // Load hbv and create inversion
  HyperBasevector hbv;
  BinaryReader::readFile(fn, &hbv);
  vec<int> inv;
  inv.clear();
  hbv.Involution(inv);

  // Create the paths and invert them
  PacbioPather pbp(&reads, &hbv);
  pbp.Hbv2Map(&hbv);

  ReadPathVec pathsr = pbp.mapReads();
  VecULongVec invPaths;


  invert(pathsr, invPaths, hbv.EdgeObjectCount());

  // pathfinders
//  PathFinder_pb(hbv, inv, pathsr, invPaths).unroll_loops(800);
//  PathFinder_pb(hbv,inv,pathsr,invPaths).untangle_pins();
  PathFinder_pb(hbv,inv,pathsr,invPaths).untangle_complex_in_out_choices(700, true);

  RemoveUnneededVertices2(hbv, inv, pathsr);
  Cleanup(hbv, inv, pathsr);

  BinaryWriter::writeFile("/Users/ggarcia/Documents/test_dataset/test_ecoli_pb/pf_after_loops.hbv", hbv);
  WriteReadPathVec(pathsr, "/Users/ggarcia/Documents/test_dataset/test_ecoli_pb/pf_after_loops.paths");
}