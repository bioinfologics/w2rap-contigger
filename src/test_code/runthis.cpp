//
// Created by Gonzalo Garcia (TGAC) on 28/10/2016.
//

#include <iostream>
#include <kmers/kmatch/KMatch.h>
#include "paths/long/large/ExtractReads.h"

#include "pacbio/pacbio_pather.h"

//#include "TenX/TenX_pather.h"
//#include "TenX/PathFinder_tx.h"


int main(int argc, char *argv[]){
  std::cout << "Loading reads..." << std::endl;
  InputDataMag dataMag(argv[1], argv[2]);
//  InputDataMag dataMag("/Users/ggarcia/Documents/arabidopsis_test_ds/configuration_file.config", "/Users/ggarcia/Documents/arabidopsis_test_ds/k_200");



  std::cout << "Loading hbv file..." << std::endl;
//  std::string fn = "/Users/ggarcia/Documents/arabidopsis_test_ds/k_200/athal_k200.large_K.clean.hbv";
  std::string fn = argv[3];

  // Load hbv and create inversion
  HyperBasevector hbv;
  BinaryReader::readFile(fn, &hbv);
  std::cout << "Objetos en el hbv: " << hbv.EdgeObjectCount() << std::endl;

  vec<int> inv;
  inv.clear();
  hbv.Involution(inv);

  std::cout << "despues del inv 1" << std::endl;
  ReadPathVec pathsr;
  VecULongVec paths_inv;
  std::cout << "despues del inv 2" << std::endl;
//  LoadReadPathVec(pathsr, argv[4]);
  std::cout << "despues del inv 3" << std::endl;
//  paths_inv.clear();
  std::cout << "despues del inv 4" << std::endl;
//  invert(pathsr, paths_inv, hbv.EdgeObjectCount());

  std::cout << "despues del inv" << std::endl;
  auto edges = hbv.Edges();
  std::cout << "Size of the paths vector" << pathsr.size() <<" , inverse: " << paths_inv.size() << std::endl;

//  /* ----- TenXpather part ----- */
//  auto reads = dataMag.mag["TEX"]->rReads;
//  std::cout << "Reads already loaded..." << std::endl;
//  std::cout << reads.size() << " Reads in the vector" << std::endl;
//
//  // Create the paths and invert them
//  std::cout << "Starting tenxPather..." << std::endl;
//  TenXPather txp (reads, hbv, inv, 5, edges, pathsr, paths_inv);
//
//  std::cout<< Date() << " Map creation." << std::endl;
//  txp.createEmptyMap(&hbv);
//
//  std::cout<< Date() << " Map filling with reads..." << std::endl;
//  txp.reads2kmerTagMap();
//
//  std::cout<< Date() << " Map filling with reads done..." << std::endl;
//  txp.kmerTagDensity();
//
//  // Pathfinder
//  std::cout<< Date() << " Starting pathfinder..." << std::endl;
//  std::cout<< Date() << " done pathfinder..." << std::endl;
//  txp.solve_region_using_TenX(1000, true);
//  /* ----- TenXpather part ----- */


  auto reads = dataMag.mag["PB1"]->bases;

  PacbioPather pbp(reads, hbv, inv, 5, edges, pathsr, paths_inv);
  pbp.Hbv2Map(&hbv);
  pbp.mapReads();
  pbp.untangle_complex_in_out_choices(1000, true);

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