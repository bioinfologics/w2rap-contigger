//
// Created by Gonzalo Garcia (TGAC) on 28/10/2016.
//

#include <iostream>
#include <kmers/kmatch/KMatch.h>
#include "paths/long/large/ExtractReads.h"

#include "pacbio/LongRead_pather.h"

#include <paths/long/large/Simplify.h>

#include "TenX/TenX_pather.h"
#include "TenX/PathFinder_tx.h"


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

//  std::ofstream inversions_file("/Users/ggarcia/Documents/notebooks/edgeinvolutions.txt", std::ios::out);
//  int cont = 0;
//  std::cout << "Edges ad involutions" << std::endl;
//  for (auto a: inv){
////    std::cout << cont << ":" << a << std::endl;
//    inversions_file << cont << ":" << a << std::endl;
//    cont++;
//  }
//  inversions_file.close();
//  std::cout << "Edges ad involutions" << std::endl;

  ReadPathVec pathsr;
  VecULongVec paths_inv;
  invert(pathsr, paths_inv, hbv.EdgeObjectCount());
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
//  txp.solve_region_using_TenX(5000, true);

//  /* ----- TenXpather part ----- */

  auto reads = dataMag.mag["PB1"]->bases;

  LongReadPather pbp(reads, hbv, inv, 5, edges, pathsr, paths_inv);
  std::cout << Date() << ": Object created" << std::endl;

  std::cout << Date() << ": Creating hbv map" << std::endl;
  pbp.Hbv2Map(hbv);
  std::cout << Date() << ": Hbv2Map ready ->("<< pbp.edgeMap.size() << " keys)" << std::endl;

  std::cout << Date() << ": Mapping reads" << std::endl;
  pbp.mapReads();
  std::cout << Date() << ": Reads mapped." << std::endl;


  pbp.solve_using_long_read(1000, true);

  std::cout<<"Removing Unneeded Vertices & Cleanup"<<std::endl;
  RemoveUnneededVertices2(hbv,inv,pathsr);
  Cleanup(hbv, inv, pathsr );

  std::string out_dir = "/Users/ggarcia/Documents/ecoli_test_dataset/test";
  auto bases = dataMag.mag["PE1"]->bases;
  auto quals = dataMag.mag["PE1"]->quals;

  //==Simplify
  int MAX_SUPP_DEL = 0;
  bool TAMP_EARLY_MIN = True;
  int MIN_RATIO2 = 8;
  int MAX_DEL2 = 200;
  bool ANALYZE_BRANCHES_VERBOSE2 = False;
  const String TRACE_SEQ = "";
  bool DEGLOOP = True;
  bool EXT_FINAL = True;
  int EXT_FINAL_MODE = 1;
  bool PULL_APART_VERBOSE = False;
  //const String PULL_APART_TRACE="{}";
  const vec<int> PULL_APART_TRACE;
  int DEGLOOP_MODE = 1;
  float DEGLOOP_MIN_DIST = 2.5;
  bool IMPROVE_PATHS = True;
  bool IMPROVE_PATHS_LARGE = False;
  bool FINAL_TINY = True;
  bool UNWIND3 = True;

  bool run_pathfinder = true;
  bool dump_pf = true;

  Simplify(out_dir, hbv, inv, pathsr, bases, quals, MAX_SUPP_DEL, TAMP_EARLY_MIN, MIN_RATIO2, MAX_DEL2,
           ANALYZE_BRANCHES_VERBOSE2, TRACE_SEQ, DEGLOOP, EXT_FINAL, EXT_FINAL_MODE,
           PULL_APART_VERBOSE, PULL_APART_TRACE, DEGLOOP_MODE, DEGLOOP_MIN_DIST, IMPROVE_PATHS,
           IMPROVE_PATHS_LARGE, FINAL_TINY, UNWIND3, run_pathfinder, dump_pf);



//  pathsr = pbp.mapReads();
//  VecULongVec invPaths;
//
//  invert(pathsr, invPaths, hbv.EdgeObjectCount());

//   pathfinders
//  PathFinder_tx(hbv, inv, pathsr, invPaths).unroll_loops(800);
//  PathFinder_tx(hbv,inv,pathsr,invPaths).untangle_pins();
//  PathFinder_tx(hbv,inv,pathsr,invPaths).untangle_complex_in_out_choices(700, true);

//  RemoveUnneededVertices2(hbv, inv, pathsr);
//  Cleanup(hbv, inv, pathsr);

//  BinaryWriter::writeFile("/Users/ggarcia/Documents/arabidopsis_testrun/pf_after_loops.hbv", hbv);
//  WriteReadPathVec(pathsr, "/Users/ggarcia/Documents/arabidopsis_testrun/pf_after_loops.paths");
}