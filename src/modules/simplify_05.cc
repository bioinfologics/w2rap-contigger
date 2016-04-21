#include "FastIfstream.h"
#include "FetchReads.h"
#include "Intvector.h"
#include "MainTools.h"
#include "PairsManager.h"
#include "ParallelVecUtilities.h"
#include "feudal/ObjectManager.h"
#include "feudal/PQVec.h"
#include "lookup/LookAlign.h"
#include "paths/HyperBasevector.h"
#include "paths/RemodelGapTools.h"
#include "paths/long/BuildReadQGraph.h"
#include "paths/long/PlaceReads0.h"
#include "paths/long/SupportedHyperBasevector.h"
#include "paths/long/large/AssembleGaps.h"
#include "paths/long/large/Simplify.h"


int simplify(const String work_dir, const string prefix, uint NUM_THREADS, int MAX_MEM_GB){
  /*Do the patching steps 
   * XXX TODO: rename this part, patching it's not correct */
  
  // Setup system resources 
  std::cout << "Set computational resources " << std::endl;
  PrintSysInfo();
  SetThreads(NUM_THREADS, False);
  int64_t max_bytes = int64_t(round(MAX_MEM_GB * 1024.0 * 1024.0 * 1024.0));
  SetMaxMemory(max_bytes);

  // Proceed
  // Load necesary files
  vec<String> subsam_names;
  vec<int64_t> subsam_starts( subsam_names.size( ), 0 );
  BinaryReader::readFile(work_dir + "/subsam.starts", &subsam_starts);

  vecbvec bases;
  bases.ReadAll( work_dir + "/frag_reads_orig.fastb" );
  ObjectManager<MasterVec<PQVec>> quals ( work_dir + "/" + "frag_reads_orig.qualp");
  
  // Load K=200 HBV and reda paths
  HyperBasevector hb;
  vec<int> inv2;
  ReadPathVec paths2;
  BinaryReader::readFile( work_dir + "/" + prefix + ".patched.hbv", &hb);
  BinaryReader::readFile( work_dir + "/" + prefix + ".patched.inv", &inv2 );
  paths2.ReadAll( work_dir + "/" + prefix + ".patched.paths" );

 
  // variables to run Simplify 
  int MAX_SUPP_DEL=0;
  bool TAMP_EARLY_MIN=True;
  int MIN_RATIO2=8;
  int MAX_DEL2=200;
  bool PLACE_PARTNERS=False;
  bool ANALYZE_BRANCHES_VERBOSE2=False;
  const String TRACE_SEQ="";
  bool DEGLOOP=True;
  bool EXT_FINAL=True;
  int EXT_FINAL_MODE=1;
  bool PULL_APART_VERBOSE=False;
  //const String PULL_APART_TRACE="{}";
  const vec<int> PULL_APART_TRACE;
  int DEGLOOP_MODE=1;
  float DEGLOOP_MIN_DIST=2.5;
  bool IMPROVE_PATHS=True;
  bool IMPROVE_PATHS_LARGE=False;
  bool FINAL_TINY=True;
  bool UNWIND3=True;
  const String fin_dir = work_dir;

  Simplify( fin_dir, hb, inv2, paths2, bases, quals.load(), MAX_SUPP_DEL, TAMP_EARLY_MIN, MIN_RATIO2, MAX_DEL2, PLACE_PARTNERS, ANALYZE_BRANCHES_VERBOSE2, TRACE_SEQ, DEGLOOP, EXT_FINAL, EXT_FINAL_MODE, PULL_APART_VERBOSE, PULL_APART_TRACE, DEGLOOP_MODE, DEGLOOP_MIN_DIST, IMPROVE_PATHS, IMPROVE_PATHS_LARGE, FINAL_TINY, UNWIND3 );
  TestInvolution( hb, inv2 );
 
  // Write pre-scaffolded final assembly
  bases.destroy(); //XXXTODO: revisar si es realmente necesario
  BinaryWriter::writeFile( work_dir + "/" + prefix + ".fin.hbv", hb );
  BinaryWriter::writeFile( work_dir + "/" + prefix + ".fin.hbx", HyperBasevectorX(hb) );
   
  vecbasevector afin;
  for (int e=0; e < hb.EdgeObjectCount(); e++){ // XXX TODO: change this int for uint 32
    afin.push_back(hb.EdgeObject(e));
  }
  afin.WriteAll(work_dir + "/" + prefix + ".fin.fastb");
  
  // For now, fix paths and write the and their inverse
  for( int i = 0; i < (int) paths2.size(); i++){ //XXX TODO: change this int for uint 32
    Bool bad=False;
    for (int j=0; j < (int) paths2[i].size(); j++)
      if (paths2[i][j] < 0) bad = True;  
    if (bad) paths2[i].resize(0);
  }
  paths2.WriteAll( work_dir + "/" + prefix + ".fin.paths");
  VecULongVec invPaths;
  invert (paths2, invPaths, hb.EdgeObjectCount());
  invPaths.WriteAll ( work_dir + "/" + prefix + ".fin.paths.inv");
  BinaryWriter::writeFile( work_dir + "/" + prefix + ".fin.inv", inv2 );
  hb.DumpFasta( work_dir + "/" +prefix+ ".fin.sequence.fasta", False ); 
 
  // Find lines and write files.
  vec<vec<vec<vec<int>>>> lines;
  int MAX_CELL_PATHS=50;
  int MAX_DEPTH=10;
  FindLines( hb, inv2, lines, MAX_CELL_PATHS, MAX_DEPTH );
  BinaryWriter::writeFile( work_dir + "/" +prefix+ ".fin.lines", lines );
  
  // XXX TODO: Solve the {} thingy, check if has any influence in the new code to run that integrated
  {  
  vec<int> llens, npairs;
  GetLineLengths( hb, lines, llens );
  GetLineNpairs( hb, inv2, paths2, lines, npairs );
  BinaryWriter::writeFile( work_dir + "/" +prefix+ ".fin.lines.npairs", npairs );
  
  vec<vec<covcount>> covs;
  ComputeCoverage( hb, inv2, paths2, lines, subsam_starts, covs );
  BinaryWriter::writeFile( work_dir + "/" +prefix+ ".fin.covs", covs );
  WriteLineStats( work_dir, lines, llens, npairs, covs ); 
  
  // Report CN stats
  double cn_frac_good = CNIntegerFraction(hb, covs);
  std::cout << "CN fraction good = " << cn_frac_good << std::endl;
  PerfStatLogger::log("cn_frac_good",ToString(cn_frac_good,2), "fraction of edges with CN near integer" );
  }
  
  // TestLineSymmetry( lines, inv2 );
  // Compute fragment distribution.
  FragDist( hb, inv2, paths2, work_dir + "/" +prefix+ ".fin.frags.dist" );
  
  return 0;
}

int main(int argc, const char* argv[]){
  /*
  bool TEST=False;
  if (TEST){
    std::cout << "Arrancando" << std::endl;
    const String work_dir = "/tgac/workarea/Research-Groups/RG-Bernardo-Clavijo/wheat/discovar_complete_run/sources/A-team/A-discovar/src/modules/testcase_def";
    simplify(work_dir, "testrun", 1, 16);
  } else {
    const String work_dir = "/tgac/workarea/Research-Groups/RG-Bernardo-Clavijo/cadenza/discovar_splited_uv2k2";
    simplify(work_dir, "cadenza_2st_run", 64, 6500);
  }
  */

  // ./simplify /tgac/workarea/Research-Groups/RG-Bernardo-Clavijo/cadenza/discovar_splited prefix 64 1500
  const String work_dir = argv[1];
  const string prefix = argv[2];
  int ncpus = atoi(argv[3]);
  int mem_mb = atoi(argv[4]);
  simplify( work_dir, prefix, ncpus, mem_mb );
  
  return 0;
}
