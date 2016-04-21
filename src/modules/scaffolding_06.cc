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
#include "paths/long/large/FinalFiles.h"
#include "paths/long/large/MakeGaps.h"
#include "paths/long/large/Samples.h"


int scaffolding(const String work_dir, const string prefix, uint NUM_THREADS, int MAX_MEM_GB){
  /*Do the patching steps 
   * XXX TODO: rename this part, patching it's not correct */
  
  // Setup system resources 
  std::cout << "Set computational resources " << std::endl;
  PrintSysInfo();
  SetThreads(NUM_THREADS, False);
  int64_t max_bytes = int64_t(round(MAX_MEM_GB * 1024.0 * 1024.0 * 1024.0));
  SetMaxMemory(max_bytes);

  // Proceed
  HyperBasevector hb;
  vec<int> inv2;
  ReadPathVec paths2;

  // Load necesary files
  vecbvec bases;
  vec<String> subsam_names;
  vec<int64_t> subsam_starts( subsam_names.size( ), 0 );
  std::cout << "Loading fastb, subsam and qualp files" << std::endl;
  bases.ReadAll( work_dir + "/frag_reads_orig.fastb" );
  BinaryReader::readFile(work_dir + "/subsam.starts", &subsam_starts);
  ObjectManager<MasterVec<PQVec>> quals ( work_dir + "/" + "frag_reads_orig.qualp");

  //XXX TODO: Load quals
  // Load K=200 HBV and reda paths
  std::cout << "loading hbv, inv and paths files" << std::endl;
  BinaryReader::readFile( work_dir + "/" + prefix + ".fin.hbv", &hb);
  BinaryReader::readFile( work_dir + "/" + prefix + ".fin.inv", &inv2 );
  paths2.ReadAll( work_dir + "/" + prefix + ".fin.paths" );

  // Scaffold.
  //{    
  VecULongVec invPaths;
  invert( paths2, invPaths, hb.EdgeObjectCount( ) );
  int MIN_LINE=5000;
  int MIN_LINK_COUNT=3; //XXX TODO: this variable is the same as -w in soap??
  string FIN=prefix; //XXX TODO: i'm putting the prefix in the fin parameter for convenience, fix this
  bool SCAFFOLD_VERBOSE=False;
  bool GAP_CLEANUP=True;
  MakeGaps( hb, inv2, paths2, invPaths, MIN_LINE, MIN_LINK_COUNT, work_dir, FIN, SCAFFOLD_VERBOSE, GAP_CLEANUP );
  //}
  
  // Carry out final analyses and write final assembly files.
  String final_dir = work_dir;
  int MAX_CELL_PATHS=50;
  int MAX_DEPTH=10;
  bool ALIGN_TO_GENOME=True;
  String EVALUATE="";
  bool EVALUATE_VERBOSE=False;
  String X="all";
  std::map<String,GapToyResults> res;
  string SAMPLE="";
  String species;
  String SELECT_FRAC="";
  int READS_TO_USE=-1;
  string DATASET="1";
  String READS="";
    
  Samples( species, SAMPLE, X, EVALUATE, SELECT_FRAC, READS_TO_USE, DATASET, READS, subsam_names );
  int PAD=30000;
  bool all=False;
  vec<String> regions;
  vec<int> fosmids;
  String F;
  //DefineRegions( X, fosmids, regions, res, PAD, all, SAMPLE, EVALUATE, F );

  vecbasevector G;
  bool SAVE_FASTA=True;
  FinalFiles( hb, inv2, paths2, subsam_names, subsam_starts, work_dir, final_dir, MAX_CELL_PATHS, MAX_DEPTH, ALIGN_TO_GENOME, EVALUATE, EVALUATE_VERBOSE, X, res, SAMPLE, species, fosmids, G, SAVE_FASTA );
                                            
  // Done.
  std::cout << "peak mem usage = " << PeakMemUsageGBString( ) << ", ";
  std::cout << "final checksum = " << hb.CheckSum( ) << std::endl;
                                                    
  return 0;
}

int main(int argc, const char* argv[]){
  /*
  //std::cout << "Arrancando" << std::endl;
  //const String work_dir = "/tgac/workarea/Research-Groups/RG-Bernardo-Clavijo/wheat/discovar_complete_run/sources/A-team/A-discovar/src/modules/testcase_def";
  //scaffolding( work_dir, "testrun", 1, 16);
  bool TEST=False;
  if (TEST){
    std::cout << "Arrancando" << std::endl;
    const String work_dir = "/tgac/workarea/Research-Groups/RG-Bernardo-Clavijo/wheat/discovar_complete_run/sources/A-team/A-discovar/src/modules/testcase_def";
    scaffolding(work_dir, "testrun", 1, 16);
  } else {
    const String work_dir = "/tgac/workarea/Research-Groups/RG-Bernardo-Clavijo/cadenza/discovar_splited";
    scaffolding(work_dir, "cadenza_2st_run", 64, 6500);
  }
  */

  // ./scaffolding /tgac/workarea/Research-Groups/RG-Bernardo-Clavijo/cadenza/discovar_splited prefix 64 1500
  const String work_dir = argv[1];
  const string prefix = argv[2];
  int ncpus = atoi(argv[3]);
  int mem_mb = atoi(argv[4]);
  scaffolding( work_dir, prefix, ncpus, mem_mb );

  return 0;
}
