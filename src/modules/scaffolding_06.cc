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
#include "tclap/CmdLine.h"


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

  // Load necessary files
  vec<String> subsam_names;
  vec<int64_t> subsam_starts( subsam_names.size( ), 0 );
  std::cout << "Loading subsam file" << std::endl;
  BinaryReader::readFile(work_dir + "/subsam.starts", &subsam_starts);

  //XXX TODO: Load quals
  // Load K=200 HBV and reda paths
  std::cout << "loading hbv, inv and paths files" << std::endl;
  BinaryReader::readFile( work_dir + "/" + prefix + ".fin.hbv", &hb);
  BinaryReader::readFile( work_dir + "/" + prefix + ".fin.inv", &inv2 );
  paths2.ReadAll( work_dir + "/" + prefix + ".fin.paths" );

  // Scaffold.

  VecULongVec invPaths;
  invert( paths2, invPaths, hb.EdgeObjectCount( ) );
  int MIN_LINE=5000;
  int MIN_LINK_COUNT=3; //XXX TODO: this variable is the same as -w in soap??
  string FIN=prefix; //XXX TODO: i'm putting the prefix in the fin parameter for convenience, fix this
  bool SCAFFOLD_VERBOSE=False;
  bool GAP_CLEANUP=True;
  MakeGaps( hb, inv2, paths2, invPaths, MIN_LINE, MIN_LINK_COUNT, work_dir, FIN, SCAFFOLD_VERBOSE, GAP_CLEANUP );

  
  // Carry out final analyses and write final assembly files.
  int MAX_CELL_PATHS=50;
  int MAX_DEPTH=10;


  vecbasevector G;
  FinalFiles( hb, inv2, paths2, subsam_names, subsam_starts, work_dir, MAX_CELL_PATHS, MAX_DEPTH, G);
                                            
  // Done.
#ifdef __linux
  std::cout << "peak mem usage = " << PeakMemUsageGBString( ) << ", ";
#endif
  std::cout << "final checksum = " << hb.CheckSum( ) << std::endl;
                                                    
  return 0;
}

int main(int argc, const char* argv[]){
  std::string out_prefix;
  std::string out_dir;
  unsigned int threads;
  int max_mem;

  //========== Command Line Option Parsing ==========

  std::cout<<"Welcome to w2rap-contigger::06_scaffolding"<<std::endl;
  try {
    TCLAP::CmdLine cmd("", ' ', "0.1 alpha");

    TCLAP::ValueArg<std::string> out_dirArg     ("o","out_dir",     "Output dir path",           true,"","string",cmd);
    TCLAP::ValueArg<std::string> out_prefixArg     ("p","prefix",     "Prefix for the output files",           true,"","string",cmd);
    TCLAP::ValueArg<unsigned int>         threadsArg        ("t","threads",        "Number of threads on parallel sections (default: 4)", false,4,"int",cmd);
    TCLAP::ValueArg<unsigned int>         max_memArg       ("m","max_mem",       "Maximum memory in GB (soft limit, impacts performance, default 10000)", false,10000,"int",cmd);

    cmd.parse( argc, argv );

    // Get the value parsed by each arg.
    out_dir=out_dirArg.getValue();
    out_prefix=out_prefixArg.getValue();
    threads=threadsArg.getValue();
    max_mem=max_memArg.getValue();

  } catch (TCLAP::ArgException &e)  // catch any exceptions
  { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; return 1;}



  scaffolding (out_dir, out_prefix, threads, max_mem );

  return 0;
}