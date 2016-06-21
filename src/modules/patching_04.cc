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
#include "tclap/CmdLine.h"


int patcher(const String work_dir, const string prefix, uint NUM_THREADS, int MAX_MEM_GB){
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
  vec<String> subsam_names =  { "C" };
  vec<int64_t> subsam_starts = { 0 };
  vecbvec bases;
  bases.ReadAll( work_dir + "/frag_reads_orig.fastb" );

  VecPQVec quals;
  quals.ReadAll( work_dir + "/frag_reads_orig.qualp" );

  //XXX TODO: Load quals
  // Load K=200 HBV and reda paths
  BinaryReader::readFile( work_dir + "/" + prefix + ".pc.large.hbv", &hb);
  hb.Involution(inv2);
  //BinaryReader::readFile( work_dir + "/" + prefix + ".pc.large.inv", &inv2 );
  paths2.ReadAll( work_dir + "/" + prefix + ".pc.large.paths" );
  VecULongVec paths2_index;
  invert(paths2, paths2_index, hb.EdgeObjectCount());

  //int nedges = hb.EdgeObkectCount();

  // Assemble gaps, scope new stuff
  vecbvec new_stuff;  

 
  bool EXTEND=False;
  bool ANNOUNCE=False;
  bool KEEP_ALL_LOCAL=False;
  bool CONSERVATIVE_KEEP=False;
  bool INJECT=False;
  bool LOCAL_LAYOUT=False;
  const String DUMP_LOCAL="";
  int K2_FLOOR=0;
  int DUMP_LOCAL_LROOT=-1;
  int DUMP_LOCAL_RROOT=-1;
  bool CYCLIC_SAVE=True;
  int A2V=5;
  int GAP_CAP=-1;
  int MAX_PROX_LEFT=400;
  int MAX_PROX_RIGHT=400;
  int MAX_BPATHS=100000;

  AssembleGaps2( hb, inv2, paths2, paths2_index, bases, quals, work_dir, EXTEND, ANNOUNCE, KEEP_ALL_LOCAL, CONSERVATIVE_KEEP, INJECT, LOCAL_LAYOUT, DUMP_LOCAL, K2_FLOOR, DUMP_LOCAL_LROOT, DUMP_LOCAL_RROOT, new_stuff, CYCLIC_SAVE, A2V, GAP_CAP, MAX_PROX_LEFT, MAX_PROX_RIGHT, MAX_BPATHS );
  
  int MIN_GAIN=5;
  //const String TRACE_PATHS="{}";
  const vec<int> TRACE_PATHS;
  int EXT_MODE=1;
  
  AddNewStuff( new_stuff, hb, inv2, paths2, bases, quals, MIN_GAIN, TRACE_PATHS, work_dir, EXT_MODE );
  PartnersToEnds( hb, paths2, bases, quals );
  TestInvolution( hb, inv2 );

  // Write patched files to disk XXX TODO: corregir los paths para que no se superpongan entre ellos y con los anteriores
  //BinaryWriter::writeFile( work_dir + "/new_stuff", new_stuff);
  //vecbvec(hb.Edges().begin(), hb.Edges().end()).WriteAll( work_dir + "/" + prefix + ".patched.fastb");
  BinaryWriter::writeFile( work_dir + "/" + prefix + ".patched.hbv", hb );
  //BinaryWriter::writeFile( work_dir + "/" + prefix + ".patched.hbx", HyperBasevectorX(hb) );
  //BinaryWriter::writeFile( work_dir + "/" + prefix + ".patched.inv", inv2 );
  paths2.WriteAll( work_dir + "/" + prefix + ".patched.paths");
 
  /*VecULongVec invPaths;
  invert ( paths2, invPaths, hb.EdgeObjectCount() );
  invPaths.WriteAll( work_dir + "/" + prefix + ".patched.paths.inv" );
  
  Validate(hb, inv2, paths2); 

  hb.DumpFasta( work_dir + "/" + prefix + ".patched.fasta", False );*/
  return 0;
}

int main(int argc, const char* argv[]){
  std::string out_prefix;
  std::string out_dir;
  unsigned int threads;
  int max_mem;

  //========== Command Line Option Parsing ==========

  std::cout<<"Welcome to w2rap-contigger::04_patching"<<std::endl;
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



  patcher(out_dir, out_prefix, threads, max_mem );

  return 0;
}
