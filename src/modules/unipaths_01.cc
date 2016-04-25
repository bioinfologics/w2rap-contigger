#include "FastIfstream.h"
#include "FetchReads.h"
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
#include "paths/long/large/ExtractReads.h"
#include "paths/long/large/Samples.h"
#include "tclap/CmdLine.h"

int create_unipaths(const string work_dir, const string prefix, const string READS, uint NUM_THREADS, int MAX_MEM_GB){
  /* Create unipaths from the reads 
   * XXX TODO: Document the function variables */

  // Set computational limits (XXX TODO: putin a separate source to import in different code)
  //std::cout << "Set computational resources " << std::endl;
  PrintSysInfo();
  SetThreads(NUM_THREADS, False);
  int64_t max_bytes = int64_t(round(MAX_MEM_GB * 1024.0 * 1024.0 * 1024.0));
  SetMaxMemory(max_bytes);

  // Create directory tree (XXX TODO: Create a program that setups the dir tree)

  // Extraxt reads, convert reads to HyperBasevector and paths
  //std::cout << "Object manager " << std::endl;
  vecbvec bases;
  ObjectManager<MasterVec<PQVec>> quals (work_dir + "/frag_reads_orig.qualp"); // XXX TODO: path to frag_reads_orig.qualp // VecPQvec is replaced for MasterVec<PQVec> using directive in feudal/PQvec.h file check
  //std::cout << "Final object manager " << std::endl;
   
  // Variables to run ExtractReads function
  const string blank_string_arg = "";
  String blank_string_feudal = ""; // Feudal string
  String X = "all";
  String READS_F = READS;
  int ZZ = 0;
  //std::cout << "Start to extract reads " << std::endl;
  vec<String> regions;
  vec<String> subsam_names;

  Samples( blank_string_feudal, "", X, blank_string_feudal, blank_string_feudal, ZZ, blank_string_arg, READS_F, subsam_names );
  vec<int64_t> subsam_starts( subsam_names.size( ), 0 );
  
  ExtractReads(blank_string_arg, blank_string_arg, READS, blank_string_feudal, -1, regions, work_dir, work_dir, False, False, False, subsam_names, subsam_starts, &bases, quals );
  //std::cout << "Final to extract reads " << std::endl;

  BinaryWriter::writeFile( work_dir + "/subsam.starts", subsam_starts );
  BinaryWriter::writeFile( work_dir + "/subsam.names", subsam_names );

  // XXX TODO: Internal error check in original source check also here (?)
  // XXX TODO: add stats calculation and output like in original source if necesary
  
  return 0; 
}


int main(const int argc, const char * argv[]){

  std::string out_prefix;
  std::string read_files;
  std::string out_dir;
  unsigned int threads;
  int max_mem;

  //========== Command Line Option Parsing ==========

  std::cout<<"Welcome to w2rap-contigger::01_unipaths"<<std::endl;
  try {
    TCLAP::CmdLine cmd("", ' ', "0.1 alpha");

    TCLAP::ValueArg<std::string> out_dirArg     ("o","out_dir",     "Output dir path",           true,"","string",cmd);
    TCLAP::ValueArg<std::string> out_prefixArg     ("p","prefix",     "Prefix for the output files",           true,"","string",cmd);

    TCLAP::ValueArg<std::string> read_filesArg  ("r","read_files",    "Input sequences (reads) files ", true,"","file1.fastq,file2.fastq",cmd);
    TCLAP::ValueArg<unsigned int>         threadsArg        ("t","threads",        "Number of threads on parallel sections (default: 4)", false,4,"int",cmd);
    TCLAP::ValueArg<unsigned int>         max_memArg       ("m","max_mem",       "Maximum memory in GB (soft limit, impacts performance, default 10000)", false,10000,"int",cmd);

    cmd.parse( argc, argv );

    // Get the value parsed by each arg.
    out_dir=out_dirArg.getValue();
    out_prefix=out_prefixArg.getValue();
    read_files=read_filesArg.getValue();
    threads=threadsArg.getValue();
    max_mem=max_memArg.getValue();

  } catch (TCLAP::ArgException &e)  // catch any exceptions
  { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; return 1;}

  //Check directory exists:
#include <sys/types.h>
#include <sys/stat.h>

  struct stat info;

  if( stat( out_dir.c_str(), &info ) != 0 || !( info.st_mode & S_IFDIR )) {
    std::cout<<"Output directory doesn't exist, or is not a directory: "<<out_dir<<std::endl;
    return 1;
  }

  create_unipaths( out_dir, out_prefix, read_files, threads, max_mem );

  return 0;
}
