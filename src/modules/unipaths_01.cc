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


int create_unipaths(const String work_dir, const string prefix, const string READS, uint NUM_THREADS, int MAX_MEM_GB){
  /* Create unipaths from the reads 
   * XXX TODO: Document the function variables */

  // Set computational limits (XXX TODO: putin a separate source to import in different code)
  std::cout << "Set computational resources " << std::endl;
  PrintSysInfo();
  SetThreads(NUM_THREADS, False);
  int64_t max_bytes = int64_t(round(MAX_MEM_GB * 1024.0 * 1024.0 * 1024.0));
  SetMaxMemory(max_bytes);

  // Create directory tree (XXX TODO: Create a program that setups the dir tree)

  // Extraxt reads, convert reads to HyperBasevector and paths
  std::cout << "Object manager " << std::endl;
  vecbvec bases;
  ObjectManager<MasterVec<PQVec>> quals (work_dir + "/frag_reads_orig.qualp"); // XXX TODO: path to frag_reads_orig.qualp // VecPQvec is replaced for MasterVec<PQVec> using directive in feudal/PQvec.h file check
  std::cout << "Final object manager " << std::endl;
   
  // Variables to run ExtractReads function
  const string blank_string_arg = "";
  String blank_string_feudal = ""; // Feudal string
  String X = "all";
  String READS_F = READS;
  int ZZ = 0;
  std::cout << "Start to extract reads " << std::endl;
  vec<String> regions;
  vec<String> subsam_names;

  Samples( blank_string_feudal, "", X, blank_string_feudal, blank_string_feudal, ZZ, blank_string_arg, READS_F, subsam_names );
  vec<int64_t> subsam_starts( subsam_names.size( ), 0 );
  
  ExtractReads(blank_string_arg, blank_string_arg, READS, blank_string_feudal, -1, regions, work_dir, work_dir, False, False, False, subsam_names, subsam_starts, &bases, quals );
  std::cout << "Final to extract reads " << std::endl;

  BinaryWriter::writeFile( work_dir + "/subsam.starts", subsam_starts );
  BinaryWriter::writeFile( work_dir + "/subsam.names", subsam_names );

  // XXX TODO: Internal error check in original source check also here (?)
  // XXX TODO: add stats calculation and output like in original source if necesary
  
  return 0; 
}


int main(int argc, char *argv[]){
  // ./unipaths_01 <working_dir> '<read_files_R1.fastq,read_files_R2.fastq>' <prefix> <ncpus> <mem_GB>
  const String work_dir = argv[1];
  const string READS = argv[2];
  const string prefix = argv[3];
  int ncpus = atoi(argv[4]);
  int mem_mb = atoi(argv[5]);
  create_unipaths( work_dir, prefix, READS, ncpus, mem_mb );

  return 0;
}
