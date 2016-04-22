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
#include "paths/long/large/Repath.h"
#include "tclap/CmdLine.h"



int qgraph_builder(const String work_dir, const string file_prefix, uint small_k, uint large_k, uint NUM_THREADS, int MAX_MEM_GB){
  
    // ********************** Set sys resources ******************
    // Set computational limits (XXX TODO: putin a separate source to import in different code)
    std::cout << "Set computational resources " << std::endl;
    PrintSysInfo();
    SetThreads(NUM_THREADS, False);
    int64_t max_bytes = int64_t(round(MAX_MEM_GB * 1024.0 * 1024.0 * 1024.0));
    SetMaxMemory(max_bytes);

    // ********************** Load files *************************
    std::cout << "Loading files " << std::endl;
    
    ObjectManager<MasterVec<PQVec>> quals ( work_dir + "/" + "frag_reads_orig.qualp");
    vec<String> subsam_names;
    vec<int64_t> subsam_starts( subsam_names.size( ), 0 );
    
    BinaryReader::readFile (work_dir + "/subsam.starts", &subsam_starts);
    BinaryReader::readFile (work_dir + "/subsam.names", &subsam_names);
    vecbvec bases;
    bases.ReadAll( work_dir + "/frag_reads_orig.fastb" );
    quals.load();

    // variables to run buildReadQGraph XXX TODO: Document one by one
    std::cout << "Starting to build Reads qgraph" << std::endl;
    HyperBasevector hbv;
    ReadPathVec paths;
    bool FILL_JOIN=False;
    bool SHORT_KMER_READ_PATHER=False;
    bool RQGRAPHER_VERBOSE=False;
    
    buildReadQGraph(bases, quals, FILL_JOIN, FILL_JOIN, 7, 3, .75, 0, "", True, SHORT_KMER_READ_PATHER, &hbv, &paths, RQGRAPHER_VERBOSE);
    std::cout << "Finish building Reads qgraph" << std::endl;
    
    vecbvec edges(hbv.Edges().begin(), hbv.Edges().end());
    vec<int> inv;
    hbv.Involution(inv);
    
    FixPaths( hbv, paths );
    
    // Save graph kmer 60
    BinaryWriter::writeFile( work_dir +"/"+ file_prefix +".small_K.hbv", hbv );
    BinaryWriter::writeFile( work_dir +"/"+ file_prefix +".small_K.hbx", HyperBasevectorX(hbv) );
    edges.WriteAll( work_dir +"/"+ file_prefix +".small_K.fastb" );
    BinaryWriter::writeFile( work_dir +"/"+ file_prefix +".small_K.inv", inv );
    
    std::cout << "Done loading " << std::endl;
    
    int64_t checksum_60 = hbv.CheckSum( );
    PRINT(checksum_60);
    
    // Variables to run Repath XXX TODO: Document one by one
    const string run_head = work_dir + "/" + file_prefix;
    
    quals.unload();
    Repath( hbv, edges, inv, paths, hbv.K(), 200, run_head+".200", True, True, False );
    
    hbv.DumpFasta( run_head + ".after_repath.fasta", False );
    return 0;
}

int main(int argc, const char* argv[]){
    String out_prefix;
    String out_dir;
    unsigned int small_K,large_K;
    unsigned int threads;
    int max_mem;

    //========== Command Line Option Parsing ==========

    std::cout<<"Welcome to w2rap-contigger::02_qgraph"<<std::endl;
    try {
        TCLAP::CmdLine cmd("", ' ', "0.1 alpha");

        TCLAP::ValueArg<std::string> out_dirArg     ("o","out_dir",     "Output dir path",           true,"","string",cmd);
        TCLAP::ValueArg<std::string> out_prefixArg     ("p","prefix",     "Prefix for the output files",           true,"","string",cmd);
        TCLAP::ValueArg<unsigned int>         small_KArg        ("k","small_k",        "Small k (default: 60)", false,60,"int",cmd);
        TCLAP::ValueArg<unsigned int>         large_KArg        ("K","large_k",        "Large k (default: 200)", false,200,"int",cmd);
        TCLAP::ValueArg<unsigned int>         threadsArg        ("t","threads",        "Number of threads on parallel sections (default: 4)", false,4,"int",cmd);
        TCLAP::ValueArg<unsigned int>         max_memArg       ("m","max_mem",       "Maximum memory in GB (soft limit, impacts performance, default 10000)", false,10000,"int",cmd);

        cmd.parse( argc, argv );

        // Get the value parsed by each arg.
        out_dir=out_dirArg.getValue();
        out_prefix=out_prefixArg.getValue();
        small_K=small_KArg.getValue();
        large_K=large_KArg.getValue();
        threads=threadsArg.getValue();
        max_mem=max_memArg.getValue();

    } catch (TCLAP::ArgException &e)  // catch any exceptions
    { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; return 1;}



    qgraph_builder(out_dir, out_prefix, small_K, large_K, threads, max_mem );

    return 0;
}
