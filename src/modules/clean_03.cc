#include <omp.h>
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
#include "paths/long/large/Clean200.h"
#include "tclap/CmdLine.h"


int clean_graph(const String work_dir, const string file_prefix, uint NUM_THREADS, int MAX_MEM_GB) {

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

    // Clean (original comment) XXX TODO: porque redeclaration (??) revisar esto!!
    HyperBasevector hb;
    vec<int> inv;
    BinaryReader::readFile( work_dir + "/"+ file_prefix +".large.hbv", &hb );
    BinaryReader::readFile( work_dir + "/"+ file_prefix +".large.inv", &inv );
    ReadPathVec paths( work_dir + "/"+ file_prefix +".large.paths" );

    // select the cleaning version (default value in discovar is 2) we always use the Clean200 option !!
    std::cout << "Cleaning!! " << std::endl;
    int CLEAN_200_VERBOSITY=0;
    int CLEAN_200V=2;
    Bool REMOVE_TINY=False;
    if ( CLEAN_200V <= 2 ) {
        Clean200( hb, inv, paths, bases, quals.load(), CLEAN_200_VERBOSITY, CLEAN_200V, REMOVE_TINY );
    } else {
        Clean200x( hb, inv, paths, bases, quals.load(), CLEAN_200_VERBOSITY, CLEAN_200V, REMOVE_TINY );
    }
    std::cout << "Done cleaning!! " << std::endl;

    std::cout << "Writing output files " << std::endl;
    // Write files.
    BinaryWriter::writeFile( work_dir + "/" + file_prefix + ".pc.large.hbv", hb );
    BinaryWriter::writeFile( work_dir + "/" + file_prefix + ".pc.large.hbx", HyperBasevectorX(hb) );
    BinaryWriter::writeFile( work_dir + "/" + file_prefix + ".pc.large.inv", inv );
    paths.WriteAll( work_dir + "/" + file_prefix + ".pc.large.paths" );

    vecbasevector edges( hb.EdgeObjectCount( ) );
    for ( int e = 0; e < hb.EdgeObjectCount( ); e++ ) {
        edges[e] = hb.EdgeObject(e);
    }
    edges.WriteAll( work_dir + "/" + file_prefix + ".pc.large.fastb" );
    std::cout << "Done writing output files " << std::endl;

    if ( true ) {  // scope paths_index (original comment ???)
        VecULongVec paths_index;
        invert( paths, paths_index, hb.EdgeObjectCount( ) );
        paths_index.WriteAll( work_dir + "/" + file_prefix + ".pc.large.paths.inv" );
    }

    // Unload quals
    quals.unload();

    hb.DumpFasta( work_dir + "/" + file_prefix + ".pc.fasta", False );
    return 0;
}

int main(int argc, const char* argv[]) {
    std::string out_prefix;
    std::string out_dir;
    unsigned int threads;
    int max_mem;

    //========== Command Line Option Parsing ==========

    std::cout<<"Welcome to w2rap-contigger::03_clean"<<std::endl;
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

    } catch (TCLAP::ArgException &e) { // catch any exceptions
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
        return 1;
    }



    clean_graph(out_dir, out_prefix, threads, max_mem );

    return 0;
}
