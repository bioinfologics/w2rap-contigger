//
// Created by Bernardo Clavijo (TGAC) on 21/06/2016.
//
#include <paths/long/large/Repath.h>
#include <paths/long/large/Clean200.h>
#include <paths/long/large/Simplify.h>
#include <paths/long/large/MakeGaps.h>
#include <paths/long/large/FinalFiles.h>
#include "FastIfstream.h"
#include "FetchReads.h"
#include "MainTools.h"
#include "PairsManager.h"
#include "ParallelVecUtilities.h"
#include "feudal/PQVec.h"
#include "lookup/LookAlign.h"
#include "paths/HyperBasevector.h"
#include "paths/RemodelGapTools.h"
#include "paths/long/BuildReadQGraph.h"
#include "paths/long/PlaceReads0.h"
#include "paths/long/SupportedHyperBasevector.h"
#include "paths/long/large/AssembleGaps.h"
#include "paths/long/large/ExtractReads.h"
#include "tclap/CmdLine.h"
#include <sys/types.h>
#include <sys/stat.h>

int main(const int argc, const char * argv[]) {

    std::string out_prefix;
    std::string read_files;
    std::string out_dir;
    unsigned int threads;
    int max_mem;
    unsigned int small_K, large_K, min_size;
    std::vector<unsigned int> allowed_k = {60, 64, 72, 80, 84, 88, 96, 100, 108, 116, 128, 136, 144, 152, 160, 168, 172,
                                           180, 188, 192, 196, 200, 208, 216, 224, 232, 240, 260, 280, 300, 320, 368,
                                           400, 440, 460, 500, 544, 640};
    bool extend_paths,dump_all;

    //========== Command Line Option Parsing ==========

    std::cout << "Welcome to w2rap-contigger" << std::endl;
    try {
        TCLAP::CmdLine cmd("", ' ', "0.1");
        TCLAP::ValueArg<unsigned int> threadsArg("t", "threads",
             "Number of threads on parallel sections (default: 4)", false, 4, "int", cmd);
        TCLAP::ValueArg<unsigned int> max_memArg("m", "max_mem",
             "Maximum memory in GB (soft limit, impacts performance, default 10000)", false, 10000, "int", cmd);

        TCLAP::ValueArg<std::string> read_filesArg("r", "read_files",
             "Input sequences (reads) files ", true, "", "file1.fastq,file2.fastq", cmd);

        TCLAP::ValueArg<std::string> out_dirArg("o", "out_dir", "Output dir path", true, "", "string", cmd);
        TCLAP::ValueArg<std::string> out_prefixArg("p", "prefix",
             "Prefix for the output files", true, "", "string", cmd);

        TCLAP::ValuesConstraint<unsigned int> largeKconst(allowed_k);
        TCLAP::ValueArg<unsigned int> large_KArg("K", "large_k",
             "Large k (default: 200)", false, 200, &largeKconst, cmd);
        TCLAP::ValueArg<unsigned int> small_KArg("k", "small_k",
                                                 "Small k (default: 60)", false, 60, &largeKconst, cmd);
        TCLAP::ValueArg<unsigned int> minSizeArg("s", "min_size",
             "Min size of disconnected elements on large_k graph (in kmers, default: 0=no min)", false, 0, "int", cmd);
        TCLAP::ValueArg<bool>         pathExtensionArg        ("","extend_paths",
             "Enable extend paths on repath (experimental)", false,false,"bool",cmd);
        TCLAP::ValueArg<bool>         dumpAllArg        ("","dump_all",
                                                               "Dump all intermediate files", false,false,"bool",cmd);
        cmd.parse(argc, argv);

        // Get the value parsed by each arg.
        out_dir = out_dirArg.getValue();
        out_prefix = out_prefixArg.getValue();
        read_files = read_filesArg.getValue();
        threads = threadsArg.getValue();
        max_mem = max_memArg.getValue();
        large_K = large_KArg.getValue();
        small_K = small_KArg.getValue();
        min_size = minSizeArg.getValue();
        extend_paths=pathExtensionArg.getValue();
        dump_all=dumpAllArg.getValue();

    } catch (TCLAP::ArgException &e)  // catch any exceptions
    {
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
        return 1;
    }

    //Check directory exists:


    struct stat info;

    if (stat(out_dir.c_str(), &info) != 0 || !(info.st_mode & S_IFDIR)) {
        std::cout << "Output directory doesn't exist, or is not a directory: " << out_dir << std::endl;
        return 1;
    }

    //========== Main Program Begins ======

    //== Set computational resources ===
    SetThreads(threads, False);
    SetMaxMemory(int64_t(round(max_mem * 1024.0 * 1024.0 * 1024.0)));


    //== Load reads (and saves in binary format) ======
    vecbvec bases;
    VecPQVec quals;

    vec<String> subsam_names = {"C"};
    vec<int64_t> subsam_starts = {0};

    //create_unipaths( out_dir, out_prefix, read_files, threads, max_mem );
    std::cout << "--== Reading input files ==--" << std::endl;
    ExtractReads(read_files, out_dir, subsam_names, subsam_starts, &bases, &quals);
    std::cout << "Reading input files DONE!" << std::endl<< std::endl<< std::endl;
    //TODO: add an option to dump the reads
    if (dump_all) {
        std::cout << "Dumping reads in fastb/qualp format..." << std::endl;
        bases.WriteAll(out_dir + "/frag_reads_orig.fastb");
        quals.WriteAll(out_dir + "/frag_reads_orig.qualp");
        std::cout << "   DONE!" << std::endl;
    }


    //== Read QGraph, and repath (k=60, k=200 (and saves in binary format) ======


    bool FILL_JOIN = False;
    bool SHORT_KMER_READ_PATHER = False;
    bool RQGRAPHER_VERBOSE = False;
    vec<int> inv;
    HyperBasevector hbvr;
    ReadPathVec pathsr;
    {//This scope-trick to invalidate old data is dirty
        HyperBasevector hbv;
        ReadPathVec paths;
        std::cout << "--== Building first graph ==--" << std::endl;
        buildReadQGraph(bases, quals, FILL_JOIN, FILL_JOIN, 7, 3, .75, 0, "", True, SHORT_KMER_READ_PATHER, &hbv, &paths,
                        small_K);
        FixPaths(hbv, paths); //TODO: is this even needed?
        std::cout << "Building first graph DONE!" << std::endl<< std::endl<< std::endl;

        std::cout << "--== Repathing to second graph ==--" << std::endl;
        vecbvec edges(hbv.Edges().begin(), hbv.Edges().end());
        hbv.Involution(inv);

        const string run_head = out_dir + "/" + out_prefix;

        pathsr.resize(paths.size());
        RepathInMemory(hbv, edges, inv, paths, hbv.K(), large_K, hbvr, pathsr, True, True, extend_paths);
        std::cout << "Repathing to second graph DONE!" << std::endl<< std::endl<< std::endl;
    }

    //== Clean ======
    std::cout << "--== Cleaning graph ==--" << std::endl;
    inv.clear();
    hbvr.Involution(inv);
    int CLEAN_200_VERBOSITY=0;
    int CLEAN_200V=3;
    Clean200x( hbvr, inv, pathsr, bases, quals, CLEAN_200_VERBOSITY, CLEAN_200V, min_size );
    if (dump_all) {
        BinaryWriter::writeFile(out_dir + "/" + out_prefix + ".pc.large.hbv", hbvr);
        pathsr.WriteAll(out_dir + "/" + out_prefix + ".pc.large.paths");
    }

    std::cout << "Cleaning graph DONE!" << std::endl<< std::endl<< std::endl;
    //== Patching ======

    std::cout << "--== Assembling gaps ==--" << std::endl;
    VecULongVec paths_inv;
    invert(pathsr, paths_inv, hbvr.EdgeObjectCount());

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

    AssembleGaps2( hbvr, inv, pathsr, paths_inv, bases, quals, out_dir, EXTEND, ANNOUNCE, KEEP_ALL_LOCAL, CONSERVATIVE_KEEP, INJECT, LOCAL_LAYOUT, DUMP_LOCAL, K2_FLOOR, DUMP_LOCAL_LROOT, DUMP_LOCAL_RROOT, new_stuff, CYCLIC_SAVE, A2V, GAP_CAP, MAX_PROX_LEFT, MAX_PROX_RIGHT, MAX_BPATHS );

    int MIN_GAIN=5;
    //const String TRACE_PATHS="{}";
    const vec<int> TRACE_PATHS;
    int EXT_MODE=1;

    AddNewStuff( new_stuff, hbvr, inv, pathsr, bases, quals, MIN_GAIN, TRACE_PATHS, out_dir, EXT_MODE );
    PartnersToEnds( hbvr, pathsr, bases, quals );
    std::cout << "Assembling gaps DONE!" << std::endl<< std::endl<< std::endl;

    std::cout << "--== Simplifying ==--" << std::endl;
    //==Simplify
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

    Simplify( out_dir, hbvr, inv, pathsr, bases, quals, MAX_SUPP_DEL, TAMP_EARLY_MIN, MIN_RATIO2, MAX_DEL2, PLACE_PARTNERS, ANALYZE_BRANCHES_VERBOSE2, TRACE_SEQ, DEGLOOP, EXT_FINAL, EXT_FINAL_MODE, PULL_APART_VERBOSE, PULL_APART_TRACE, DEGLOOP_MODE, DEGLOOP_MIN_DIST, IMPROVE_PATHS, IMPROVE_PATHS_LARGE, FINAL_TINY, UNWIND3 );

    // For now, fix paths and write the and their inverse
    for( int i = 0; i < (int) pathsr.size(); i++){ //XXX TODO: change this int for uint 32
        Bool bad=False;
        for (int j=0; j < (int) pathsr[i].size(); j++)
            if (pathsr[i][j] < 0) bad = True;
        if (bad) pathsr[i].resize(0);
    }
    // TODO: this is "bj making sure the inversion still works", but shouldn't be required
    paths_inv.clear();
    invert(pathsr,paths_inv,hbvr.EdgeObjectCount());
    std::cout << "Simplifying DONE!" << std::endl<< std::endl<< std::endl;

    std::cout << "--== Finding lines ==--" << std::endl;
    // Find lines and write files.
    vec<vec<vec<vec<int>>>> lines;
    int MAX_CELL_PATHS=50;
    int MAX_DEPTH=10;
    FindLines( hbvr, inv, lines, MAX_CELL_PATHS, MAX_DEPTH );
    BinaryWriter::writeFile( out_dir + "/" + out_prefix + ".fin.lines", lines );

    // XXX TODO: Solve the {} thingy, check if has any influence in the new code to run that integrated
    {
        vec<int> llens, npairs;
        GetLineLengths( hbvr, lines, llens );
        GetLineNpairs( hbvr, inv, pathsr, lines, npairs );
        BinaryWriter::writeFile( out_dir + "/" + out_prefix+ ".fin.lines.npairs", npairs );

        vec<vec<covcount>> covs;
        ComputeCoverage( hbvr, inv, pathsr, lines, subsam_starts, covs );
        //BinaryWriter::writeFile( work_dir + "/" +prefix+ ".fin.covs", covs );
        //WriteLineStats( work_dir, lines, llens, npairs, covs );

        // Report CN stats
        double cn_frac_good = CNIntegerFraction(hbvr, covs);
        std::cout << "CN fraction good = " << cn_frac_good << std::endl;
        PerfStatLogger::log("cn_frac_good",ToString(cn_frac_good,2), "fraction of edges with CN near integer" );
    }

    // TestLineSymmetry( lines, inv2 );
    // Compute fragment distribution.
    FragDist( hbvr, inv, pathsr, out_dir + "/" +out_prefix+ ".fin.frags.dist" );

    std::cout << "Finding lines DONE!" << std::endl<< std::endl<< std::endl;


    //== Scaffolding
    std::cout << "--== Scaffolding ==--" << std::endl;
    int MIN_LINE=5000;
    int MIN_LINK_COUNT=3; //XXX TODO: this variable is the same as -w in soap??

    bool SCAFFOLD_VERBOSE=False;
    bool GAP_CLEANUP=True;
    MakeGaps( hbvr, inv, pathsr, paths_inv, MIN_LINE, MIN_LINK_COUNT, out_dir, out_prefix, SCAFFOLD_VERBOSE, GAP_CLEANUP );

    std::cout << "--== Scaffolding DONE!" << std::endl<< std::endl<< std::endl;
    // Carry out final analyses and write final assembly files.

    vecbasevector G;
    FinalFiles( hbvr, inv, pathsr, subsam_names, subsam_starts, out_dir, MAX_CELL_PATHS, MAX_DEPTH, G);


    return 0;
}

