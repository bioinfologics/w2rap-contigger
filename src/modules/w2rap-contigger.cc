//
// Created by Bernardo Clavijo (TGAC) on 21/06/2016.
//
#include <paths/long/large/Repath.h>
#include <paths/long/large/Clean200.h>
#include <paths/long/large/Simplify.h>
#include <paths/long/large/MakeGaps.h>
#include <paths/long/large/FinalFiles.h>
#include "FastIfstream.h"
#include "CoreTools.h"
#include "system/ParsedArgs.h"
#include "system/SysConf.h"
#include "PairsManager.h"
#include "ParallelVecUtilities.h"
#include "feudal/PQVec.h"
#include "paths/HyperBasevector.h"
#include "paths/RemodelGapTools.h"
#include "paths/long/BuildReadQGraph.h"
#include "paths/long/SupportedHyperBasevector.h"
#include "paths/long/large/AssembleGaps.h"
#include "paths/long/large/ExtractReads.h"
#include "deps/cxxopts/cxxopts.hpp"
#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>
#include <sys/time.h>
#include <paths/PathFinder.h>
#include <paths/long/large/ImprovePath.h>
#include <paths/long/large/GraphImprover.h>
#include <paths/long/large/ConsensusChecker.h>
#include <kmers/BigKMer.h>
#include "GFADump.h"
#include "util/OutputLog.h"
#include "SMR.h"
#include "KMerFreqFactory.h"
#include "kmers/KMer.h"
#include <omp.h>

 // Dummy defines, this should come from CMakeLists.txt
#ifndef GIT_COMMIT_HASH
#define GIT_COMMIT_HASH "dummy"
#endif
#ifndef GIT_BRANCH
#define GIT_BRANCH "dummy"
#endif

//TODO: stupid globals!

int MAX_CELL_PATHS = 50;
int MAX_DEPTH = 10;
const uint64_t GB = 1024*1024*1024;

void step_7DV(HyperBasevector &hbv,
            vec<int> &hbvinv,
            vec<vec<vec<vec<int>>>> &lines,
            vec<int> &npairs,
            ReadPathVec &paths,
            vecbvec &bases,
            VecPQVec &quals,
            std::string out_dir,
            std::string out_prefix){
    OutputLog(2)<<"DISCOVAR-like heuristics being run"<<std::endl;


    SimplifyDV(out_dir, hbv, hbvinv, paths, bases, quals);

    // For now, fix paths and write the and their inverse
    for (size_t i = 0; i <  paths.size(); i++) { //XXX TODO: change this int for uint 32
        Bool bad = False;
        for (size_t j = 0; j <  paths[i].size(); j++)
            if (paths[i][j] < 0) bad = True;
        if (bad) paths[i].resize(0);
    }
    VecULongVec pathsinv;
    OutputLog(2)<<"creating path-to-edge mapping"<<std::endl;
    invert(paths,pathsinv,hbv.EdgeObjectCount());

    FindLines(hbv, hbvinv, lines, MAX_CELL_PATHS, MAX_DEPTH);
    BinaryWriter::writeFile(out_dir + "/" + out_prefix + ".fin.lines", lines);

    // XXX TODO: Solve the {} thingy, check if has any influence in the new code to run that integrated
    vec<int> llens;
    GetLineLengths(hbv, lines, llens);
    GetLineNpairs(hbv, hbvinv, paths, lines, npairs);
    BinaryWriter::writeFile(out_dir + "/" + out_prefix + ".fin.lines.npairs", npairs);

    vec<vec<covcount>> covs;
    vec<int64_t> subsam_starts={0};
    ComputeCoverage(hbv, hbvinv, paths, lines, subsam_starts, covs);

    //TODO: maybe Report some similar to CN stats ???
    //double cn_frac_good = CNIntegerFraction(hbv, covs);
    //std::cout << "CN fraction good = " << cn_frac_good << std::endl;
    //PerfStatLogger::log("cn_frac_good", ToString(cn_frac_good, 2), "fraction of edges with CN near integer");

}

void step_7(HyperBasevector &hbv,
            vec<int> &hbvinv,
            vec<vec<vec<vec<int>>>> &lines,
            vec<int> &npairs,
            ReadPathVec &paths,
            vecbvec &bases,
            VecPQVec &quals,
            unsigned int min_input_reads,
            std::string out_dir,
            std::string out_prefix){

    int MAX_SUPP_DEL = min_input_reads;//was 0
    bool TAMP_EARLY_MIN = True;
    int MIN_RATIO2 = 8;
    int MAX_DEL2 = 200;
    bool ANALYZE_BRANCHES_VERBOSE2 = False;
    const String TRACE_SEQ = "";
    bool DEGLOOP = True;
    bool EXT_FINAL = True;
    int EXT_FINAL_MODE = 1;
    bool PULL_APART_VERBOSE = False;
    const vec<int> PULL_APART_TRACE;
    int DEGLOOP_MODE = 1;
    float DEGLOOP_MIN_DIST = 2.5;
    bool IMPROVE_PATHS = True;
    bool IMPROVE_PATHS_LARGE = False;
    bool FINAL_TINY = True;
    bool UNWIND3 = True;

    Simplify(out_dir, hbv, hbvinv, paths, bases, quals, MAX_SUPP_DEL, TAMP_EARLY_MIN, MIN_RATIO2, MAX_DEL2,
             ANALYZE_BRANCHES_VERBOSE2, TRACE_SEQ, DEGLOOP, EXT_FINAL, EXT_FINAL_MODE,
             PULL_APART_VERBOSE, PULL_APART_TRACE, DEGLOOP_MODE, DEGLOOP_MIN_DIST, IMPROVE_PATHS,
             IMPROVE_PATHS_LARGE, FINAL_TINY, UNWIND3);//TODO: the last 3 Falses disable pathfinder

    // For now, fix paths and write the and their inverse
    for (size_t i = 0; i < paths.size(); i++) { //XXX TODO: change this int for uint 32
        Bool bad = False;
        for (size_t j = 0; j < paths[i].size(); j++)
            if (paths[i][j] < 0) bad = True;
        if (bad) paths[i].resize(0);
    }
    VecULongVec pathsinv;
    OutputLog(2)<<"creating path-to-edge mapping"<<std::endl;
    invert(paths,pathsinv,hbv.EdgeObjectCount());

    // Find lines and write files.

    FindLines(hbv, hbvinv, lines, MAX_CELL_PATHS, MAX_DEPTH);
    BinaryWriter::writeFile(out_dir + "/" + out_prefix + ".fin.lines", lines);

    // XXX TODO: Solve the {} thingy, check if has any influence in the new code to run that integrated
      vec<int> llens;
      GetLineLengths(hbv, lines, llens);
      GetLineNpairs(hbv, hbvinv, paths, lines, npairs);
      BinaryWriter::writeFile(out_dir + "/" + out_prefix + ".fin.lines.npairs", npairs);

      vec<vec<covcount>> covs;
      vec<int64_t> subsam_starts={0};
      ComputeCoverage(hbv, hbvinv, paths, lines, subsam_starts, covs);

      //TODO: maybe Report some similar to CN stats ???
      //double cn_frac_good = CNIntegerFraction(hbv, covs);
      //std::cout << "CN fraction good = " << cn_frac_good << std::endl;
      //PerfStatLogger::log("cn_frac_good", ToString(cn_frac_good, 2), "fraction of edges with CN near integer");

}

void step_7EXP(HyperBasevector &hbv,
            vec<int> &hbvinv,
            vec<vec<vec<vec<int>>>> &lines,
            vec<int> &npairs,
            ReadPathVec &paths,
            vecbvec &bases,
            VecPQVec &quals,
            unsigned int min_input_reads,
            std::string out_dir,
            std::string out_prefix){
    OutputLog(2)<<"EXPERIMENTAL heuristics being run"<<std::endl;
    int MAX_SUPP_DEL = min_input_reads;//was 0
    bool TAMP_EARLY_MIN = True;
    int MIN_RATIO2 = 8;
    int MAX_DEL2 = 200;
    bool ANALYZE_BRANCHES_VERBOSE2 = False;
    const String TRACE_SEQ = "";
    bool DEGLOOP = True;
    bool EXT_FINAL = True;
    int EXT_FINAL_MODE = 1;
    bool PULL_APART_VERBOSE = False;
    const vec<int> PULL_APART_TRACE;
    int DEGLOOP_MODE = 1;
    float DEGLOOP_MIN_DIST = 2.5;
    bool IMPROVE_PATHS = True;
    bool IMPROVE_PATHS_LARGE = False;
    bool FINAL_TINY = True;
    bool UNWIND3 = True;

    SimplifyEXP(out_dir, hbv, hbvinv, paths, bases, quals, MAX_SUPP_DEL, TAMP_EARLY_MIN, MIN_RATIO2, MAX_DEL2,
                ANALYZE_BRANCHES_VERBOSE2, TRACE_SEQ, DEGLOOP, EXT_FINAL, EXT_FINAL_MODE,
                PULL_APART_VERBOSE, PULL_APART_TRACE, DEGLOOP_MODE, DEGLOOP_MIN_DIST, IMPROVE_PATHS,
                IMPROVE_PATHS_LARGE, FINAL_TINY, UNWIND3, False, False, False);//TODO: the last 3 Falses disable pathfinder

    // For now, fix paths and write the and their inverse
    for (size_t i = 0; i < paths.size(); i++) { //XXX TODO: change this int for uint 32
        Bool bad = False;
        for (size_t j = 0; j < paths[i].size(); j++)
            if (paths[i][j] < 0) bad = True;
        if (bad) paths[i].resize(0);
    }
    VecULongVec pathsinv;
    OutputLog(2)<<"creating path-to-edge mapping"<<std::endl;
    invert(paths,pathsinv,hbv.EdgeObjectCount());

    FindLines(hbv, hbvinv, lines, MAX_CELL_PATHS, MAX_DEPTH);
    BinaryWriter::writeFile(out_dir + "/" + out_prefix + ".fin.lines", lines);

    // XXX TODO: Solve the {} thingy, check if has any influence in the new code to run that integrated
    vec<int> llens;
    GetLineLengths(hbv, lines, llens);
    GetLineNpairs(hbv, hbvinv, paths, lines, npairs);
    BinaryWriter::writeFile(out_dir + "/" + out_prefix + ".fin.lines.npairs", npairs);

    vec<vec<covcount>> covs;
    vec<int64_t> subsam_starts={0};
    ComputeCoverage(hbv, hbvinv, paths, lines, subsam_starts, covs);

    //TODO: maybe Report some similar to CN stats ???
    //double cn_frac_good = CNIntegerFraction(hbv, covs);
    //std::cout << "CN fraction good = " << cn_frac_good << std::endl;
    //PerfStatLogger::log("cn_frac_good", ToString(cn_frac_good, 2), "fraction of edges with CN near integer");

}


struct cmdline_args {
    std::vector<std::string> r1_files;
    std::vector<std::string> r2_files;
    std::string out_dir;
    std::string prefix;
    std::string tmp_dir;
    unsigned int threads=4;
    unsigned int minFreq=4;
    unsigned int minQual=0;
    int max_mem=10;
    uint64_t count_batch_size=0;
    unsigned int small_K=60;
    unsigned int large_K=200;
    unsigned int min_size=0;
    unsigned int from_step=1;
    unsigned int to_step=8;
    unsigned int pair_sample=200;
    unsigned int disk_batches=0;
    unsigned int min_input_reads=3;
    bool paired_repath=true;
    bool solve_complex_repeats=false;
    bool dv_repeat_solving=false;
    bool dump_all=false;
    unsigned int log_level=3;

};


struct cmdline_args parse_cmdline_args( int argc,  char* argv[]) {
    struct cmdline_args parsed_args;
    try {
        cxxopts::Options options(argv[0], "");

        options.add_options()
                ("1", "r1 file(s)", cxxopts::value(parsed_args.r1_files))
                ("2", "r2 file(s)", cxxopts::value(parsed_args.r2_files))
                ("o,output_dir","output dir", cxxopts::value(parsed_args.out_dir))
                ("p,prefix","output prefix", cxxopts::value(parsed_args.prefix))
                ("tmp_dir","temporary files dir (default: output_dir)", cxxopts::value(parsed_args.tmp_dir))
                ("t,threads","parallel threads (default: 4)", cxxopts::value(parsed_args.threads))
                ("max_memory","memory soft limit, in GB (default: 10)",cxxopts::value(parsed_args.max_mem))
                ("min_freq","minimum frequency for k-mers on first DBG",cxxopts::value(parsed_args.minFreq))
                ("min_qual","quality to trim read ends on first DBG (default: 0 - don't trim)",cxxopts::value(parsed_args.minQual))
                ("small_k","k for first DBG",cxxopts::value(parsed_args.small_K))
                ("large_K","k for second DBG",cxxopts::value(parsed_args.large_K))
                ("paired_repath","(EXPERIMENTAL) use pairs when repathing from first to second DBG",cxxopts::value(parsed_args.paired_repath))
                ("dv_repeat_solving","(OBSOLETE) use DISCOVAR's repeat resolution parameters",cxxopts::value(parsed_args.dv_repeat_solving))
                ("solve_complex_repeats","(EXPERIMENTAL) enable full PathFinder heuristics",cxxopts::value(parsed_args.solve_complex_repeats))
                ("from_step","starting step (default: 1)",cxxopts::value(parsed_args.from_step))
                ("to_step","ending step (default: 8)",cxxopts::value(parsed_args.to_step))
                ("dump_all","dump all intermediate graphs and status files",cxxopts::value(parsed_args.dump_all))
                ("h,help","show help message")
                ;

        auto result = options.parse(argc, argv);

        if (result.count("help")) {
            std::cout << options.help({""}) << std::endl;
            exit(0);
        }

        if (result.count("1")!=result.count("2")) {
            std::cout << "Please specify -1 and -2 in pairs of read files" << std::endl;
            exit(0);
        }

        if (result.count("output_dir")!=1 or result.count("prefix")!=1) {
            std::cout << "Please specify output dir and prefix" << std::endl;
            exit(0);
        }

        if (result.count("tmp_dir")==0) {
            parsed_args.tmp_dir=parsed_args.out_dir;
        }


        return parsed_args;

    } catch (const cxxopts::OptionException& e) {
        std::cout << "error parsing options: " << e.what() << std::endl;
        exit(1);
    }
}

std::string step_names[8]={
        "Processing read files",
        "Kmer counting",
        "Small K graph construction",
        "Large K graph construction",
        "Large K graph cleanup",
        "Local Assembly",
        "Final graph cleanup",
        "Paired-end scaffolding"
};
std::string step_outputg_prefix[8]={
        "",
        "",
        "small_K",
        "large_K",
        "large_K_clean",
        "large_K_patched",
        "contigs",
        "assembly"
};

std::string step_inputg_prefix[8]={
        "",
        "",
        "",
        "small_K",
        "large_K",
        "large_K_clean",
        "large_K_patched",
        "contigs"
};

int main( int argc,  char * argv[]) {
    std::cout << "w2rap-contigger - a paired short read assembler" << std::endl << std::endl;
    std::cout << "Git origin: " << GIT_ORIGIN_URL << " -> "  << GIT_BRANCH << std::endl;
    std::cout << "Git commit: " << GIT_COMMIT_HASH << std::endl<<std::endl;
    std::cout << "Command: ";
    for (auto i=0;i<argc;i++) std::cout<<argv[i]<<" ";
    std::cout<<std::endl<<std::endl;
    //Parse and validate args into args structure (DO NOT PASS THE STRUCTURE AROUND!)

    auto args=parse_cmdline_args(argc,argv);
    OutputLogLevel=args.log_level;

    //===== Check directory exists =====

    struct stat info{};
    if (stat(args.out_dir.c_str(), &info) != 0 || !(info.st_mode & S_IFDIR)) {
        std::cout << "Output directory doesn't exist, or is not a directory: " << args.out_dir << std::endl;
        return 1;
    }



#if __GNUC__ > 4 || \
              (__GNUC__ == 4 && (__GNUC_MINOR__ > 9 || \
                                 (__GNUC_MINOR__ == 2 && \
                                  __GNUC_PATCHLEVEL__ > 0)))

    if (omp_get_proc_bind()==omp_proc_bind_false) std::cout<< "WARNING: you are running the code with omp_proc_bind_false, parallel performance may suffer"<<std::endl;
    if (omp_get_proc_bind()==omp_proc_bind_master) std::cout<< "WARNING: you are running the code with omp_proc_bind_master, parallel performance may suffer"<<std::endl;
#endif


    //===== Set computational resources =====
    SetThreads(args.threads, false);
    SetMaxMemory(int64_t(round(args.max_mem * 1024.0 * 1024.0 * 1024.0)));
    //TODO: try to find out max memory on the system to default to.

    //===== Main variables, used across 2 or more steps =====

    vecbvec bases;
    VecPQVec quals;
    std::shared_ptr<KmerList> kmercounts=std::make_shared<KmerList>();
    HyperBasevector hbv;
    vec<int> hbvinv;
    ReadPathVec paths;
    vec<vec<vec<vec<int>>>> lines;
    vec<int> npairs;

    //==== Step-by-step execution loop =====

    for (auto step=args.from_step; step <=args.to_step; ++step){
        int ostep = step-1;

        //===== pre-step: Make sure all data is there =====
        if ( (2==step or 3==step or 5==step or 6==step or 7==step) and (quals.size()==0 or bases.size()==0)){
            if (bases.size()==0) {
                OutputLog(2) << "Loading bases..." << std::endl;
                bases.ReadAll(args.out_dir + "/pe_data.fastb");
            }
            if (quals.size()==0) {
                OutputLog(2) << "Loading quals..." << std::endl;
                load_quals(quals, args.out_dir + "/pe_data.cqual");
            }
            OutputLog(2) << "Read data loaded" << std::endl << std::endl;
        }
        if ( 3==step and 0==kmercounts->size){
            OutputLog(2) << "Loading kmer counts..." << std::endl;
            kmercounts->load(args.out_dir+"/raw_kmers.data");
            OutputLog(2) << "kmer count data loaded" << std::endl << std::endl;
        }
        if ( 4==step ) kmercounts.reset(); //cleanup just in case

        if ( 8==step ) {
            if (lines.empty()) {
                BinaryReader::readFile( args.out_dir + "/" + args.prefix + ".fin.lines", &lines );
            }
            if (npairs.empty()) {
                BinaryReader::readFile( args.out_dir + "/" + args.prefix + ".fin.lines.npairs", &npairs );
            }
        }

        //steps that require a graph
        if (step_inputg_prefix[ostep]!="" and hbv.N()==0) {
            //Load hbv
            OutputLog(2) <<"Loading graph..." << std::endl;
            BinaryReader::readFile(args.out_dir + "/" + args.prefix + "." + step_inputg_prefix[ostep] + ".hbv", &hbv);
            //Create inversion
            OutputLog(4) <<"Creating graph involution..." << std::endl;
            hbvinv.clear();
            hbv.Involution(hbvinv);
            //load paths
            OutputLog(2) <<"Loading paths..." << std::endl;
            LoadReadPathVec(paths,(args.out_dir + "/" + args.prefix + "." + step_inputg_prefix[ostep] + ".paths").c_str());
            //create path inversion
            OutputLog(2) << "Graph and paths loaded" << std::endl << std::endl;
            graph_status(hbv);
            path_status(paths);

        }

        //Print step header and start timers and memory metrics


        OutputLog(1,false) << std::endl << "--== Step " << step << ": " << step_names[ostep] <<" ";
        OutputLog(1)<< " ==--"<< std::endl << std::endl;
        auto step_time=WallClockTime();
        //TODO: reset memory metrics?

        //===== STEP SWITCH =====
        switch (step) {
            //===== STEP 1 (reads to binary) =====
            case 1: {
                OutputLog(2) << "processing reads into bases and quals" << std::endl;
                for (auto i = 0; i < args.r1_files.size(); ++i) {
                    OutputLog(2) << "R1: " << args.r1_files[i] << "  /  R2: " << args.r2_files[i] << std::endl;
                    ExtractPairedReads(args.r1_files[i], args.r2_files[i], bases, quals);
                }
                break;
            }
            //===== STEP 2 (kmer counting) =====
            case 2: {
                if (args.minQual > 0)
                    OutputLog(2) << "Trimming " << quals.size() << " reads at quality >= " << args.minQual << std::endl;
                vec<uint16_t> rlen;
                create_read_lengths(rlen, quals, args.minQual);
                OutputLog(2) << "Unloading quals" << std::endl;
                quals.clear();
                quals.shrink_to_fit();
                kmercounts = buildKMerCount(bases, rlen, args.minFreq, args.out_dir, args.tmp_dir, args.disk_batches,
                                            args.count_batch_size);
                kmercounts->dump(args.out_dir + "/raw_kmers.data");
                OutputLog(2) << "Kmer counting done and dumped with " << kmercounts->size << " kmers" << std::endl;
                break;
            }
            //===== STEP 3 (kmers -> small_k graph) =====
            case 3: {
                bool FILL_JOIN = False;
                buildReadQGraph(bases, quals, kmercounts, FILL_JOIN, FILL_JOIN, args.minFreq, .75, 0, &hbv, &paths,
                                args.small_K);

                kmercounts.reset();
                OutputLog(2) << "computing graph involution and fragment sizes" << std::endl;
                hbvinv.clear();
                hbv.Involution(hbvinv);
                FragDist(hbv, hbvinv, paths, args.out_dir + "/small_K.frags.dist");
                break;
            }
            //===== STEP 4 (small_k graph -> large_k graph) =====
            case 4: {
                //Swap the old graph and such to variables private to this context
                ReadPathVec old_paths;
                std::swap(old_paths, paths);
                paths.resize(old_paths.size());
                HyperBasevector old_hbv;
                std::swap(hbv, old_hbv);
                vec<int> old_hbvinv;
                std::swap(hbvinv, old_hbvinv);

                if (!args.paired_repath) {
                    vecbvec old_edges(old_hbv.Edges().begin(), old_hbv.Edges().end()); //TODO: why do we even need this?
                    //Produce the new graph and such in the argument variables
                    RepathInMemory(old_hbv, old_edges, old_hbvinv, old_paths, old_hbv.K(), args.large_K, hbv, paths,
                                   True,
                                   True, False);
                } else {
                    RepathInMemoryEXP(old_hbv, old_hbvinv, old_paths, args.large_K, hbv, paths);

                }
                hbvinv.clear();
                hbv.Involution(hbvinv);
                FragDist(hbv, hbvinv, paths, args.out_dir + "/large_K.frags.dist");
                break;
            }
            //===== STEP 5 (large_k graph cleaning) =====
            case 5: {
                OutputLog(2) << "cleaning graph" << std::endl;
                int CLEAN_200_VERBOSITY = 0;
                int CLEAN_200V = 3;
                Clean200x(hbv, hbvinv, paths, bases, quals, CLEAN_200_VERBOSITY, CLEAN_200V, args.min_size);
                break;
            }
            //===== STEP 6 (local assemblies) =====
            case 6:{
                vecbvec new_stuff;
                //TODO: Hardcoded parameters
                bool CYCLIC_SAVE = True;
                int A2V = 5;
                int MAX_PROX_LEFT = 400;
                int MAX_PROX_RIGHT = 400;
                int MAX_BPATHS = 100000;
                std::vector<int> k2floor_sequence = {0, 100, 128, 144, 172, 200};
                if (hbv.K() >= 224) k2floor_sequence.push_back(224);
                if (hbv.K() >= 240) k2floor_sequence.push_back(240);
                if (hbv.K() >= 260) k2floor_sequence.push_back(260);
                AssembleGaps2(hbv, hbvinv, paths, bases, quals, args.out_dir, k2floor_sequence,
                              new_stuff, CYCLIC_SAVE, A2V, MAX_PROX_LEFT, MAX_PROX_RIGHT, MAX_BPATHS, args.pair_sample);
                int MIN_GAIN = 5;
                int EXT_MODE = 1;

                AddNewStuff(new_stuff, hbv, hbvinv, paths, bases, quals, MIN_GAIN, EXT_MODE);
                PartnersToEnds(hbv, paths, bases, quals);
                break;
            }
            //===== STEP 7 (repeat resolution, 3 fairly complex approaches, keeping this in separated functions) =====
            case 7: {
                if (args.dv_repeat_solving) step_7DV(hbv, hbvinv, lines, npairs, paths, bases, quals, args.out_dir, args.prefix);
                else if (args.solve_complex_repeats)
                    step_7EXP(hbv, hbvinv, lines, npairs, paths, bases, quals, args.min_input_reads, args.out_dir, args.prefix);
                else step_7(hbv, hbvinv, lines, npairs, paths, bases, quals, args.min_input_reads, args.out_dir, args.prefix);
                break;
            }
            case 8: {
                int MIN_LINE = 5000;
                int MIN_LINK_COUNT = 3; //XXX TODO: this variable is the same as -w in soap??

                bool SCAFFOLD_VERBOSE = False;
                bool GAP_CLEANUP = True;

                VecULongVec pathsinv;
                OutputLog(2)<<"creating path-to-edge mapping"<<std::endl;
                invert(paths,pathsinv,hbv.EdgeObjectCount());

                MakeGaps(hbv, hbvinv, lines, npairs, paths, pathsinv, MIN_LINE, MIN_LINK_COUNT, args.out_dir, args.prefix, SCAFFOLD_VERBOSE, GAP_CLEANUP);

                // Carry out final analyses and write final assembly files.

                vecbasevector G;
                vec<int64_t> subsam_starts={0};
                vec<String> subsam_names={"C"};
                FinalFiles(hbv, hbvinv, paths, subsam_names, subsam_starts, args.out_dir, args.prefix+ ".assembly", MAX_CELL_PATHS, MAX_DEPTH, G);
                break;
            }
            default:
                std::cout << "ERROR: This is not a valid step " << step << std::endl;
                exit(1);
                break;

        }


        //Cleanup


        //TODO: Report time, memory, etc


        if (1==step) {
            //TODO: dump reads
            OutputLog(2) << "Dumping reads in fastb/cqual format..." << std::endl;
            bases.WriteAll(args.out_dir + "/pe_data.fastb");
            save_quals(quals,args.out_dir + "/pe_data.cqual");
            OutputLog(2) << "DONE!" << std::endl;
        } else {
            if (step_outputg_prefix[ostep]!="") {
                if (args.dump_all or step == args.to_step) {
                    //TODO: dump graph and paths
                    OutputLog(2) << "Dumping graph and paths..." << std::endl;
                    BinaryWriter::writeFile(args.out_dir + "/" + args.prefix + "." + step_outputg_prefix[ostep] + ".hbv",
                                            hbv);
                    WriteReadPathVec(paths, (args.out_dir + "/" + args.prefix + "." + step_outputg_prefix[ostep] +
                                             ".paths").c_str());
                    GFADump(std::string(args.out_dir + "/" + args.prefix + "." + step_outputg_prefix[ostep]), hbv, hbvinv,
                            paths, 0, 0, false);
                    OutputLog(2) << "DONE!" << std::endl;
                }
            }

        }
        OutputLog(1) << "Step "<< step << " completed in "<<TimeSince(step_time)<<std::endl<<std::endl;
        if (step_outputg_prefix[ostep]!=""){
            graph_status(hbv);
            path_status(paths);
            OutputLog(2,false)<<std::endl;
        }
    }
    //===== STEP 1: read formating and loading =====



    //===== STEP 1: read formating and loading =====
}

