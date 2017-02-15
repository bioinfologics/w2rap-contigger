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
#include "paths/HyperBasevector.h"
#include "paths/RemodelGapTools.h"
#include "paths/long/BuildReadQGraph.h"
//#include "paths/long/PlaceReads0.h"
#include "paths/long/SupportedHyperBasevector.h"
#include "paths/long/large/AssembleGaps.h"
#include "paths/long/large/ExtractReads.h"
#include "tclap/CmdLine.h"
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
#include <omp.h>

//TODO: stupid globals!

int MAX_CELL_PATHS = 50;
int MAX_DEPTH = 10;


void step_1(vecbvec & bases,
            VecPQVec & quals,
            std::string out_dir,
            std::string read_files){
    OutputLog(2) << "processing reads into bases and quals"<<std::endl;
    ExtractReads(read_files, out_dir, &bases, &quals);
}

void step_2(std::shared_ptr<KmerList> & kmercounts, vecbvec const& reads, VecPQVec &quals,
            unsigned int minQual, unsigned minCount,
            std::string workdir, std::string tmpdir, unsigned char disk_batches, uint64_t count_batch_size) {
    vec<uint16_t> rlen;
    create_read_lengths(rlen,quals,minQual);
    OutputLog(2)<<"Unloading quals"<<std::endl;
    quals.clear();
    quals.shrink_to_fit();
    kmercounts=buildKMerCount( reads, rlen, minCount, workdir, tmpdir, disk_batches, count_batch_size );
    kmercounts->dump(workdir+"/raw_kmers.data");
    OutputLog(2) << "Kmer counting done and dumped with "<<kmercounts->size<< " kmers" <<std::endl;
}

void step_3(HyperBasevector &hbv,
            vec<int> &hbvinv,
            ReadPathVec &paths,
            vecbvec &bases,
            VecPQVec &quals,
            std::shared_ptr<KmerList> kmercounts,
            unsigned int minFreq,
            unsigned int small_K,
            std::string out_dir) {
    bool FILL_JOIN = False;

    buildReadQGraph(bases, quals, kmercounts, FILL_JOIN, FILL_JOIN, minFreq, .75, 0, &hbv, &paths, small_K);

    kmercounts.reset();
    OutputLog(2)<<"computing graph involution and fragment sizes"<<std::endl;
    hbvinv.clear();
    hbv.Involution(hbvinv);
    FragDist(hbv, hbvinv, paths, out_dir + "/small_K.frags.dist");
}

void step_4DV(HyperBasevector &hbv,
            vec<int> &hbvinv,
            ReadPathVec &paths,
            unsigned int large_K,
            std::string out_dir){

    //Swap the old graph and such to variables private to this context
    ReadPathVec old_paths;
    std::swap(old_paths,paths);
    paths.resize(old_paths.size());
    HyperBasevector old_hbv;
    std::swap(hbv,old_hbv);
    vec<int> old_hbvinv;
    std::swap(hbvinv,old_hbvinv);

    vecbvec old_edges(old_hbv.Edges().begin(), old_hbv.Edges().end()); //TODO: why do we even need this?

    //Produce the new graph and such in the argument variables
    RepathInMemory(old_hbv, old_edges, old_hbvinv, old_paths, old_hbv.K(), large_K, hbv, paths, True, True, False);
    OutputLog(2)<<"computing graph involution and fragment sizes"<<std::endl;
    hbvinv.clear();
    hbv.Involution(hbvinv);
    FragDist(hbv, hbvinv, paths, out_dir + "/large_K.frags.dist");
}

void step_4(HyperBasevector &hbv,
            vec<int> &hbvinv,
            ReadPathVec &paths,
            unsigned int large_K,
            std::string out_dir){
    step_4DV(hbv,hbvinv,paths,large_K,out_dir);
}



void step_4EXP(HyperBasevector &hbv,
               vec<int> &hbvinv,
               ReadPathVec &paths,
               unsigned int large_K,
               std::string out_dir){
    OutputLog(2)<<"EXPERIMENTAL heuristics being run"<<std::endl;
    /*if (clean_smallk_graph) {
        std::cout << "--== Step 3a: improving small_K graph ==--" << std::endl;
        inv.clear();
        hbv.Involution(inv);
        VecULongVec paths_inv;
        invert(paths, paths_inv, hbv.EdgeObjectCount());
        GraphImprover gi(hbv, inv, paths, paths_inv);
        gi.improve_graph();
        GraphImprover gi2(hbv, inv, paths, paths_inv);
        gi2.expand_cannonical_repeats(2,2);


        std::cout << "--== Step 3b: Repathing to second (large K) graph ==--" << std::endl;
    } else {
        std::cout << "--== Step 3: Repathing to second (large K) graph ==--" << std::endl;
    }*/
    ReadPathVec old_paths;
    std::swap(old_paths,paths);
    paths.resize(old_paths.size());
    HyperBasevector old_hbv;
    std::swap(hbv,old_hbv);
    vec<int> old_hbvinv;
    std::swap(hbvinv,old_hbvinv);



    //Produce the new graph and such in the argument variables
    RepathInMemoryEXP(old_hbv, old_hbvinv, old_paths, large_K, hbv, paths);
    OutputLog(2)<<"computing graph involution and fragment sizes"<<std::endl;
    hbvinv.clear();
    hbv.Involution(hbvinv);
    FragDist(hbv, hbvinv, paths, out_dir + "/large_K.frags.dist");

}

void step_5(HyperBasevector &hbv,
            vec<int> &hbvinv,
            ReadPathVec &paths,
            vecbvec &bases,
            VecPQVec &quals,
            unsigned int min_size){
    OutputLog(2)<<"cleaning graph"<<std::endl;
    int CLEAN_200_VERBOSITY = 0;
    int CLEAN_200V = 3;
    Clean200x(hbv, hbvinv, paths, bases, quals, CLEAN_200_VERBOSITY, CLEAN_200V, min_size);

}

void step_6(HyperBasevector &hbv,
            vec<int> &hbvinv,
            ReadPathVec &paths,
            vecbvec &bases,
            VecPQVec &quals,
            unsigned int pair_sample,
            std::string out_dir){
    vecbvec new_stuff;
    //TODO: Hardcoded parameters
    bool CYCLIC_SAVE = True;
    int A2V = 5;
    int MAX_PROX_LEFT = 400;
    int MAX_PROX_RIGHT = 400;
    int MAX_BPATHS = 100000;
    std::vector<int> k2floor_sequence={0, 100, 128, 144, 172, 200};
    if (hbv.K()>=224) k2floor_sequence.push_back(224);
    if (hbv.K()>=240) k2floor_sequence.push_back(240);
    if (hbv.K()>=260) k2floor_sequence.push_back(260);
    AssembleGaps2(hbv, hbvinv, paths, bases, quals, out_dir, k2floor_sequence,
                  new_stuff, CYCLIC_SAVE, A2V, MAX_PROX_LEFT, MAX_PROX_RIGHT, MAX_BPATHS, pair_sample);
    int MIN_GAIN = 5;
    int EXT_MODE = 1;

    AddNewStuff(new_stuff, hbv, hbvinv, paths, bases, quals, MIN_GAIN, EXT_MODE);
    PartnersToEnds(hbv, paths, bases, quals);

}

void step_7DV(HyperBasevector &hbv,
            vec<int> &hbvinv,
            ReadPathVec &paths,
            vecbvec &bases,
            VecPQVec &quals,
            std::string out_dir,
            std::string out_prefix){
    OutputLog(2)<<"DISCOVAR-like heuristics being run"<<std::endl;


    SimplifyDV(out_dir, hbv, hbvinv, paths, bases, quals);

    // For now, fix paths and write the and their inverse
    for (int i = 0; i < (int) paths.size(); i++) { //XXX TODO: change this int for uint 32
        Bool bad = False;
        for (int j = 0; j < (int) paths[i].size(); j++)
            if (paths[i][j] < 0) bad = True;
        if (bad) paths[i].resize(0);
    }
    VecULongVec pathsinv;
    OutputLog(2)<<"creating path-to-edge mapping"<<std::endl;
    invert(paths,pathsinv,hbv.EdgeObjectCount());

    // Find lines and write files.
    vec<vec<vec<vec<int>>>> lines;

    FindLines(hbv, hbvinv, lines, MAX_CELL_PATHS, MAX_DEPTH);
    BinaryWriter::writeFile(out_dir + "/" + out_prefix + ".fin.lines", lines);

    // XXX TODO: Solve the {} thingy, check if has any influence in the new code to run that integrated
    {
        vec<int> llens, npairs;
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

}

void step_7(HyperBasevector &hbv,
            vec<int> &hbvinv,
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
    for (int i = 0; i < (int) paths.size(); i++) { //XXX TODO: change this int for uint 32
        Bool bad = False;
        for (int j = 0; j < (int) paths[i].size(); j++)
            if (paths[i][j] < 0) bad = True;
        if (bad) paths[i].resize(0);
    }
    VecULongVec pathsinv;
    OutputLog(2)<<"creating path-to-edge mapping"<<std::endl;
    invert(paths,pathsinv,hbv.EdgeObjectCount());

    // Find lines and write files.
    vec<vec<vec<vec<int>>>> lines;

    FindLines(hbv, hbvinv, lines, MAX_CELL_PATHS, MAX_DEPTH);
    BinaryWriter::writeFile(out_dir + "/" + out_prefix + ".fin.lines", lines);

    // XXX TODO: Solve the {} thingy, check if has any influence in the new code to run that integrated
    {
        vec<int> llens, npairs;
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

}

void step_7EXP(HyperBasevector &hbv,
            vec<int> &hbvinv,
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
    for (int i = 0; i < (int) paths.size(); i++) { //XXX TODO: change this int for uint 32
        Bool bad = False;
        for (int j = 0; j < (int) paths[i].size(); j++)
            if (paths[i][j] < 0) bad = True;
        if (bad) paths[i].resize(0);
    }
    VecULongVec pathsinv;
    OutputLog(2)<<"creating path-to-edge mapping"<<std::endl;
    invert(paths,pathsinv,hbv.EdgeObjectCount());

    // Find lines and write files.
    vec<vec<vec<vec<int>>>> lines;

    FindLines(hbv, hbvinv, lines, MAX_CELL_PATHS, MAX_DEPTH);
    BinaryWriter::writeFile(out_dir + "/" + out_prefix + ".fin.lines", lines);

    // XXX TODO: Solve the {} thingy, check if has any influence in the new code to run that integrated
    {
        vec<int> llens, npairs;
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

}

void step_8(HyperBasevector &hbv,
            vec<int> &hbvinv,
            ReadPathVec &paths,
            std::string out_dir,
            std::string out_prefix){
    int MIN_LINE = 5000;
    int MIN_LINK_COUNT = 3; //XXX TODO: this variable is the same as -w in soap??

    bool SCAFFOLD_VERBOSE = False;
    bool GAP_CLEANUP = True;

    VecULongVec pathsinv;
    OutputLog(2)<<"creating path-to-edge mapping"<<std::endl;
    invert(paths,pathsinv,hbv.EdgeObjectCount());

    MakeGaps(hbv, hbvinv, paths, pathsinv, MIN_LINE, MIN_LINK_COUNT, out_dir, out_prefix, SCAFFOLD_VERBOSE, GAP_CLEANUP);

    // Carry out final analyses and write final assembly files.

    vecbasevector G;
    vec<int64_t> subsam_starts={0};
    vec<String> subsam_names={"C"};
    FinalFiles(hbv, hbvinv, paths, subsam_names, subsam_starts, out_dir, out_prefix+ "_assembly", MAX_CELL_PATHS, MAX_DEPTH, G);
    GFADump(out_prefix+ "_assembly",hbv,hbvinv,paths,MAX_CELL_PATHS,MAX_DEPTH,true);
}



int main(const int argc, const char * argv[]) {

    std::string out_prefix;
    std::string read_files;
    std::string out_dir;
    std::string tmp_dir;
    unsigned int threads;
    unsigned int minFreq,minCount;
    unsigned int minQual;
    int max_mem;
    uint64_t count_batch_size;
    unsigned int small_K, large_K, min_size,from_step,to_step, pair_sample, disk_batches, min_input_reads;
    std::vector<unsigned int> allowed_k = {60, 64, 72, 80, 84, 88, 96, 100, 108, 116, 128, 136, 144, 152, 160, 168, 172,
                                           180, 188, 192, 196, 200, 208, 216, 224, 232, 240, 260, 280, 300, 320, 368,
                                           400, 440, 460, 500, 544, 640};
    std::vector<unsigned int> allowed_steps = {1,2,3,4,5,6,7,8};
    bool dump_detailed_gfa,dump_all,run_dv,run_exp;

    //========== Command Line Option Parsing ==========
    for (auto i=0;i<argc;i++) std::cout<<argv[i]<<" ";
    std::cout<<std::endl<<std::endl;
    std::cout << "Welcome to w2rap-contigger" << std::endl << std::endl;

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
        //TCLAP::ValueArg<unsigned int> small_KArg("k", "small_k",
        //                                         "Small k (default: 60)", false, 60, &largeKconst, cmd);

        TCLAP::ValuesConstraint<unsigned int> steps(allowed_steps);
        TCLAP::ValueArg<unsigned int> fromStep_Arg("", "from_step",
                                                 "Start on step (default: 1)", false, 1, &steps, cmd);

        TCLAP::ValueArg<unsigned int> toStep_Arg("", "to_step",
                                                   "Stop after step (default: 8)", false, 8, &steps, cmd);

        TCLAP::ValueArg<unsigned int> disk_batchesArg("d", "disk_batches",
                                                 "number of disk batches for step2 (default: 0, 0->in memory)", false, 0, "int", cmd);

        TCLAP::ValueArg<uint32_t> count_batchArg("", "count_batch",
                                                     "number of reads to count on each parallel iteration (default: total/(4*threads) )", false, 0, "int", cmd);

        TCLAP::ValueArg<std::string> tmp_dirArg("", "tmp_dir",
                                                      "tmp dir for step2 disk batches (default: workdir)", false, "", "string", cmd);

        TCLAP::ValueArg<unsigned int> minQualArg("", "min_qual",
                                                 "minimum quality for small k-mers (default: 7)", false, 7, "int", cmd);

        TCLAP::ValueArg<unsigned int> minCountArg("", "min_count",
                                                 "minimum frequency for k-mers counting (default: 4)", false, 4, "int", cmd);

        TCLAP::ValueArg<unsigned int> minFreqArg("", "min_freq",
                                                 "minimum frequency for small k-mer graph (default: 4)", false, 4, "int", cmd);

        TCLAP::ValueArg<unsigned int> minSizeArg("s", "min_size",
             "Min size of disconnected elements on large_k graph (in kmers, default: 0=no min)", false, 0, "int", cmd);

        TCLAP::ValueArg<unsigned int> pairSampleArg("", "pair_sample",
                                                    "max number of read pairs to use in local assemblies (default: 200)", false, 200, "int", cmd);

        TCLAP::ValueArg<unsigned int> minInputArg("", "min_input",
                                                    "min number of read entering an edge at step 6(default: 3)", false, 3, "int", cmd);

        TCLAP::ValueArg<unsigned int> logLevelArg("", "log_level",
                                                  "verbosity level (default: 3)", false, 3, "1-4", cmd);

        TCLAP::ValueArg<bool>         runDiscovarCompatibleArg        ("","dv_like",
                                                               "Run with discovar-like heuristics (default: 0)", false,false,"bool",cmd);

        TCLAP::ValueArg<bool>         runExperimentalArg        ("","experimental",
                                                                    "Run latest EXPERIMENTAL heuristics (default: 0)", false,false,"bool",cmd);

        TCLAP::ValueArg<bool>         dumpAllArg        ("","dump_all",
                                                               "Dump all intermediate files (default: 0)", false,false,"bool",cmd);
        TCLAP::ValueArg<bool>         dumpDetailedGFAArg        ("","dump_detailed_gfa",
                                                         "Dump detailed GFA for every graph (default: 0)", false,false,"bool",cmd);

        cmd.parse(argc, argv);
        // Get the value parsed by each arg.
        out_dir = out_dirArg.getValue();
        out_prefix = out_prefixArg.getValue();
        read_files = read_filesArg.getValue();
        threads = threadsArg.getValue();
        max_mem = max_memArg.getValue();
        large_K = large_KArg.getValue();
        small_K = 60;//small_KArg.getValue();
        min_size = minSizeArg.getValue();
        dump_detailed_gfa=dumpDetailedGFAArg.getValue();
        dump_all=dumpAllArg.getValue();
        from_step=fromStep_Arg.getValue();
        to_step=toStep_Arg.getValue();
        pair_sample=pairSampleArg.getValue();
        minFreq=minFreqArg.getValue();
        minCount=minCountArg.getValue();
        minQual=minQualArg.getValue();
        disk_batches=disk_batchesArg.getValue();
        count_batch_size=count_batchArg.getValue();
        tmp_dir=tmp_dirArg.getValue();
        min_input_reads=minInputArg.getValue();
        OutputLogLevel=logLevelArg.getValue();
        run_dv=runDiscovarCompatibleArg.getValue();
        run_exp=runExperimentalArg.getValue();

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

#if __GNUC__ > 4 || \
              (__GNUC__ == 4 && (__GNUC_MINOR__ > 9 || \
                                 (__GNUC_MINOR__ == 2 && \
                                  __GNUC_PATCHLEVEL__ > 0)))

    if (omp_get_proc_bind()==omp_proc_bind_false) std::cout<< "WARNING: you are running the code with omp_proc_bind_false, parallel performance may suffer"<<std::endl;
    if (omp_get_proc_bind()==omp_proc_bind_master) std::cout<< "WARNING: you are running the code with omp_proc_bind_master, parallel performance may suffer"<<std::endl;
#endif


    //== Set computational resources ===
    SetThreads(threads, False);
    SetMaxMemory(int64_t(round(max_mem * 1024.0 * 1024.0 * 1024.0)));
    //TODO: try to find out max memory on the system to default to.
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
            "large_K_expanded",
            "large_K_final"
    };

    std::string step_inputg_prefix[8]={
            "",
            "",
            "",
            "small_K",
            "large_K",
            "large_K_clean",
            "large_K_patched",
            "large_K_expanded"
    };

    //========== Main Program Begins ======
    vecbvec bases;
    VecPQVec quals;


    //double wtimer,cputimer;

    HyperBasevector hbv;
    vec<int> hbvinv;
    ReadPathVec paths;
    std::shared_ptr<KmerList> kmercounts=std::make_shared<KmerList>();



    //Step-by-step execution loop
    for (auto step=from_step; step <=to_step; ++step){
        //First make sure all needed data is there.
        if ( (2==step or 3==step or 5==step or 6==step or 7==step) and (quals.size()==0 or bases.size()==0)){
            if (bases.size()==0) {
                OutputLog(2) << "Loading bases..." << std::endl;
                bases.ReadAll(out_dir + "/pe_data.fastb");
            }
            if (quals.size()==0) {
                OutputLog(2) << "Loading quals..." << std::endl;
                load_quals(quals, out_dir + "/pe_data.cqual");
            }
            OutputLog(2) << "Read data loaded" << std::endl << std::endl;
        }
        if ( 3==step and 0==kmercounts->size){
            OutputLog(2) << "Loading kmer counts..." << std::endl;
            kmercounts->load(out_dir+"/raw_kmers.data");
            OutputLog(2) << "Kmer count data loaded" << std::endl << std::endl;
        }
        if ( 4==step ) kmercounts.reset(); //cleanup just in case

        //steps that require a graph
        if (step_inputg_prefix[step-1]!="" and hbv.N()==0) {
            //Load hbv
            OutputLog(2) <<"Loading graph..." << std::endl;
            BinaryReader::readFile(out_dir + "/" + out_prefix + "." + step_inputg_prefix[step-1] + ".hbv", &hbv);
            //Create inversion
            OutputLog(4) <<"Creating graph involution..." << std::endl;
            hbvinv.clear();
            hbv.Involution(hbvinv);
            //load paths
            OutputLog(2) <<"Loading paths..." << std::endl;
            LoadReadPathVec(paths,(out_dir + "/" + out_prefix + "." + step_inputg_prefix[step-1] + ".paths").c_str());
            //create path inversion
            OutputLog(2) << "Graph and paths loaded" << std::endl << std::endl;
            graph_status(hbv);
            path_status(paths);

        }

        //Print step header and start timers and memory metrics


        OutputLog(1,false) << std::endl << "--== Step " << step << ": " << step_names[step-1] <<" ";
        OutputLog(1)<< " ==--"<< std::endl << std::endl;
        auto step_time=WallClockTime();
        //TODO: reset memory metrics?

        //Run step
        switch (step){
            case 1:
                step_1(bases, quals, out_dir, read_files);
                break;
            case 2:
                step_2(kmercounts,bases, quals, minQual, minCount, out_dir, tmp_dir,disk_batches,count_batch_size);
                break;
            case 3:
                step_3(hbv,hbvinv,paths,bases,quals,kmercounts,minFreq,small_K,out_dir);
                break;
            case 4:
                if (run_dv) step_4DV(hbv,hbvinv,paths,large_K,out_dir);
                else if (run_exp) step_4EXP(hbv,hbvinv,paths,large_K,out_dir);
                else step_4(hbv,hbvinv,paths,large_K,out_dir);
                break;
            case 5:
                step_5(hbv, hbvinv, paths, bases, quals, min_size);
                break;
            case 6:
                step_6(hbv, hbvinv, paths, bases, quals, pair_sample, out_dir);
                break;
            case 7:
                if (run_dv) step_7DV(hbv, hbvinv, paths, bases, quals, out_dir, out_prefix);
                else if (run_exp) step_7EXP(hbv, hbvinv, paths, bases, quals, min_input_reads, out_dir, out_prefix);
                else step_7(hbv, hbvinv, paths, bases, quals, min_input_reads, out_dir, out_prefix);
                break;
            case 8:
                step_8(hbv, hbvinv, paths, out_dir, out_prefix);
                break;
        }


        //Cleanup


        //TODO: Report time, memory, etc


        if (1==step and (to_step<6 or dump_all)) {
            //TODO: dump reads
            OutputLog(2) << "Dumping reads in fastb/cqual format..." << std::endl;
            bases.WriteAll(out_dir + "/pe_data.fastb");
            save_quals(quals,out_dir + "/pe_data.cqual");
            OutputLog(2) << "DONE!" << std::endl;
        } else {
            if (step_outputg_prefix[step-1]!="" and (dump_all or step==to_step)){
                //TODO: dump graph and paths
                OutputLog(2) << "Dumping graph and paths..." << std::endl;
                BinaryWriter::writeFile(out_dir + "/" + out_prefix + "." + step_outputg_prefix[step-1] +".hbv", hbv);
                WriteReadPathVec(paths,(out_dir + "/" + out_prefix + "." + step_outputg_prefix[step-1] +".paths").c_str());
                OutputLog(2) << "DONE!" << std::endl;
            }
            if (dump_detailed_gfa){
                //TODO:
            }
        }
        OutputLog(1) << "Step "<< step << " completed in "<<TimeSince(step_time)<<std::endl<<std::endl;
        if (step_outputg_prefix[step-1]!=""){
            graph_status(hbv);
            path_status(paths);
            OutputLog(2,false)<<std::endl;
        }
    }
    return 0;
}

