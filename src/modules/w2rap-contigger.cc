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
#include "GFADump.h"

int main(const int argc, const char * argv[]) {

    std::string out_prefix;
    std::string read_files;
    std::string out_dir;
    std::string dev_run;
    std::string tmp_dir;
    unsigned int threads;
    unsigned int minFreq;
    unsigned int minQual;
    int max_mem;
    unsigned int small_K, large_K, min_size,from_step,to_step, pair_sample, disk_batches, min_input_reads;
    std::vector<unsigned int> allowed_k = {60, 64, 72, 80, 84, 88, 96, 100, 108, 116, 128, 136, 144, 152, 160, 168, 172,
                                           180, 188, 192, 196, 200, 208, 216, 224, 232, 240, 260, 280, 300, 320, 368,
                                           400, 440, 460, 500, 544, 640};
    std::vector<unsigned int> allowed_steps = {1,2,3,4,5,6,7};
    bool extend_paths,run_pathfinder,dump_detailed_gfa,dump_all,dump_pf,pf_verbose,clean_smallk_graph;

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
                                                   "Stop after step (default: 7)", false, 7, &steps, cmd);

        TCLAP::ValueArg<unsigned int> disk_batchesArg("d", "disk_batches",
                                                 "number of disk batches for step2 (default: 0, 0->in memory)", false, 8, "int", cmd);
        TCLAP::ValueArg<std::string> tmp_dirArg("", "tmp_dir",
                                                      "tmp dir for step2 disk batches (default: workdir)", false, "", "string", cmd);
        TCLAP::ValueArg<unsigned int> minSizeArg("s", "min_size",
             "Min size of disconnected elements on large_k graph (in kmers, default: 0=no min)", false, 0, "int", cmd);
        TCLAP::ValueArg<unsigned int> minFreqArg("", "min_freq",
                                                 "minimum frequency for small k-mers on step 2 (default: 4)", false, 4, "int", cmd);
        TCLAP::ValueArg<unsigned int> minQualArg("", "min_qual",
                                                 "minimum quality for small k-mers on step 2 (default: 7)", false, 7, "int", cmd);
        TCLAP::ValueArg<unsigned int> pairSampleArg("", "pair_sample",
                                                    "max number of read pairs to use in local assemblies on step 5(default: 200)", false, 200, "int", cmd);
            TCLAP::ValueArg<unsigned int> minInputArg("", "min_input",
                                                    "min number of read entering an edge at step 6(default: 3)", false, 3, "int", cmd);
        TCLAP::ValueArg<bool>         pathExtensionArg        ("","extend_paths",
                                                               "Enable extend paths on repath (experimental)", false,false,"bool",cmd);
        TCLAP::ValueArg<bool>         cleanSmallKGraphArg        ("","clean_smallk_graph",
                                                               "Enable cleaning on step3 (highly experimental)", false,false,"bool",cmd);
        TCLAP::ValueArg<bool>         pathFinderArg        ("","path_finder",
                                                            "Run PathFinder (experimental)", false,false,"bool",cmd);
        TCLAP::ValueArg<bool>         pathFinderVerboseArg        ("","pf_verbose",
                                                            "PathFinder verbose (experimental)", false,false,"bool",cmd);
        TCLAP::ValueArg<bool>         dumpAllArg        ("","dump_all",
                                                               "Dump all intermediate files", false,false,"bool",cmd);
        TCLAP::ValueArg<bool>         dumpDetailedGFAArg        ("","dump_detailed_gfa",
                                                         "Dump all intermediate files", false,false,"bool",cmd);
        TCLAP::ValueArg<bool>         dumpPFArg        ("","dump_pf",
                                                          "Dump pathfinder info (devel)", false,false,"bool",cmd);

        TCLAP::ValueArg<std::string> dev_runArg("", "dev_run_test",
                                                   "runs development tests", false, "", "devel only", cmd);

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
        extend_paths=pathExtensionArg.getValue();
        run_pathfinder=pathFinderArg.getValue();
        dump_detailed_gfa=dumpDetailedGFAArg.getValue();
        dump_all=dumpAllArg.getValue();
        from_step=fromStep_Arg.getValue();
        to_step=toStep_Arg.getValue();
        dev_run=dev_runArg.getValue();
        dump_pf=dumpPFArg.getValue();
        pair_sample=pairSampleArg.getValue();
        minFreq=minFreqArg.getValue();
        minQual=minQualArg.getValue();
        disk_batches=disk_batchesArg.getValue();
        tmp_dir=tmp_dirArg.getValue();
        pf_verbose=pathFinderVerboseArg.getValue();
        min_input_reads=minInputArg.getValue();
        clean_smallk_graph=cleanSmallKGraphArg.getValue();

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

    if (omp_get_proc_bind()==omp_proc_bind_false) std::cout<< "WARNING: you are running the code with omp_proc_bind_false, parallel performance may suffer"<<std::endl;
    if (omp_get_proc_bind()==omp_proc_bind_master) std::cout<< "WARNING: you are running the code with omp_proc_bind_master, parallel performance may suffer"<<std::endl;

    //========== Main Program Begins ======
    vecbvec bases;
    VecPQVec quals;

    vec<String> subsam_names = {"C"};
    vec<int64_t> subsam_starts = {0};
    //double wtimer,cputimer;
    vec<int> inv;
    HyperBasevector hbvr;
    ReadPathVec pathsr;


    int MAX_CELL_PATHS = 50;
    int MAX_DEPTH = 10;

    //== Set computational resources ===
    SetThreads(threads, False);
    SetMaxMemory(int64_t(round(max_mem * 1024.0 * 1024.0 * 1024.0)));
    //TODO: try to find out max memory on the system to default to.

    //== Load reads (and saves in binary format) ======

    if (from_step==1)
    {
        std::cout << "--== Step 1: Reading input files ==--" << std::endl;
        ExtractReads(read_files, out_dir, subsam_names, subsam_starts, &bases, &quals);
        //TODO: add an option to dump the reads
        if (dump_all || to_step<6) {
            std::cout << Date() << ": Dumping reads in fastb/qualp format..." << std::endl;
            bases.WriteAll(out_dir + "/frag_reads_orig.fastb");
            quals.WriteAll(out_dir + "/frag_reads_orig.qualp");
        }
        std::cout << Date() << ": Reading input files DONE!" << std::endl << std::endl << std::endl;
    }

    //== Read QGraph, and repath (k=60, k=200 (and saves in binary format) ======

    if (from_step>1 && from_step<7 and not (from_step==3 and to_step==3)){
        std::cout << Date() << ": Loading reads in fastb/qualp format..." << std::endl;
        bases.ReadAll(out_dir + "/frag_reads_orig.fastb");
        quals.ReadAll(out_dir + "/frag_reads_orig.qualp");
        std::cout << Date() << ": Loading reads DONE!" << std::endl << std::endl;
    }
    {//This scope-trick to invalidate old data is dirty

        HyperBasevector hbv;
        ReadPathVec paths;
        if (from_step<=2 and to_step>=2) {
            bool FILL_JOIN = False;
            std::cout << "--== Step 2: Building first (small K) graph ==--" << std::endl;
            buildReadQGraph(bases, quals, FILL_JOIN, FILL_JOIN, minQual, minFreq, .75, 0, &hbv, &paths, small_K, out_dir,tmp_dir,disk_batches);
            //FixPaths(hbv, paths); //TODO: is this even needed?
            if(dump_detailed_gfa) GFADumpDetail(out_dir + "/" + out_prefix + ".small_K",hbv,inv);
            if (dump_all || to_step ==2){
                std::cout << Date() << ": Dumping small_K graph and paths..." << std::endl;
                BinaryWriter::writeFile(out_dir + "/" + out_prefix + ".small_K.hbv", hbv);
                WriteReadPathVec(paths,(out_dir + "/" + out_prefix + ".small_K.paths").c_str());
                graph_status(hbv);
                path_status(paths);

                std::cout << Date() << ": Dumping small_K graph and paths DONE!" << std::endl;
            }
            std::cout << Date() << ": Building first graph DONE!" << std::endl << std::endl << std::endl;
        }

        if (from_step==3){
            std::cout << Date() << ": Reading small_K graph and paths..." << std::endl;
            BinaryReader::readFile(out_dir + "/" + out_prefix + ".small_K.hbv", &hbv);
            LoadReadPathVec(paths,(out_dir + "/" + out_prefix + ".small_K.paths").c_str());
            graph_status(hbv);
            path_status(paths);
            std::cout << Date() << ": Reading small_K graph and paths DONE!" << std::endl << std::endl;
        }
        if (from_step<=3 and to_step>=3) {
            if (clean_smallk_graph) {
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
            }
            vecbvec edges(hbv.Edges().begin(), hbv.Edges().end());
            inv.clear();
            hbv.Involution(inv);
            FragDist(hbv, inv, paths, out_dir + "/" + out_prefix + ".first.frags.dist");
            const string run_head = out_dir + "/" + out_prefix;

            pathsr.resize(paths.size());

            RepathInMemory(hbv, edges, inv, paths, hbv.K(), large_K, hbvr, pathsr, True, True, extend_paths);
            if(dump_detailed_gfa) GFADumpDetail(out_dir + "/" + out_prefix + ".large_K",hbvr,inv);
            if (dump_all || to_step ==3){
                std::cout << Date() << ": Dumping large_K graph and paths..." << std::endl;
                BinaryWriter::writeFile(out_dir + "/" + out_prefix + ".large_K.hbv", hbvr);
                WriteReadPathVec(pathsr,(out_dir + "/" + out_prefix + ".large_K.paths").c_str());
                graph_status(hbv);
                path_status(paths);
                std::cout << Date() << ": Dumping large_K graph and paths DONE!" << std::endl;
            }
            std::cout << Date() << ": Repathing to second graph DONE!" << std::endl << std::endl;
        }

    }

    //== Clean ======
    if (from_step==4){
        std::cout << Date() << ": Reading large_K graph and paths..." << std::endl;
        BinaryReader::readFile(out_dir + "/" + out_prefix + ".large_K.hbv", &hbvr);
        LoadReadPathVec(pathsr,(out_dir + "/" + out_prefix + ".large_K.paths").c_str());
        
        graph_status(hbvr);
        path_status(pathsr);
        std::cout << Date() << ": Reading large_K graph and paths DONE!" << std::endl << std::endl;
    }
    if (from_step<=4 and to_step>=4) {
        std::cout << "--== Step 4: Cleaning graph ==--" << std::endl;
        inv.clear();
        hbvr.Involution(inv);
        int CLEAN_200_VERBOSITY = 0;
        int CLEAN_200V = 3;
        Clean200x(hbvr, inv, pathsr, bases, quals, CLEAN_200_VERBOSITY, CLEAN_200V, min_size);
        if(dump_detailed_gfa) GFADumpDetail(out_dir + "/" + out_prefix + ".large_K.clean",hbvr,inv);
        if (dump_all || to_step ==4){
            std::cout << Date() << ": Dumping large_K clean graph and paths..." << std::endl;
            BinaryWriter::writeFile(out_dir + "/" + out_prefix + ".large_K.clean.hbv", hbvr);
            WriteReadPathVec(pathsr,(out_dir + "/" + out_prefix + ".large_K.clean.paths").c_str());
            graph_status(hbvr);
            path_status(pathsr);
            std::cout << Date() << ": Dumping large_K clean graph and paths DONE!" << std::endl;
        }
        std::cout << Date() << ": Cleaning graph DONE!" << std::endl<< std::endl;
    }

    //== Patching ======

    VecULongVec paths_inv;

    if (from_step==5){
        std::cout << Date() << ": Reading large_K clean graph and paths..." << std::endl;
        BinaryReader::readFile(out_dir + "/" + out_prefix + ".large_K.clean.hbv", &hbvr);
        LoadReadPathVec(pathsr,(out_dir + "/" + out_prefix + ".large_K.clean.paths").c_str());
        inv.clear();
        hbvr.Involution(inv);
        graph_status(hbvr);
        path_status(pathsr);
        std::cout << Date() << ": Reading large_K clean graph and paths DONE!" << std::endl << std::endl;
    }
    if (from_step<=5 and to_step>=5) {
        std::cout << "--== Step 5: Assembling gaps ==--" << std::endl;
        std::cout << Date() <<": creating edge-to-path index"<<std::endl;
        invert(pathsr, paths_inv, hbvr.EdgeObjectCount());
       /* ConsensusChecker cc(hbvr,inv,pathsr,bases);
        for (auto e=0;e<hbvr.EdgeObjectCount();++e){
            if (not cc.consensus_OK(e)){
                std::cout<<"=========================================================================="<<e<<std::endl;
                std::cout<<"Consensus problem on edge "<<e<<std::endl;
                cc.print_detail();
                std::cout<<"=========================================================================="<<e<<std::endl;

            }
        }*/

        vecbvec new_stuff;
        //TODO: Hardcoded parameters
        bool CYCLIC_SAVE = True;
        int A2V = 5;
        int MAX_PROX_LEFT = 400;
        int MAX_PROX_RIGHT = 400;
        int MAX_BPATHS = 100000;
        std::vector<int> k2floor_sequence={0, 100, 128, 144, 172, 200};
        if (hbvr.K()>=224) k2floor_sequence.push_back(224);
        if (hbvr.K()>=240) k2floor_sequence.push_back(240);
        if (hbvr.K()>=260) k2floor_sequence.push_back(260);

        AssembleGaps2(hbvr, inv, pathsr, paths_inv, bases, quals, out_dir, k2floor_sequence,
                      new_stuff, CYCLIC_SAVE, A2V, MAX_PROX_LEFT, MAX_PROX_RIGHT, MAX_BPATHS, pair_sample);
        int MIN_GAIN = 5;
        int EXT_MODE = 1;

        AddNewStuff(new_stuff, hbvr, inv, pathsr, bases, quals, MIN_GAIN, EXT_MODE);
        PartnersToEnds(hbvr, pathsr, bases, quals);
        if(dump_detailed_gfa) GFADumpDetail(out_dir + "/" + out_prefix + ".large_K.final",hbvr,inv);
        if (dump_all || to_step ==5){
            std::cout << Date() << ": Dumping large_K final graph and paths..." << std::endl;
            BinaryWriter::writeFile(out_dir + "/" + out_prefix + ".large_K.final.hbv", hbvr);
            WriteReadPathVec(pathsr,(out_dir + "/" + out_prefix + ".large_K.final.paths").c_str());
            graph_status(hbvr);
            path_status(pathsr);
            std::cout << Date() << ": Dumping large_K final graph and paths DONE!" << std::endl;
        }
        std::cout << Date() << ": Assembling gaps DONE!" << std::endl << std::endl ;

    }



    if (from_step==6){
        std::cout << Date() << ": Reading large_K final graph and paths..." << std::endl;
        BinaryReader::readFile(out_dir + "/" + out_prefix + ".large_K.final.hbv", &hbvr);
        LoadReadPathVec(pathsr,(out_dir + "/" + out_prefix + ".large_K.final.paths").c_str());
        inv.clear();
        hbvr.Involution(inv);
        graph_status(hbvr);
        path_status(pathsr);
        std::cout << Date() << ": Reading large_K final graph and paths DONE!" << std::endl << std::endl;
    }
    if (from_step<=6 and to_step>=6) {
        std::cout << "--== Step 6: Graph simplification and path finding ==--" << std::endl;

        //==Simplify
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

        Simplify(out_dir, hbvr, inv, pathsr, bases, quals, MAX_SUPP_DEL, TAMP_EARLY_MIN, MIN_RATIO2, MAX_DEL2,
                 ANALYZE_BRANCHES_VERBOSE2, TRACE_SEQ, DEGLOOP, EXT_FINAL, EXT_FINAL_MODE,
                 PULL_APART_VERBOSE, PULL_APART_TRACE, DEGLOOP_MODE, DEGLOOP_MIN_DIST, IMPROVE_PATHS,
                 IMPROVE_PATHS_LARGE, FINAL_TINY, UNWIND3, run_pathfinder, dump_pf, pf_verbose);

        // For now, fix paths and write the and their inverse
        for (int i = 0; i < (int) pathsr.size(); i++) { //XXX TODO: change this int for uint 32
            Bool bad = False;
            for (int j = 0; j < (int) pathsr[i].size(); j++)
                if (pathsr[i][j] < 0) bad = True;
            if (bad) pathsr[i].resize(0);
        }
        // TODO: this is "bj making sure the inversion still works", but shouldn't be required
        paths_inv.clear();
        invert(pathsr, paths_inv, hbvr.EdgeObjectCount());

        // Find lines and write files.
        vec<vec<vec<vec<int>>>> lines;

        FindLines(hbvr, inv, lines, MAX_CELL_PATHS, MAX_DEPTH);
        BinaryWriter::writeFile(out_dir + "/" + out_prefix + ".fin.lines", lines);

        // XXX TODO: Solve the {} thingy, check if has any influence in the new code to run that integrated
        {
            vec<int> llens, npairs;
            GetLineLengths(hbvr, lines, llens);
            GetLineNpairs(hbvr, inv, pathsr, lines, npairs);
            BinaryWriter::writeFile(out_dir + "/" + out_prefix + ".fin.lines.npairs", npairs);

            vec<vec<covcount>> covs;
            ComputeCoverage(hbvr, inv, pathsr, lines, subsam_starts, covs);
            //BinaryWriter::writeFile( work_dir + "/" +prefix+ ".fin.covs", covs );
            //WriteLineStats( work_dir, lines, llens, npairs, covs );

            // Report CN stats
            double cn_frac_good = CNIntegerFraction(hbvr, covs);
            std::cout << "CN fraction good = " << cn_frac_good << std::endl;
            PerfStatLogger::log("cn_frac_good", ToString(cn_frac_good, 2), "fraction of edges with CN near integer");
        }

        // TestLineSymmetry( lines, inv2 );
        // Compute fragment distribution.
        FragDist(hbvr, inv, pathsr, out_dir + "/" + out_prefix + ".fin.frags.dist");
        if(dump_detailed_gfa) GFADumpDetail(out_dir + "/" + out_prefix + ".contig",hbvr,inv);
        //TODO: add contig fasta dump.
        if (dump_all || to_step == 6){
            std::cout << Date() << ": Dumping contig graph and paths..." << std::endl;
            BinaryWriter::writeFile(out_dir + "/" + out_prefix + ".contig.hbv", hbvr);
            WriteReadPathVec(pathsr,(out_dir + "/" + out_prefix + ".contig.paths").c_str());
            graph_status(hbvr);
            path_status(pathsr);
            std::cout << Date() << ": Dumping contig graph and paths DONE!" << std::endl;
        }
        //vecbasevector G;
        //FinalFiles(hbvr, inv, pathsr, subsam_names, subsam_starts, out_dir, out_prefix + "_contigs", MAX_CELL_PATHS, MAX_DEPTH, G);
        GFADump(out_dir +"/"+ out_prefix + "_contigs", hbvr, inv, pathsr, MAX_CELL_PATHS, MAX_DEPTH, true);
        std::cout << Date() << ": Contigging DONE!" << std::endl << std::endl;

    }
    if (from_step==7){
        std::cout << Date() << ": Reading contig graph and paths..." << std::endl;
        BinaryReader::readFile(out_dir + "/" + out_prefix + ".contig.hbv", &hbvr);
        LoadReadPathVec(pathsr,(out_dir + "/" + out_prefix + ".contig.paths").c_str());
        inv.clear();
        hbvr.Involution(inv);
        paths_inv.clear();
        invert(pathsr, paths_inv, hbvr.EdgeObjectCount());
        graph_status(hbvr);
        path_status(pathsr);
        std::cout << Date() << ": Reading contig graph and paths DONE!" << std::endl << std::endl;
    }
    if (from_step<=7 and to_step>=7) {
        //== Scaffolding
        std::cout << "--== Step 7: PE-Scaffolding ==--" << std::endl;
        int MIN_LINE = 5000;
        int MIN_LINK_COUNT = 3; //XXX TODO: this variable is the same as -w in soap??

        bool SCAFFOLD_VERBOSE = False;
        bool GAP_CLEANUP = True;


        MakeGaps(hbvr, inv, pathsr, paths_inv, MIN_LINE, MIN_LINK_COUNT, out_dir, out_prefix, SCAFFOLD_VERBOSE, GAP_CLEANUP);

        // Carry out final analyses and write final assembly files.

        vecbasevector G;
        FinalFiles(hbvr, inv, pathsr, subsam_names, subsam_starts, out_dir, out_prefix+ "_assembly", MAX_CELL_PATHS, MAX_DEPTH, G);
        graph_status(hbvr);
        path_status(pathsr);
        GFADump(out_dir +"/"+ out_prefix + "_assembly", hbvr, inv, pathsr, MAX_CELL_PATHS, MAX_DEPTH, true);
        if(dump_detailed_gfa) GFADumpDetail(out_dir + "/" + out_prefix + ".assembly",hbvr,inv);
        std::cout << Date() << ": PE-Scaffolding DONE!" << std::endl << std::endl << std::endl;

    }
    return 0;
}

