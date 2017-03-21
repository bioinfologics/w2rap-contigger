//
// Created by Bernardo Clavijo (TGAC) on 21/06/2016.
//
#include <paths/long/large/Repath.h>
#include <paths/long/large/FinalFiles.h>
#include "ParallelVecUtilities.h"
#include "tclap/CmdLine.h"
#include "GFADump.h"

int main(const int argc, const char * argv[]) {

    std::string out_prefix;
    std::string in_prefix;
    bool find_lines, stats_only;

    uint64_t genome_size;
    //========== Command Line Option Parsing ==========

    std::cout << "hbv2gfa from w2rap-contigger" << std::endl;
    try {
        TCLAP::CmdLine cmd("", ' ', "0.1");
        TCLAP::ValueArg<std::string> out_prefixArg("o", "out_prefix",
             "Prefix for the output files", true, "", "string", cmd);
        TCLAP::ValueArg<std::string> in_prefixArg("i", "in_prefix",
                                                   "Prefix for the input files", true, "", "string", cmd);
        TCLAP::ValueArg<uint32_t> genomeSize_Arg("g", "genome_size",
                                                 "Genome size for NGXX stats in Kbp (default: 0, no NGXX stats)", false, 0,"int", cmd);
        TCLAP::ValueArg<bool>         find_linesArg        ("l","find_lines",
                                                            "Find lines", false,false,"bool",cmd);
        TCLAP::ValueArg<bool>         statsOnly_Arg        ("","stats_only",
                                                            "Compute stats only (do not dump GFA)", false,false,"bool",cmd);
        cmd.parse(argc, argv);

        // Get the value parsed by each arg.
        out_prefix = out_prefixArg.getValue();
        in_prefix = in_prefixArg.getValue();
        find_lines = find_linesArg.getValue();
        genome_size = 1000UL * genomeSize_Arg.getValue();
        stats_only = statsOnly_Arg.getValue();

    } catch (TCLAP::ArgException &e)  // catch any exceptions
    {
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
        return 1;
    }
    HyperBasevector hbv;
    ReadPathVec paths;
    vec<int> inv;

    std::cout << "Reading graph and paths..." << std::endl;
    BinaryReader::readFile(in_prefix + ".hbv", &hbv);
    hbv.Involution(inv);
    TestInvolution(hbv,inv);
    LoadReadPathVec(paths,(in_prefix + ".paths").c_str());
    std::cout << "   DONE!" << std::endl;

    std::cout<<"=== Graph stats === "<<std::endl;
    std::vector<uint64_t> e_sizes;
    uint64_t total_size,canonical_size;

    for (int e=0;e<hbv.EdgeObjectCount();e++){
        auto eo=hbv.EdgeObject(e);
        total_size+=eo.size();
        if (eo.getCanonicalForm()==CanonicalForm::FWD or eo.getCanonicalForm()==CanonicalForm::PALINDROME){
            canonical_size+=eo.size();
            e_sizes.push_back(eo.size());
        }
    }
    std::sort(e_sizes.begin(),e_sizes.end());
    auto ns=e_sizes.rbegin();
    int64_t cs=0;
    std::cout<<"Canonical graph sequences size: "<<canonical_size<<std::endl;
    for (auto i=10;i<100;i+=10){
        while (((double)(cs * 100.0))/ canonical_size < i)
            cs+=*ns++;
        std::cout<<"N"<<i<<": "<<*(ns-1)<<std::endl;
    }
    if (genome_size) {
        ns=e_sizes.rbegin();
        cs=0;
        std::cout<<std::endl<<"User provided size: "<<genome_size<<std::endl;
        for (auto i = 10; i < 100; i += 10) {
            while (((double) (cs * 100.0)) / genome_size < i and ns!=e_sizes.rend())
                cs += *ns++;
            if (ns==e_sizes.rend())
                std::cout << "NG" << i << ": n/a" << std::endl;
            else
                std::cout<<"NG" << i << ": " << *(ns-1) << std::endl;
        }
    }

    int MAX_CELL_PATHS = 50;
    int MAX_DEPTH = 10;
    if (!stats_only) {
        std::cout << "Dumping gfa" << std::endl;
        GFADump(out_prefix, hbv, inv, paths, MAX_CELL_PATHS, MAX_DEPTH, find_lines);
    }

    return 0;
}

