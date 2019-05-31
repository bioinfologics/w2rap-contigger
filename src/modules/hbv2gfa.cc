//
// Created by Bernardo Clavijo (TGAC) on 21/06/2016.
//
#include <paths/long/large/Repath.h>
#include <paths/long/large/FinalFiles.h>
#include <deps/cxxopts/cxxopts.hpp>
#include "ParallelVecUtilities.h"
#include "GFADump.h"

int main( int argc, char * argv[]) {

    std::string out_prefix;
    std::string in_prefix;
    bool find_lines, stats_only;
    std::string dump_detailed_gfa;

    uint64_t genome_size;
    //========== Command Line Option Parsing ==========
    std::vector<std::string> validGFAOpts{"basic", "abyss"};


    cxxopts::Options options(argv[0], "");
    try {
        options.add_options()
                ("i,in_prefix", "Prefix for input files", cxxopts::value(in_prefix))
                ("o,out_prefix", "Prefix for output files", cxxopts::value(out_prefix))
                ("g,genome_size", "Genome size for NGXX stats in Kbp (default: 0, no NGXX stats)", cxxopts::value(genome_size))
                ("l,find_lines", "Find lines", cxxopts::value(find_lines))
                ("stats_only", "Compute stats only (do not dump GFA)", cxxopts::value(stats_only))
                ("dump_detailed_gfa", "Dump detailed GFA for every graph (default: basic)", cxxopts::value(dump_detailed_gfa))
                ("h,help","show help message")
                ;

        auto result = options.parse(argc, argv);

        genome_size*=1000UL;
        if (result.count("help")) {
            std::cout << options.help({""}) << std::endl;
            exit(0);
        }

        if (!result.count("i")) {
            std::cerr << "Required option i,in_prefix missing" << std::endl;
            std::cout << options.help({""}) << std::endl;
            exit(1);
        }
        if (!result.count("o")) {
            std::cerr << "Required option o,out_prefix missing" << std::endl;
            std::cout << options.help({""}) << std::endl;
            exit(1);
        }

        if (std::find(validGFAOpts.begin(),validGFAOpts.end(), dump_detailed_gfa) == validGFAOpts.cend()) {
            std::cerr << "Invalid argument for dump_detailed_gfa. Please choose between options: ";
            std::copy(validGFAOpts.begin(), validGFAOpts.end(), std::ostream_iterator<std::string>(std::cerr, ", "));
            std::cerr << std::endl;
            std::cout << options.help({""}) << std::endl;
            exit(1);
        }

    } catch (const cxxopts::OptionException& e) {
        std::cerr << "error parsing options: " << e.what() << std::endl;
        std::cout << options.help({""}) << std::endl;
        exit(1);
    }

    std::cout << "hbv2gfa from w2rap-contigger" << std::endl;
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
        if (validGFAOpts[0] == dump_detailed_gfa) {GFADump(out_prefix, hbv, inv, paths, MAX_CELL_PATHS, MAX_DEPTH, find_lines);}
        else if (validGFAOpts[1] == dump_detailed_gfa) {GFADumpAbyss(out_prefix, hbv, inv, paths, MAX_CELL_PATHS, MAX_DEPTH, find_lines);}
    }

    return 0;
}

