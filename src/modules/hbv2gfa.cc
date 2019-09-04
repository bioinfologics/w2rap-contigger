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
    bool find_lines=false, stats_only=false, translated_paths=false;
    std::string dump_detailed_gfa;

    uint64_t genome_size=0;
    //========== Command Line Option Parsing ==========
    std::vector<std::string> validGFAOpts{"basic", "abyss"};


    cxxopts::Options options(argv[0], "");
    try {
        options.add_options()
                ("i,in_prefix", "Prefix for input files", cxxopts::value(in_prefix))
                ("o,out_prefix", "Prefix for output files", cxxopts::value(out_prefix))
                ("g,genome_size", "Genome size for NGXX stats in Kbp (default: 0, no NGXX stats)", cxxopts::value(genome_size)->default_value("0"))
                ("l,find_lines", "Find lines", cxxopts::value(find_lines))
                ("t,translate_paths", "Write translated paths", cxxopts::value(translated_paths))
                ("stats_only", "Compute stats only (do not dump GFA)", cxxopts::value(stats_only))
                ("dump_detailed_gfa", "Dump detailed GFA for every graph (default: basic)", cxxopts::value(dump_detailed_gfa)->default_value("basic"))
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
    uint64_t total_size=0,canonical_size=0;

    const auto num_edges(hbv.EdgeObjectCount());
    for (int e=0; e<num_edges ;e++){
        const auto& eo=hbv.EdgeObject(e);
        total_size+=eo.size();
        if (eo.getCanonicalForm()==CanonicalForm::FWD or eo.getCanonicalForm()==CanonicalForm::PALINDROME){
            canonical_size+=eo.size();
            e_sizes.push_back(eo.size());
        }
    }
    std::sort(e_sizes.begin(),e_sizes.end());
    auto ns=e_sizes.rbegin();
    int64_t cs=0;
    std::cout << "Number of edges in the graph: " << num_edges << std::endl;
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

    if (/*translated_paths*/ true) {
        std::vector<int> edge_id(num_edges);

        for (int64_t edge=0; edge < num_edges; edge++){
            auto eo = hbv.EdgeObject(edge);
            if (eo.getCanonicalForm()==CanonicalForm::REV){
                edge_id[edge] = -inv[edge];
            } else {
                edge_id[edge] = edge;
            }
        }
        std::string filename(out_prefix+".tpaths");

        std::ofstream f(filename, std::ios::out | std::ios::trunc | std::ios::binary);
        uint64_t pathcount=paths.size();
        f.write((const char *) &pathcount, sizeof(pathcount));
        uint16_t ps;
        int mOffset;
        for (int64_t rpidx=0; rpidx < paths.size(); rpidx++){
            const auto &rp = paths[rpidx];
            mOffset=rp.getOffset();
            ps=rp.size();
            f.write((const char *) &mOffset, sizeof(mOffset));
            f.write((const char *) &ps, sizeof(ps));
            for (int i = 0; i < ps; i++) {
                f.write((const char *) &edge_id[rp[i]], sizeof(int));
            }
        }
        f.close();

    }

    return 0;
}

