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
#include "GFADump.h"

int main(const int argc, const char * argv[]) {

    std::string out_prefix;
    std::string in_prefix;
    bool find_lines;

    //========== Command Line Option Parsing ==========

    std::cout << "hbv2gfa from w2rap-contigger" << std::endl;
    try {
        TCLAP::CmdLine cmd("", ' ', "0.1");
        TCLAP::ValueArg<std::string> out_prefixArg("o", "prefix",
             "Prefix for the output files", true, "", "string", cmd);
        TCLAP::ValueArg<std::string> in_prefixArg("i", "prefix",
                                                   "Prefix for the input files", true, "", "string", cmd);

        TCLAP::ValueArg<bool>         find_linesArg        ("l","find_lines",
                                                         "Find lines", false,false,"bool",cmd);
        cmd.parse(argc, argv);

        // Get the value parsed by each arg.
        out_prefix = out_prefixArg.getValue();
        in_prefix = out_prefixArg.getValue();
        find_lines = find_linesArg.getValue();

    } catch (TCLAP::ArgException &e)  // catch any exceptions
    {
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
        return 1;
    }
    HyperBasevector hbv;
    ReadPathVec paths;
    vec<int> inv;
    hbv.Involution(inv);
    std::cout << "Reading graph and paths..." << std::endl;
    BinaryReader::readFile(in_prefix + ".hbv", &hbv);
    paths.ReadAll(in_prefix + ".paths");
    std::cout << "   DONE!" << std::endl;
    std::cout << "Dumping gfa" << std::endl;
    int MAX_CELL_PATHS = 50;
    int MAX_DEPTH = 10;

    GFADump(out_prefix, hbv, inv, paths, MAX_CELL_PATHS, MAX_DEPTH);

    return 0;
}

