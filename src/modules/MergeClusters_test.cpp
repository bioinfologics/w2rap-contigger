#include <iostream>
#include <Vec.h>
#include <util/OutputLog.h>
#include <CoreTools.h>
#include <feudal/BinaryStream.h>
#include <util/OutputLog.h>
#include "CoreTools.h"
#include "Equiv.h"
//#include "ParallelVecUtilities.h"
#include "ParseSet.h"
#include "VecUtilities.h"
#include "paths/HyperBasevector.h"
#include "paths/long/MakeKmerStuff.h"
#include "paths/long/ReadPath.h"
#include "paths/long/large/Unsat.h"
#include "system/SortInPlace.h"

int main(int argc, char **argv) {
    OutputLogLevel = 2;

    HyperBasevector hb;
    vec<int> inv2;
    ReadPathVec paths2;

    //Load hbv
    OutputLog(2) <<"Loading graph..." << std::endl;
    BinaryReader::readFile("step6.hbv", &hb);
    //Create inversion
    OutputLog(4) <<"Creating graph involution..." << std::endl;
    inv2.clear();
    hb.Involution(inv2);
    //load paths
    OutputLog(2) <<"Loading paths..." << std::endl;
    LoadReadPathVec(paths2,"step6.paths");
    //create path inversion
    OutputLog(2) << "Graph and paths loaded" << std::endl << std::endl;

    vec<vec<std::pair<int, int> > > xs,xs2;
    int A2V = 5;

    Unsat2(hb, inv2, paths2, xs2, "./", A2V);
    Unsat(hb, inv2, paths2, xs, "./", A2V);

    return 0;
}