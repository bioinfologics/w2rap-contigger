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
    vec<vec<std::pair<int, int> > > xs;
    vec<vec<int>> n;
    BinaryReader::readFile("xs.out", &xs);
    BinaryReader::readFile("n.out", &n);
    int merge_passes(10);
    auto prev_xs(xs.size());
    OutputLogLevel = 2;

    for (int p = 1; p <= merge_passes; p++) {
        OutputLog(2) << "Merge clusters \n";
        OutputLog(2) << xs.size() << " elements in xs\n";
        prev_xs = xs.size();

        MergeClusters(xs, xs, n, n.size());

        // Check if made any change, if hasn't simply bomb out
        if (prev_xs == xs.size()) {
            OutputLog(2)<<"Completed using " << p << " / " << merge_passes << " passes, last pass made no difference" << std::endl;
            break;
        }
        OutputLog(2)<<"Completed pass " << p << " / " << merge_passes << std::endl;

    }

    return 0;
}