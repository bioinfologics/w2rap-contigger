///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <util/OutputLog.h>
#include "Basevector.h"
#include "Bitvector.h"
#include "FetchReads.h"
#include "TokenizeString.h"
#include "math/HoInterval.h"
#include "paths/long/LoadCorrectCore.h"
#include "paths/long/large/ExtractReads.h"

void ExtractReads( String reads, const String& work_dir, vecbvec* pReads, VecPQVec* quals )
{
    reads.GlobalReplaceBy(" ", "");
    vec<String> filenames;
    Tokenize(reads, ',', filenames);
    if (filenames.size()!=2) {
        std::cout<<"ERROR: 2 input files for paired reads needed."<<std::endl;
        Scram(1);
    }
    if (filenames.size()!=2) {
        std::cout<<"ERROR: 2 input files for paired reads needed."<<std::endl;
        Scram(1);
    }
    auto fn1=filenames[0];
    auto fn2=filenames[1];
    std::ifstream in1(filenames[0]);
    std::ifstream in2(filenames[1]);
    std::string line1,line2;

    basevector b1, b2;
    while (1) {
        getline(in1, line1), getline(in2, line2);
        if (in1.fail() && in2.fail()) break;
        if ((in1.fail() && !in2.fail()) || (in2.fail() && !in1.fail())) {
            std::cout << "\nThe files " << fn1 << " and " << fn2
                      << " appear to be paired, yet have "
                      << "different numbers of records.\n" << std::endl;
            Scram(1);
        }

        // Fetch bases.  Turn Ns into As.

        getline(in1, line1), getline(in2, line2);
        if (in1.fail() || in2.fail()) {
            std::cout << "\nSee incomplete record in " << fn1
                      << " or " << fn2 << ".\n" << std::endl;
            Scram(1);
        }
        for (int i = 0; i < line1.size(); i++)
            if (line1[i] == 'N') line1[i] = 'A';
        for (int i = 0; i < line2.size(); i++)
            if (line2[i] == 'N') line2[i] = 'A';
        b1.SetFromString(line1);
        b2.SetFromString(line2);

        // Skip line.

        getline(in1, line1), getline(in2, line2);
        if (in1.fail() || in2.fail()) {
            std::cout << "\nSee incomplete record in " << fn1
                      << " or " << fn2 << ".\n" << std::endl;
            Scram(1);
        }

        // Fetch quals.

        getline(in1, line1), getline(in2, line2);
        if (in1.fail() || in2.fail()) {
            std::cout << "\nSee incomplete record in " << fn1
                      << " or " << fn2 << ".\n" << std::endl;
            Scram(1);
        }
        if (b1.size() != line1.size()
            || b2.size() != line2.size()) {
            std::cout << "\n1: " << b1.size() << " bases "
                      << ", " << line1.size() << " quals" << std::endl;
            std::cout << "2: " << b2.size() << " bases "
                      << ", " << line2.size() << " quals" << std::endl;
            std::cout << "See inconsistent base/quality lengths "
                      << "in " << fn1 << " or " << fn2 << std::endl;
            Scram(1);
        }

        // Save.

        QualVec q1;
        QualVec q2;
        q1.resize(line1.size()), q2.resize(line2.size());
        for (int i = 0; i < line1.size(); i++)
            q1[i] = line1[i] - 33;
        for (int i = 0; i < line2.size(); i++)
            q2[i] = line2[i] - 33;

        pReads->push_back(b1), pReads->push_back(b2);
        quals->emplace_back(PQVec(q1));
        quals->emplace_back(PQVec(q2));
    }
}
