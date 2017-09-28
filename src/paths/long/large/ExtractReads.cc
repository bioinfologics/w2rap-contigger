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
#include "util/kseq.hpp"

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

    kseq seq1,seq2;
    int l1=0;
    int c1=0;

    int l2=0;
    int c2=0;

    gzFile fp1 = gzopen(fn1.c_str(), "r");
    gzFile fp2 = gzopen(fn2.c_str(), "r");
    FunctorZlib gzr1, gzr2;
    kstream<gzFile, FunctorZlib> ks1(fp1,gzr1);
    kstream<gzFile, FunctorZlib> ks2(fp2,gzr2);
    l1 = ks1.read(seq1);
    l2 = ks2.read(seq2);
    int finished_early=0;
    basevector b1, b2;
    while ( l1 >= 0 and l2 >= 0 ){

        if (seq1.seq.empty() ) {
            std::cout << "Error " << std::string(fn1.c_str()) << " on read " << c1 << " is invalid" << std::endl;
            finished_early=1;
            break;
        }
        if (seq2.seq.empty()) {
            std::cout << "Error " << std::string(fn2.c_str()) << " on read " << c2 << " is invalid" << std::endl;
            finished_early=1;
            break;
        }

        //Create seq
        for (char &i : seq1.seq)
            if (i == 'N') i = 'A';
        for (char &i : seq2.seq)
            if (i == 'N') i = 'A';
        b1.SetFromString(seq1.seq);
        b2.SetFromString(seq2.seq);
        //Create qual
        QualVec q1;
        QualVec q2;
        q1.resize(seq1.qual.size()), q2.resize(seq2.qual.size());
        for (int i = 0; i < seq1.qual.size(); i++)
            q1[i] = seq1.qual[i] - 33;
        for (int i = 0; i < seq2.qual.size(); i++)
            q2[i] = seq2.qual[i] - 33;
        //Store
        pReads->push_back(b1), pReads->push_back(b2);
        quals->emplace_back(PQVec(q1));
        quals->emplace_back(PQVec(q2));

        //std::cout << seq1.name << " " << seq1.comment << "\t" << seq2.name << " " << seq2.comment << std::endl;
        //std::cout << seq1.qual << "\t" << seq2.qual << std::endl;
        c1++;
        c2++;
        l1 = ks1.read(seq1);
        l2 = ks2.read(seq2);
    }

    if (!finished_early) {
        if (l1 != l2) {
            std::cout << "Error on files after " << c1 << " reads, check they are correctly paired" << std::endl;
        }
        if (l1 == -2) {
            std::cout << std::string(fn1.c_str()) << " is invalid at " << c1 << std::endl;
        }
        if (l2 == -2) {
            std::cout << std::string(fn2.c_str()) << " is invalid at " << c2 << std::endl;
        }
    }

    gzclose(fp1);
    gzclose(fp2);
}
