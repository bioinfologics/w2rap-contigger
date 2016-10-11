///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "Basevector.h"
#include "Bitvector.h"
#include "FastIfstream.h"
#include "FetchReads.h"
#include "PairsManager.h"
#include "TokenizeString.h"
#include "bam/ReadBAM.h"
#include "math/HoInterval.h"
#include "paths/long/LoadCorrectCore.h"
#include "paths/long/large/ExtractReads.h"

void ExtractReads( String reads, const String& work_dir, vecbvec* pReads, VecPQVec* quals )
{
     double lclock = WallClockTime( );


     // Parse, check and load files.  This is the main extraction path.

     std::cout << Date() << ": finding input files" << std::endl;
     vecbasevector &xbases = (*pReads);
     VecPQVec &xquals = (*quals);
     std::vector<std::string> infiles_pairs;

     auto infiles = tokenize(reads.c_str(), ',');

     for (auto a: infiles){
          std::cout << a << std::endl;
     }

     // Check that files are OK.

     for (auto fn: infiles){
          if (!IsRegularFile(fn)) {
               std::cout << "\nFailed to find file " << fn << ".\n" << std::endl;
               Scram(1);
          }
          std::vector<std::string> suf = {".bam", ".fastq", ".fastq.gz", ".fq", ".fq.gz", ".fastb"};
          Bool ok = False;
          for (auto s : suf)
               if (fn.find(s) != std::string::npos) ok = True;
          if (!ok) {
               std::cout << "\nInput file " << fn << " has unsupported type.\n"
               << std::endl;
               Scram(1);
          }
     }

     // Read the files.

     std::cout << Date() << ": reading " << infiles.size()
     << " files (which may take a while)" << std::endl;

     for (auto fn: infiles){

          // Parse bam files.

          if (std::string::npos != fn.find(".bam")) {
               bool const UNIQUIFY_NAMES = true;
               vecString *pxnames = 0;
               BAMReader bamReader(False /*USE_PF_ONLY*/,
                                   UNIQUIFY_NAMES,
                                   1, //infiles_meta[g].frac, //TODO: check this fraction, i think is the fraction of reads to use
                                   long(-1/*READS_TO_USE*/)); //TODO: Check this also, i think this is not necesary, never used
               bamReader.readBAM(fn, &xbases, &xquals, pxnames);
          }

          // Parse fastb/qualb/qualp files.

          else if (fn.find(".fastb") != std::string::npos) {
               std::string fn2b = fn.substr(0, fn.find(".fastb")) + ".qualb";
               std::string fn2p = fn.substr(0, fn.find(".fastb")) + ".qualp";
               if (IsRegularFile(fn2b)) {
                    xbases.ReadAll(fn, True);
                    vecqualvector q;
                    q.ReadAll(fn2b);
                    convertAppendParallel(q.begin(), q.end(), xquals);
               }
               else if (IsRegularFile(fn2p)) {
                    xbases.ReadAll(fn, True);
                    xquals.ReadAll(fn2p, True);
               }
          }

          // Parse paired fastq files.

          else if (fn.find(".fastq") != std::string::npos) {
               infiles_pairs.push_back(fn);
               if (infiles_pairs.size() == 2) {
                    const std::string fn1 = infiles_pairs[0];
                    const std::string fn2 = infiles_pairs[1];
                    std::string command1 = "cat " + fn1;
                    std::string command2 = "cat " + fn2;
                    if (fn1.find(".gz") != std::string::npos) command1 = "z" + command1;
                    if (fn2.find(".gz") != std::string::npos) command2 = "z" + command2;
                    fast_pipe_ifstream in1(command1);
                    fast_pipe_ifstream in2(command2);
                    String line1;
                    String line2;

                    // Buffer for quality score compression in batches.

                    const int qbmax = 10000000;
                    vec<qvec> qualsbuf;
                    MempoolOwner<char> alloc;
                    for (int i = 0; i < qbmax; i++)
                         qualsbuf.emplace_back(alloc);
                    int qbcount = 0;

                    // Go through the input files.

                    basevector b1;
                    basevector b2;
                    while (1) {
                         getline(in1, line1), getline(in2, line2);
                         if (in1.fail() && in2.fail()) break;
                         if ((in1.fail() && !in2.fail())
                             || (in2.fail() && !in1.fail())) {
                              std::cout << "\nThe files " << fn1 << " and " << fn2
                              << " appear to be paired, yet have "
                              << "different numbers of records.\n" << std::endl;
                              Scram(1);
                         }

                         getline(in1, line1), getline(in2, line2);
                         if (in1.fail() || in2.fail()) {
                              std::cout << "\nSee incomplete record in " << fn1
                              << " or " << fn2 << ".\n" << std::endl;
                              Scram(1);
                         }

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

                         qvec &q1 = qualsbuf[qbcount++];
                         qvec &q2 = qualsbuf[qbcount++];

                         //TODO: left as is for the moment but this has to be replaced for a better thing !
                         q1.resize(line1.size()), q2.resize(line2.size());
                         if (qbcount == qbmax) {

//######################################################################################################################
                              auto beg = qualsbuf.begin();
                              auto end = qualsbuf.end();

                              size_t nnn = end-beg;
                              xquals.resize(xquals.size()+nnn);
                              auto oItr = xquals.end()-nnn;

                              PQVecEncoder enc;
                              unsigned char* buf = nullptr;
                              PQVec::allocator_type alloc = oItr->get_allocator();

                              for (auto ib = beg; ib != end; ++ib) {
                                   enc.init(*ib);
                                   size_t need = enc.size();
                                   buf = alloc.allocate(need);
                                   oItr->clear().setData(buf);
                                   ++oItr;
                              }

//######################################################################################################################

                         qbcount = 0;
                         }
                         for (int i = 0; i < line1.size(); i++)
                              q1[i] = line1[i] - 33;
                         for (int i = 0; i < line2.size(); i++)
                              q2[i] = line2[i] - 33;

                         xbases.push_back(b1);
                         xbases.push_back(b2);
                    }
                    convertAppendParallel(qualsbuf.begin(),
                                          qualsbuf.begin() + qbcount, xquals);
               }
          }
     }

     // Report stats.

     std::cout << Date() << ": data extraction complete"
     #ifdef __linux
     << ", peak = " << PeakMemUsageGBString( )
     #endif
     << std::endl;
     std::cout << TimeSince(lclock) << " used extracting reads" << std::endl;
}