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
#include "feudal/PQVec.h"
#include "math/HoInterval.h"
#include "paths/long/LoadCorrectCore.h"
#include "paths/long/large/ExtractReads.h"

void GetCannedReferenceSequences( const String& sample, const String& species,
     const String& work_dir );

void ExtractReads( String reads, const String& work_dir, vec<String>& subsam_names,
     vec<int64_t>& subsam_starts, vecbvec* pReads, VecPQVec* quals )
{
     double lclock = WallClockTime( );


     // Parse, check and load files.  This is the main extraction path.

     std::cout << Date() << ": finding input files" << std::endl;
     vecbasevector &xbases = (*pReads);
     VecPQVec &xquals = (*quals);
     std::vector<std::string> infiles_pairs;

     std::string line;

     auto infiles = tokenize(reads.c_str(), ',');

     for (auto a: infiles){
          std::cout << a << std::endl;
     }

     // TODO: Make this 2 steps simpler (type and duplication check)
     // Check that files are OK.

     for (auto fn: infiles){
          if (!IsRegularFile(fn)) {
               std::cout << "\nFailed to find file " << fn << ".\n" << std::endl;
               Scram(1);
          }
          vec<String> suf = {".bam", ".fastq", ".fastq.gz", ".fq", ".fq.gz", ".fastb"};
          Bool ok = False;
          for (auto s : suf)
               if (fn.find(s) != std::string::npos) ok = True;
          if (!ok) {
               std::cout << "\nInput file " << fn << " has unsupported type.\n"
               << std::endl;
               Scram(1);
          }
     }
     // Chech that no file is included twice (in the original code)

     // Read the files.

     std::cout << Date() << ": reading " << infiles.size()
     << " files (which may take a while)" << std::endl;

     for (auto fn: infiles){

          // Parse bam files.

          if (fn.find(".bam") != std::string::npos) {
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
//                    infiles[g][j] = fn.Before(".fastb") + ".{fastb,qualb}"; //TODO: solucionar esto
               }
               else if (IsRegularFile(fn2p)) {
                    xbases.ReadAll(fn, True);
                    xquals.ReadAll(fn2p, True);
//                    infiles[g][j] = fn.Before(".fastb") + ".{fastb,qualp}"; //TODO: solucionar esto
               }
          }

          // Parse paired fastq files.

//          else if (infiles != "" && j < infiles_rn[g].isize() - 1 && infiles_rn[g][j] == infiles_rn[g][j+1]) {
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
                    std::int64_t total = 0;
                    std::int64_t taken = 0;
//                    double frac = infiles_meta[g].frac;

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

                         // Fetch bases.  Turn Ns into As. <------ TODO: TURN Ns into As check impact

                         getline(in1, line1), getline(in2, line2);
                         if (in1.fail() || in2.fail()) {
                              std::cout << "\nSee incomplete record in " << fn1
                              << " or " << fn2 << ".\n" << std::endl;
                              Scram(1);
                         }
                         for (int i = 0; i < line1.size(); ++i)
                              if (line1[i] == 'N') line1[i] = 'A';
                         for (int i = 0; i < line2.size(); ++i)
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


                         qvec &q1 = qualsbuf[qbcount];
                         qbcount++;
                         qvec &q2 = qualsbuf[qbcount];
                         qbcount++;

                         //TODO: left as is for the moment but this has to be replaced for a better thing !
                         q1.resize(line1.size()), q2.resize(line2.size());
                         if (qbcount == qbmax) {
                              convertAppendParallel(qualsbuf.begin(),
                                                    qualsbuf.begin() + qbcount, xquals);
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
