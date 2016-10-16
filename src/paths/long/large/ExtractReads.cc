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
#include "FetchReads.h"
#include "TokenizeString.h"
#include "bam/ReadBAM.h"
#include "math/HoInterval.h"
#include "paths/long/LoadCorrectCore.h"
#include "paths/long/large/ExtractReads.h"

std::tuple<std::string, std::string, std::string> InputFileReader::get_record(std::ifstream& in){

     std::string nullstr;
     std::string header;
     std::string seq;
     std::string qual;

     getline(in, header);
//     if (in.fail()) {
//          std::cout <<" Something wrong with thte record " << std::endl;
//          Scram(1);
//     }

     getline(in, seq);
//     if (in.fail()) {
//          std::cout << "\nSee incomplete record in\n" << std::endl;
//          Scram(1);
//     }

     // Skip line.
     getline(in, nullstr);
//     if (in.fail()) {
//          std::cout << "\nSee incomplete record in.\n" << std::endl;
//          Scram(1);
//     }

     // Fetch quals.
     getline(in, qual);
//     if (in.fail()) {
//          std::cout << "\nSee incomplete record in.\n" << std::endl;
//          Scram(1);
//     }
//     if (seq.size() != qual.size()) {
//          std::cout << "\n1: " << seq.size() << " bases "
//          << ", " << seq.size() << " quals" << std::endl;
//          Scram(1);
//     }
//     std::cout<<header.size()<< " "<<seq.size()<<" "<<qual.size()<<std::endl;
     return std::make_tuple(header, seq, qual);
}

InputFileReader::InputFileReader(const std::string pfilename_string) {

     // Parse, check and load files.  This is the main extraction path.

     this->filename_string = pfilename_string;

     std::vector<std::string> infiles;
     infiles = tokenize(filename_string.c_str(), ',');

     std::cout << Date() << ": finding input files" << std::endl;
     // Check that files are OK.
     for (auto fn: infiles) {
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

     for (auto a: infiles) {
          // TODO: [GONZA] add check for proper pairing here, headers, numbers and pairs!!!!
          if (a.find(".fastq") != std::string::npos) {
               std::cout << a << std::endl;
               this->infiles_pairs.push_back(a);
          }
     }
}

int InputFileReader::read_file(vecbvec *Reads, VecPQVec *Quals){
     // Read the files.
     double lclock = WallClockTime();

     vecbasevector& pReads = (*Reads);
     VecPQVec& pQuals = (*Quals);

     std::cout << Date() << ": reading " << this->infiles_pairs.size()
     << " files (which may take a while)" << std::endl;

          // Parse bam files.

//          if (std::string::npos != fn.find(".bam")) {
//               bool const UNIQUIFY_NAMES = true;
//               vecString *pxnames = 0;
//               BAMReader bamReader(False /*USE_PF_ONLY*/,
//                                   UNIQUIFY_NAMES,
//                                   1, //infiles_meta[g].frac, //TODO: check this fraction, i think is the fraction of reads to use
//                                   long(-1/*READS_TO_USE*/)); //TODO: Check this also, i think this is not necesary, never used
//               bamReader.readBAM(fn, &xbases, &xquals, pxnames);
//          }

          // Parse paired fastq files.

     if (this->infiles_pairs.size() == 2) {
          const std::string fn1 = this->infiles_pairs[0];
          const std::string fn2 = this->infiles_pairs[1];

          std::ifstream in1(fn1);
          std::ifstream in2(fn2);

          // Buffer for quality score compression in batches.

          const int qbmax = 10000000;
          std::vector<qvec> qualsbuf; // TODO: [GONZA] replace this for a standard vector!!!
          MempoolOwner<char> alloc;
          for (int i = 0; i < qbmax; i++)
               qualsbuf.emplace_back(alloc);
          int qbcount = 0;

          // Go through the input files.

          basevector b1;
          basevector b2;
          bool more_records = true;
          while (more_records) {
               // TODO: [GONZA] add the checks
               auto record1 = this->get_record(in1);
               auto record2 = this->get_record(in2);

               if (std::get<1>(record1).size() < 10){
                    break;
               }

               b1.SetFromString(std::get<1>(record1));
               b2.SetFromString(std::get<1>(record2));


               pReads.push_back(b1);
               pReads.push_back(b2);

               // Save.
               auto qc1 = std::get<2>(record1);
               auto qc2 = std::get<2>(record2);

               qvec &q1 = qualsbuf[qbcount++];
               qvec &q2 = qualsbuf[qbcount++];

               q1.resize(qc1.size()), q2.resize(qc2.size());
               if (qbcount == qbmax) {
                    convertAppendParallel(qualsbuf.begin(),
                                          qualsbuf.begin() + qbcount, pQuals);
                    qbcount = 0;
               }
               for (int i = 0; i < qc1.size(); i++)
                    q1[i] = qc1[i] - 33;
               for (int i = 0; i < qc2.size(); i++)
                    q2[i] = qc2[i] - 33;
          }
          convertAppendParallel(qualsbuf.begin(),
                                qualsbuf.begin() + qbcount, pQuals);
     }

     // Report stats.

     std::cout << Date() << ": data extraction complete"
     #ifdef __linux
     << ", peak = " << PeakMemUsageGBString( )
     #endif
     << std::endl;
     std::cout << TimeSince(lclock) << " used extracting reads" << std::endl;
     return 0;
}