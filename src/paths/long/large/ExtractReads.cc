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

InputFileReader::InputFileReader() {
}

bool InputFileReader::get_fastq_record(std::ifstream& in, std::tuple<std::string, std::string, std::string> *record){

     std::string nullstr;
     std::string header;
     std::string seq;
     std::string qual;

     getline(in, header);
     if (in.fail()) {
          return false;
     }

     getline(in, seq);

     // Skip line.
     getline(in, nullstr);

     // Fetch quals.
     getline(in, qual);

     std::get<0>(*record) = header;
     std::get<1>(*record) = seq;
     std::get<2>(*record) = qual;
     return true;
}

bool InputFileReader::get_bam_record(){
//               bool const UNIQUIFY_NAMES = true;
//               vecString *pxnames = 0;
//               BAMReader bamReader(False /*USE_PF_ONLY*/,
//                                   UNIQUIFY_NAMES,
//                                   1, //infiles_meta[g].frac, //TODO: check this fraction, i think is the fraction of reads to use
//                                   long(-1/*READS_TO_USE*/)); //TODO: Check this also, i think this is not necesary, never used
//               bamReader.readBAM(fn, &xbases, &xquals, pxnames);
     return true;
}

int InputFileReader::read_binary(std::string out_dir, std::string prefix){
     // thead the fastb &qualp pair
     this->bases.ReadAll(out_dir + "/" + prefix + "frag_reads_orig.fastb");
     this->quals.ReadAll(out_dir + "/" + prefix + "frag_reads_orig.qualp");
}

int InputFileReader::write_binary(std::string out_dir, std::string prefix){
     // Write files fast & qualb
     this->bases.WriteAll(out_dir + "/" + prefix + "frag_reads_orig.fastb");
     this->quals.WriteAll(out_dir + "/" + prefix + "frag_reads_orig.qualp");
}

PeData::PeData(std::string preads_filename){
     //
     //TODO: check zipped files (byte1 == 0x1f) && (byte2 == 0x8b)

     this->filename_string = preads_filename;
     std::vector<std::string> infiles;
     infiles = tokenize(filename_string.c_str(), ',');

     std::cout << Date() << ": finding input files" << std::endl;
     // Check that files are OK.
     for (auto fn: infiles) {
          if (!IsRegularFile(fn)) {
               std::cout << "\nFailed to find file " << fn << ".\n" << std::endl;
               Scram(1);
          }
          std::vector<std::string> suf = {".bam", ".fastq", ".fastq.gz", ".fq", ".fq.gz"};
          Bool ok = False;
          for (auto s : suf)
               if (fn.find(s) != std::string::npos) ok = True;
          if (!ok) {
               std::cout << "\nInput file " << fn << " has unsupported type.\n" << std::endl;
               Scram(1);
          }
     }
     std::cout << Date() << ": reading " << this->infiles_pairs.size() << " files (which may take a while)" << std::endl;

     for (auto a: infiles) {
          // TODO: [GONZA] add check for proper pairing here, headers, numbers and pairs!!!!
          if (a.find(".fastq") != std::string::npos) {
               std::cout << a << std::endl;
               this->infiles_pairs.push_back(a);
          }
     }

     auto ec = this->read_files(&bases, &quals);
}

int PeData::read_files(vecbvec *Reads, VecPQVec *Quals){
     // Read the files.
     double lclock = WallClockTime();

     vecbasevector& pReads = (*Reads);
     VecPQVec& pQuals = (*Quals);

     // Parse bam files.

     if (this->infiles_pairs.size() == 2) {
          const std::string fn1 = this->infiles_pairs[0];
          const std::string fn2 = this->infiles_pairs[1];

          std::ifstream in1(fn1);
          std::ifstream in2(fn2);

          // Buffer for quality score compression in batches.

          const int qbmax = 10000000;
          std::vector<qvec> qualsbuf;
          MempoolOwner<char> alloc;
          for (int i = 0; i < qbmax; i++)
               qualsbuf.emplace_back(alloc);
          int qbcount = 0;

          // Go through the input files.
          std::tuple<std::string, std::string, std::string> record1;
          std::tuple<std::string, std::string, std::string> record2;

          basevector b1;
          basevector b2;
          while (1) {
               // TODO: [GONZA] add the checks

               bool r1_status = this->get_fastq_record(in1, &record1);
               bool r2_status = this->get_fastq_record(in2, &record2);

               if (!r1_status || !r2_status){
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