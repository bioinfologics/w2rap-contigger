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
#include "system/gzstream/gzstream.h"

InputFileReader::InputFileReader() {
}

bool InputFileReader::get_fastq_record(std::basic_istream<char>& in, std::tuple<std::string, std::string, std::string> *record){
     // Get a record from the fastq file

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

     if (seq.length() != qual.length()){
          std::cout << "Sequence and quality have different sizes... broken input files" << std::endl;
          Scram(1);
     }
     return true;
}

bool InputFileReader::get_bam_record(){
     // This is to read bam files (TODO)
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
     // Read the fastb &qualp pair
     this->bases.ReadAll(out_dir + "/" + prefix + "frag_reads_orig.fastb");
     this->quals.ReadAll(out_dir + "/" + prefix + "frag_reads_orig.qualp");
}

int InputFileReader::write_binary(std::string out_dir, std::string prefix){
     // Write files fastb & qualb
     this->bases.WriteAll(out_dir + "/" + prefix + "frag_reads_orig.fastb");
     this->quals.WriteAll(out_dir + "/" + prefix + "frag_reads_orig.qualp");
}

bool InputFileReader::FilesExist(std::string infiles){
     // Check if single file exist and is accessible
     if (!IsRegularFile(infiles)) {
          std::cout << "\nFailed to find file " << infiles << ".\n" << std::endl;
          return false;
     }
     std::vector<std::string> suf = {".bam", ".fastq", ".fastq.gz", ".fq", ".fq.gz"};
     Bool ok = False;
     for (auto s : suf)
          if (infiles.find(s) != std::string::npos) ok = True;
     if (!ok) {
          std::cout << "\nInput file " << infiles << " has unsupported type.\n" << std::endl;
          return false;
     }
     return true;
}

bool InputFileReader::FilesExist(std::vector<std::string> infiles){
     // Check if file exist, is read type and are accessible
     for (auto fn: infiles) {
          if (!IsRegularFile(fn)) {
               std::cout << "\nFailed to find file " << fn << ".\n" << std::endl;
               return false;
          }
          std::vector<std::string> suf = {".bam", ".fastq", ".fastq.gz", ".fq", ".fq.gz"};
          Bool ok = False;
          for (auto s : suf)
               if (fn.find(s) != std::string::npos) ok = True;
          if (!ok) {
               std::cout << "\nInput file " << fn << " has unsupported type.\n" << std::endl;
               return false;
          }
     }
     return true;
}

bool InputFileReader::IsGz(std::string filename){
     // Return true if the file is .gz
     if(filename.substr( filename.length() - 3 ) == ".gz"){
          return true;
     }
     return false;
}

bool InputFileReader::ProduceValidPair(std::string filename_string){
     // Get the filename string from the arguments and produce a valid pair of files to process

     std::vector<std::string> infiles;
     infiles = tokenize(filename_string.c_str(), ',');

     if (!this->FilesExist(infiles)) Scram(1);

     std::cout << Date() << ": reading " << this->infiles_pair.size() << " files (which may take a while)" << std::endl;
     for (auto a: infiles) {
          if (a.find(".fastq") != std::string::npos) {
               std::cout << a << std::endl;
               this->infiles_pair.push_back(a);
          }
     }

     if (this->infiles_pair.size() != 2) {
          std::cout << "Error there are " << this->infiles_pair.size() << "files in the list but there should be 2, exiting" << std::endl;
          Scram(1);
     }

     return true;
}

// -------------- Pe Data -------------
PeData::PeData(std::string reads_filename){
     //

     this->filename_string = reads_filename;


     std::cout << Date() << ": finding input files" << std::endl;

     if (!ProduceValidPair(this->filename_string)) Scram(1);

     const std::string fn1 = this->infiles_pair[0];
     const std::string fn2 = this->infiles_pair[1];

     // check if the file is gzip
     if (this->IsGz(fn1) && this->IsGz(fn2)){
          std::cout << "File1 is gzipped changing stream" << std::endl;
          igzstream in1(fn1.c_str());
          igzstream in2(fn2.c_str());
          auto ec = this->read_files(in1, in2, &bases, &quals);
     } else {
          std::ifstream in1(fn1);
          std::ifstream in2(fn2);
          auto ec = this->read_files(in1, in2, &bases, &quals);
     }
}

int PeData::read_files(std::basic_istream<char>& in1, std::basic_istream<char>& in2, vecbvec *Reads, VecPQVec *Quals){
     // Read the files.
     double lclock = WallClockTime();

     vecbasevector& pReads = (*Reads);
     VecPQVec& pQuals = (*Quals);

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

     // Report stats.

     std::cout << Date() << ": data extraction complete"
     #ifdef __linux
     << ", peak = " << PeakMemUsageGBString( )
     #endif
     << std::endl;
     std::cout << TimeSince(lclock) << " used extracting reads" << std::endl;
     return 0;
}

// -------------- Mp Data -------------
MpData::MpData(std::string reads_filename){

     this->filename_string = reads_filename;
     if (!ProduceValidPair(this->filename_string)) Scram(1);

     const std::string fn1 = this->infiles_pair[0];
     const std::string fn2 = this->infiles_pair[1];

     // check if the file is gzip
     if (this->IsGz(fn1) && this->IsGz(fn2)){
          std::cout << "File1 is gzipped changing stream" << std::endl;
          igzstream in1(fn1.c_str());
          igzstream in2(fn2.c_str());
          auto ec = this->read_files(in1, in2, &bases, &quals);
     } else {
          std::ifstream in1(fn1);
          std::ifstream in2(fn2);
          auto ec = this->read_files(in1, in2, &bases, &quals);
     }

}

int MpData::read_files(std::basic_istream<char>& in1, std::basic_istream<char>& in2, vecbvec *Reads, VecPQVec *Quals){
     // Read the files.
     double lclock = WallClockTime();

     vecbasevector& pReads = (*Reads);
     VecPQVec& pQuals = (*Quals);

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

     // Report stats.

     std::cout << Date() << ": data extraction complete"
     #ifdef __linux
     << ", peak = " << PeakMemUsageGBString( )
     #endif
     << std::endl;
     std::cout << TimeSince(lclock) << " used extracting reads" << std::endl;
     return 0;
}

// -------------- 10X Data -------------
TenXData::TenXData(std::string reads_filename){
     //

     this->filename_string = reads_filename;
     if (!ProduceValidPair(this->filename_string)) Scram(1);

     const std::string fn1 = this->infiles_pair[0];
     const std::string fn2 = this->infiles_pair[1];

     // check if the file is gzip
     if (this->IsGz(fn1) && this->IsGz(fn2)){
          std::cout << "File1 is gzipped changing stream" << std::endl;
          igzstream in1(fn1.c_str());
          igzstream in2(fn2.c_str());
          auto ec = this->read_files(in1, in2, &bases, &quals, &rIndexs);
     } else {
          std::ifstream in1(fn1);
          std::ifstream in2(fn2);
          auto ec = this->read_files(in1, in2, &bases, &quals, &rIndexs);
     }
}

int TenXData::read_files(std::basic_istream<char>& in1, std::basic_istream<char>& in2, vecbvec *Reads, VecPQVec *Quals, vecbvec *rIndexs){
     // Read the files.
     double lclock = WallClockTime();

     vecbasevector& pReads = (*Reads);
     VecPQVec& pQuals = (*Quals);
     vecbasevector& prIndexs= (*rIndexs);

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
     basevector idx;
     while (1) {
          // TODO: [GONZA] add the checks

          bool r1_status = this->get_fastq_record(in1, &record1);
          bool r2_status = this->get_fastq_record(in2, &record2);

          if (!r1_status || !r2_status){
               break;
          }

          idx.SetFromString(std::get<1>(record1).substr(0,16));
          b1.SetFromString(std::get<1>(record1).substr(16));
          b2.SetFromString(std::get<1>(record2));

          pReads.push_back(b1);
          prIndexs.push_back(idx);
          pReads.push_back(b2);

          // Save.
          auto qc1 = std::get<2>(record1).substr(16);
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

     // Report stats.

     std::cout << Date() << ": data extraction complete"
     #ifdef __linux
     << ", peak = " << PeakMemUsageGBString( )
     #endif
     << std::endl;
     std::cout << TimeSince(lclock) << " used extracting reads" << std::endl;
     return 0;
}

int TenXData::read_binary(std::string out_dir, std::string prefix){
     // thead the fastb &qualp pair
     this->bases.ReadAll(out_dir + "/" + prefix + "frag_reads_orig.fastb");
     this->bases.ReadAll(out_dir + "/" + prefix + "frag_reads_orig.idxb");
     this->quals.ReadAll(out_dir + "/" + prefix + "frag_reads_orig.qualp");
}

int TenXData::write_binary(std::string out_dir, std::string prefix){
     // Write files fast & qualb
     this->bases.WriteAll(out_dir + "/" + prefix + "frag_reads_orig.fastb");
     this->bases.WriteAll(out_dir + "/" + prefix + "frag_reads_orig.idxb");
     this->quals.WriteAll(out_dir + "/" + prefix + "frag_reads_orig.qualp");
}

// -------------- PacBio -------------
PacBioData::PacBioData(std::string read_filename){
     //

     this->filename_string = read_filename;
     std::cout << Date() << ": finding input files" << std::endl;

     // Check that files are OK.
     if (!this->FilesExist(this->filename_string)) Scram(1);
     std::cout << Date() << ": reading 1 file (which may take a while)" << std::endl;
     const std::string fn = read_filename;

     // check if the file is gzip
     if (this->IsGz(fn)){
          std::cout << "File1 is gzipped changing stream" << std::endl;
          igzstream in(fn.c_str());
          auto ec = this->read_file(in, &bases, &quals);
     } else {
          std::ifstream in(fn);
          auto ec = this->read_file(in, &bases, &quals);
     }
}

int PacBioData::read_file(std::basic_istream<char>& in1, vecbvec *Reads, VecPQVec *Quals){
     // Read the files.
     double lclock = WallClockTime();

     vecbasevector& pReads = (*Reads);
     VecPQVec& pQuals = (*Quals);

     // Buffer for quality score compression in batches.
     const int qbmax = 10000000;
     std::vector<qvec> qualsbuf;
     MempoolOwner<char> alloc;
     for (int i = 0; i < qbmax; i++)
          qualsbuf.emplace_back(alloc);
     int qbcount = 0;

     // Go through the input files.
     std::tuple<std::string, std::string, std::string> record1;
     basevector b1;

     while (1) {

          bool r1_status = this->get_fastq_record(in1, &record1);
          if (!r1_status){
               break;
          }

          b1.SetFromString(std::get<1>(record1));
          pReads.push_back(b1);

          // Save.
          auto qc1 = std::get<2>(record1);
          qvec &q1 = qualsbuf[qbcount++];

          q1.resize(qc1.size());
          if (qbcount == qbmax) {
               convertAppendParallel(qualsbuf.begin(),
                                     qualsbuf.begin() + qbcount, pQuals);
               qbcount = 0;
          }
          for (int i = 0; i < qc1.size(); i++)
               q1[i] = qc1[i] - 33;
     }
     convertAppendParallel(qualsbuf.begin(),
                           qualsbuf.begin() + qbcount, pQuals);

     // Report stats.

     std::cout << Date() << ": data extraction complete"
     #ifdef __linux
     << ", peak = " << PeakMemUsageGBString( )
     #endif
     << std::endl;
     std::cout << TimeSince(lclock) << " used extracting reads" << std::endl;
     return 0;
}

int PacBioData::read_binary(std::string out_dir, std::string prefix){
     this->bases.ReadAll(out_dir + "/" + prefix + "frag_reads_orig.fastb");
}

int PacBioData::write_binary(std::string out_dir, std::string prefix){
     this->bases.WriteAll(out_dir + "/" + prefix + "frag_reads_orig.fastb");
}

// -------------- To read the data from the configuration file
InputDataMag::InputDataMag(std::string config_file_path){
     // Open and read the configuration file
     std::ifstream cf(config_file_path);

     std::string line;
     while(getline(cf, line)){
          std::cout << line << std::endl;
          std::vector<std::string> sline = tokenize(line.c_str(), ' ');
          if (sline[1] == "pe"){
               std::cout << "Lib: " << sline[0] << " is a pe library in " << sline[2] << std::endl;
               PeData ped(sline[2]);
               mag.insert(std::pair<std::string, InputFileReader>(sline[0], ped));
          } else if (sline[1] == "mp") {
               std::cout << "Lib: " << sline[0] << " is a mp library in " << sline[2] << std::endl;
               MpData ped(sline[2]);
               mag.insert(std::pair<std::string, InputFileReader>(sline[0], ped));
          } else if (sline[1] == "10x") {
               std::cout << "Lib: " << sline[0] << " is a 10x library in " << sline[2] << std::endl;
               TenXData ped(sline[2]);
               mag.insert(std::pair<std::string, InputFileReader>(sline[0], ped));
          } else if (sline[1] == "pb") {
               std::cout << "Lib: " << sline[0] << " is a pb library in " << sline[2] << std::endl;
               PacBioData ped(sline[2]);
               mag.insert(std::pair<std::string, InputFileReader>(sline[0], ped));
          }
     }
}