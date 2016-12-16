///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef FILE_READER_H
#define FILE_READER_H

#include "Basevector.h"
#include "CoreTools.h"
#include "feudal/PQVec.h"
#include "system/gzstream/gzstream.h"

typedef struct {
    basevector r1;
    basevector r2;
    basevector index;
    basevector tag;
} tenXRead;

class InputFileReader{
  public:
    InputFileReader();

    vecbvec bases;
    VecPQVec quals;
    std::vector<tenXRead> rReads;

    bool get_fastq_record(std::basic_istream<char>& in, std::tuple<std::string, std::string, std::string> *record);
    bool get_bam_record();
    int read_binary(std::string out_dir, std::string library_name);
    int write_binary(std::string out_dir, std::string library_name);
    bool FilesExist(std::string infiles);
    bool FilesExist(std::vector<std::string> infiles);
    bool IsGz(std::string filename);
    bool ProduceValidPair(std::string rfilename_string);
    bool ReadBinaryIfExist(std::string out_dir, std::string library_name);

    // To hold files
    std::string filename_string;
    std::vector<std::string> infiles_pair;
};

class PeData: public InputFileReader{
  public:
    PeData(std::string out_dir, std::string library_name, std::string reads_filename);

  private:
    int read_files(std::basic_istream<char>& in1, std::basic_istream<char>& in2, vecbvec *Reads, VecPQVec *Quals);
};

class MpData: public InputFileReader{
public:
    MpData(std::string out_dir, std::string library_name, std::string reads_filename);

private:
    int read_files(std::basic_istream<char>& in1, std::basic_istream<char>& in2, vecbvec *Reads, VecPQVec *Quals);
};



class TenXData: public InputFileReader{
  public:
    TenXData(std::string out_dir, std::string library_name, std::string reads_filename);

    int read_binary(std::string out_dir, std::string prefix);
    int write_binary(std::string out_dir, std::string prefix);

  private:
    int read_files(std::basic_istream<char>& in1, std::basic_istream<char>& in2, std::basic_istream<char>& ini);
    bool split_filenames(std::string filestring);

};

class PacBioData: public InputFileReader{
public:
    PacBioData(std::string out_dir, std::string library_name, std::string read_filename);
    int read_binary(std::string out_dir, std::string prefix);
    int write_binary(std::string out_dir, std::string prefix);

private:
    int read_file(std::basic_istream<char>& in1, vecbvec *Reads, VecPQVec *Quals);
};

void ExtractReads( String reads,
                   const String& work_dir, vecbvec* pReads, VecPQVec* quals );


class InputDataMag{
public:
    InputDataMag(std::string config_file_path, std::string out_dir);
    std::map<std::string, std::unique_ptr<InputFileReader>> mag;
};
#endif
