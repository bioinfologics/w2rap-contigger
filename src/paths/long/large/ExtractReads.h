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

class InputFileReader{
  public:
    InputFileReader();

    vecbvec bases;
    VecPQVec quals;

    bool InputFileReader::get_fastq_record(std::basic_istream<char>& in, std::tuple<std::string, std::string, std::string> *record);
    bool InputFileReader::get_bam_record();
    int InputFileReader::read_binary(std::string out_dir, std::string library_name);
    int InputFileReader::write_binary(std::string out_dir, std::string library_name);
    bool InputFileReader::FilesExist(std::string infiles);
    bool InputFileReader::FilesExist(std::vector<std::string> infiles);
    bool InputFileReader::IsGz(std::string filename);
    bool InputFileReader::ProduceValidPair(std::string rfilename_string);
    bool InputFileReader::ReadBinaryIfExist(std::string out_dir, std::string library_name);

    // To hold files
    std::string filename_string;
    std::vector<std::string> infiles_pair;
};

class PeData: public InputFileReader{
  public:
    PeData(std::string out_dir, std::string library_name, std::string reads_filename);

  private:
    int PeData::read_files(std::basic_istream<char>& in1, std::basic_istream<char>& in2, vecbvec *Reads, VecPQVec *Quals);
};

class MpData: public InputFileReader{
public:
    MpData(std::string out_dir, std::string library_name, std::string reads_filename);

private:
    int MpData::read_files(std::basic_istream<char>& in1, std::basic_istream<char>& in2, vecbvec *Reads, VecPQVec *Quals);
};

class TenXData: public InputFileReader{
  public:
    TenXData(std::string out_dir, std::string library_name, std::string reads_filename);
    vecbvec rIndexs;
    int TenXData::read_binary(std::string out_dir, std::string prefix);
    int TenXData::write_binary(std::string out_dir, std::string prefix);

  private:
    int TenXData::read_files(std::basic_istream<char>& in1, std::basic_istream<char>& in2, vecbvec *Reads, VecPQVec *Quals, vecbvec *rIndexs);

};

class PacBioData: public InputFileReader{
public:
    PacBioData(std::string out_dir, std::string library_name, std::string read_filename);
    int PacBioData::read_binary(std::string out_dir, std::string prefix);
    int PacBioData::write_binary(std::string out_dir, std::string prefix);

private:
    int PacBioData::read_file(std::basic_istream<char>& in1, vecbvec *Reads, VecPQVec *Quals);
};

void ExtractReads( String reads,
                   const String& work_dir, vecbvec* pReads, VecPQVec* quals );


class InputDataMag{
public:
    InputDataMag(std::string config_file_path, std::string out_dir);
    std::map<std::string, std::unique_ptr<InputFileReader>> mag;
};
#endif
