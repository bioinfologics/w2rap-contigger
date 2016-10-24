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

class InputFileReader{
public:
    InputFileReader();

    vecbvec bases;
    VecPQVec quals;

    bool get_fastq_record(std::ifstream& in, std::tuple<std::string, std::string, std::string> *record);
    bool get_bam_record();
    int InputFileReader::read_binary(std::string workdir, std::string prefix);
    int InputFileReader::write_binary(std::string workdir, std::string prefix);

};

class PeData: public InputFileReader{
public:
  PeData(std::string reads_filename);

private:
  std::string filename_string;
  std::vector<std::string> infiles_pairs;

  int read_files(vecbvec *Reads, VecPQVec *Quals);
};

void ExtractReads( String reads,
                   const String& work_dir, vecbvec* pReads, VecPQVec* quals );

#endif
