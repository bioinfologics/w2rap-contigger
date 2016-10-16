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
    InputFileReader(const std::string reads_filename);
    std::tuple<std::string, std::string, std::string> get_record(std::ifstream& in);
    int read_file(vecbvec *Reads, VecPQVec *Quals);

private:
    std::string filename_string;
    std::vector<std::string> infiles_pairs;
//    vecbasevector& pReads;
//    VecPQVec& pQuals;
};

void ExtractReads( String reads,
                   const String& work_dir, vecbvec* pReads, VecPQVec* quals );

#endif
