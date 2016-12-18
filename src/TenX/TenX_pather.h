//
// Created by Gonzalo Garcia (TGAC) on 16/12/2016.
//

#ifndef W2RAP_CONTIGGER_TENX_PATHER_H
#define W2RAP_CONTIGGER_TENX_PATHER_H

#include "kmers/kmatch/KMatch.h"
#include "paths/long/ReadPath.h"
#include "paths/long/large/ExtractReads.h"

typedef struct {
    int read_id;
    int read_size;
    int edge_id;
    int inv_edge_id;
    int edge_offset;
    int read_offset;
    int kmer;
} linkReg;

struct linkreg_less_than {
    inline bool operator() (const linkReg& struct1, const linkReg& struct2)
    {
      return (struct1.read_offset < struct2.read_offset);
    }
};

class TenXPather: public KMatch {
public:
    TenXPather(std::vector<tenXRead>* aseqVector, HyperBasevector* ahbv);
//    ReadPathVec mapReads();

    // Qc the run
    int readsQc();

    // Map indexing records per tag and index
    int makeIndexMap(bool to_disc=false);
    int makeTagMap(bool to_disc=false);
    std::map<std::string, std::vector<int>> indexMap;
    std::map<std::string, std::vector<int>> tagMap;

private:
    std::vector<tenXRead>* seqVector;
    HyperBasevector* hbv;





//    std::vector<std::vector<linkReg>> getReadsLinks(bool output_to_file=true);
//    std::vector<linkReg> readOffsetFilter(std::vector<linkReg> data);
//    std::vector<linkReg> matchLengthFilter(std::vector<linkReg> data);
};



#endif //W2RAP_CONTIGGER_TENX_PATHER_H
