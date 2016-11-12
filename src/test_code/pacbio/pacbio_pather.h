//
// Created by Gonzalo Garcia (TGAC) on 04/11/2016.
//

#ifndef W2RAP_CONTIGGER_PACBIO_PATHER_H
#define W2RAP_CONTIGGER_PACBIO_PATHER_H

#include "kmers/kmatch/KMatch.h"

typedef struct {
    int read_id;
    int read_size;
    int edge_id;
    int inv_edge_id;
    int edge_offset;
    int read_offset;
    int kmer;
} linkReg;

class PacbioPather: public KMatch {
public:
    PacbioPather::PacbioPather(vecbvec* seqVector, HyperBasevector* hbv);
    void PacbioPather::mapReads();

private:
    vecbvec* seqVector;
    HyperBasevector* hbv;

    std::vector<linkReg> PacbioPather::getReadsLinks(bool output_to_file=true);
    std::vector<linkReg> PacbioPather::readOffsetFilter(std::vector<linkReg> data);
    std::vector<linkReg> PacbioPather::readLinksFilter(std::vector<linkReg> data, int read_id);

};


#endif //W2RAP_CONTIGGER_PACBIO_PATHER_H
