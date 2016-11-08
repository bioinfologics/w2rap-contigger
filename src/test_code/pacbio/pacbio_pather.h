//
// Created by Gonzalo Garcia (TGAC) on 04/11/2016.
//

#ifndef W2RAP_CONTIGGER_PACBIO_PATHER_H
#define W2RAP_CONTIGGER_PACBIO_PATHER_H

#include "kmers/kmatch/KMatch.h"

class PacbioPather: public KMatch {
  public:
    PacbioPather::PacbioPather();
    void PacbioPather::mapReads(vecbvec seqVector, HyperBasevector* hbv);
};


#endif //W2RAP_CONTIGGER_PACBIO_PATHER_H
