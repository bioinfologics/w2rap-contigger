#ifndef KMATCH_INCLUDED
#define KMATCH_INCLUDED 1
#include <vector>
#include <list>
#include <algorithm>
#include <iostream>
#include <fstream>
#define __STDC_LIMIT_MACROS
#include <stdint.h>
#include <cstdlib>
#include <paths/HyperBasevector.h>
#include "matchresult.h"

#define KMATCH_MAX_FREQ 4
#define KMATCH_NUC_A 0
#define KMATCH_NUC_C 1
#define KMATCH_NUC_G 2
#define KMATCH_NUC_T 3
#define KMATCH_NOKMER INT64_MIN
#define KMATCH_POSITION_CHR_CNST 10000000000L

// TODO: include this object to replace the vector of pairs thing
typedef struct {
    uint64_t kmer;
    int offset;
} pKmer;

typedef struct {
    uint64_t kmer;
    int edge_id;
    int edge_offset;
    int read_offset;
} edgeKmerPosition;

class KMatch {
  public:
    KMatch(int K);
    void KMatch::Hbv2Map(HyperBasevector* hbv);

    std::vector<pKmer> KMatch::ProduceKmers(std::string seq);

    std::map<uint64_t, std::vector<edgeKmerPosition>> edgeMap;

//    std::vector<int> KMatch::MapReads(vecbvec seqVector, HyperBasevector *hbv);
    std::vector<edgeKmerPosition> KMatch::lookupRead(std::string read);

  private:
    uint8_t K;
};

#endif //KMATCH_INCLUDED
