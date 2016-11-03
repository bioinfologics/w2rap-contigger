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
struct {
    int edge_idx;
    int offset;
} edge_kmer;

class KMatch {
  public:
    KMatch(int K);
    void KMatch::Hbv2Map(HyperBasevector* hbv);
    std::vector<int> KMatch::MapReads(vecbvec seqVector);
    std::vector<std::pair<uint_least64_t, int>> KMatch::ProduceKmers(std::string seq);
    std::map<uint64_t, std::vector<std::pair<int, int>>> edgeMap;

  private:
    uint8_t K;
    std::vector<std::pair<uint64_t, int>> KMatch::lookupRead(std::string read);


};

#endif //KMATCH_INCLUDED
