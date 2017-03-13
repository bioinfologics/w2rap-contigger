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
#include <unordered_map>

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

class edgeKmerPosition {
public:
    edgeKmerPosition(uint64_t aKmer, int aEdge_id, int aEdge_offset, int aRead_offset):
            kmer(aKmer), edge_id(aEdge_id), edge_offset(aEdge_offset), read_offset(aRead_offset) {};

    edgeKmerPosition(uint64_t aKmer, int aEdge_id, int aEdge_offset):
            kmer(aKmer), edge_id(aEdge_id), edge_offset(aEdge_offset) {};

    uint64_t kmer;
    int edge_id;
    int edge_offset;
    int read_offset; // TODO: Check where i'm using this one value
};

class KMatch {
  public:
    KMatch(const int K);
    void Hbv2Map(const HyperBasevector &hbv);

    std::vector<pKmer> ProduceKmers(const std::string &seq) const;

    std::unordered_map<uint64_t, std::vector<edgeKmerPosition>> edgeMap;

//    std::vector<int> MapReads(vecbvec seqVector, HyperBasevector *hbv);
    std::vector<edgeKmerPosition> lookupRead(const std::string &read) const;

  private:
    uint8_t K;
};

#endif //KMATCH_INCLUDED
