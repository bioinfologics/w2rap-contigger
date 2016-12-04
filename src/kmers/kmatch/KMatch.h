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
#include <unordered_map>
//#include "matchresult.h"

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

//Warning: operators are overloaded for search/sort, only comparing the kmer value
typedef struct ekpst{
    uint64_t kmer;
    int edge_id;
    int edge_offset;
    bool operator<(const ekpst &rhs)const{return kmer<rhs.kmer;};
    bool operator==(const ekpst &rhs)const{return kmer==rhs.kmer;};
    bool operator<(const edgeKmerPosition &rhs)const{return kmer<rhs.kmer;};
    bool operator==(const edgeKmerPosition &rhs)const{return kmer==rhs.kmer;};

} edgeKmerPositionNR;


typedef struct{
    int edge_id;
    int edge_offset;
} edgeKmerPositionV;

typedef std::map<uint64_t, std::vector<edgeKmerPositionV>> KMAP_t;

class KMatch {

public:
    KMatch(int K);
    void Hbv2Map(HyperBasevector &hbv);

    void Hbv2Index(const HyperBasevector &hbv, uint step=1);

    std::vector<pKmer> ProduceKmers(const std::string & seq);

    KMAP_t edgeMap;


    std::vector<edgeKmerPositionNR> KmerIndex;

    std::vector<edgeKmerPosition> lookupRead(const std::string & read);
    std::vector<edgeKmerPosition> lookupReadInMap( const std::string & read);
    uint64_t countReadMatches(const std::string & read);

private:
    uint8_t K;

};

#endif //KMATCH_INCLUDED