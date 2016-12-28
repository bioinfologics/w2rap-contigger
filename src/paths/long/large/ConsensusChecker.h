//
// Created by Bernardo Clavijo (TGAC) on 23/12/2016.
//

#ifndef W2RAP_CONTIGGER_CONSENSUSCHECKER_H
#define W2RAP_CONTIGGER_CONSENSUSCHECKER_H


#include <paths/HyperBasevector.h>
#include <paths/long/ReadPath.h>

struct read_call {
    uint64_t A;
    uint64_t C;
    uint64_t G;
    uint64_t T;
    uint64_t N;
    uint64_t total(){ return A+C+G+T+N;};
    double support(uint8_t base){
        if (base==0) return ((double) A)/total();
        if (base==1) return ((double) C)/total();
        if (base==2) return ((double) G)/total();
        if (base==3) return ((double) T)/total();
        if (base==4) return ((double) N)/total();
    }
} ;

class ConsensusChecker {
public:
    ConsensusChecker( HyperBasevector& hbv, vec<int>& inv, ReadPathVec& paths, vecbasevector& bases) :
    mHBV(hbv),
    mInv(inv),
    mPaths(paths),
    mBases(bases)
    {

    };

    bool consensus_OK (uint64_t edge);
    void collect_read_calls(uint64_t edge);
    void collect_edge_calls(uint64_t edge);
    vecbasevector generate_alternative_consensus(uint64_t edge);
    void print_detail();

private:
    HyperBasevector &mHBV;
    vec<int> &mInv;
    ReadPathVec &mPaths;
    vecbasevector &mBases;
    std::vector<uint8_t> edge_call;
    std::vector<struct read_call> read_calls;
};


#endif //W2RAP_CONTIGGER_CONSENSUSCHECKER_H
