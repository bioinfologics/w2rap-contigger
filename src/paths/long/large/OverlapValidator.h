//
// Created by Bernardo Clavijo (TGAC) on 12/12/2016.
//

#ifndef W2RAP_CONTIGGER_OVERLAPVALIDATOR_H
#define W2RAP_CONTIGGER_OVERLAPVALIDATOR_H


#include <paths/HyperBasevector.h>
#include <paths/long/ReadPath.h>

class InformativePair{ //a path that has something to offer, we'll save them all together for convenience

public:
    InformativePair(uint64_t index, ReadPath & r1Path, ReadPath & r2Path, vec<int> &inv);
    inline bool is_combined(){return combined_path.size()>0;};
    std::vector<uint64_t> overlaps_crossed(vec<int> & toLeft, vec<int> &toRight);
    std::vector<uint64_t> overlaps_jumped(vec<int> & toLeft, vec<int> &toRight, uint8_t max_dist);
    bool crosses_transition(uint64_t e1, uint64_t e2);
    bool jumps_transition(uint64_t e1, uint64_t e2);

    //todo: add support for RF pairs
    uint64_t pathIndex; //original index of the first component of the path;
    std::vector<uint64_t> r1path,r1rpath,r2path,r2rpath,combined_path,combined_rpath;
};

class OverlapValidator {
public:
    OverlapValidator( HyperBasevector& hbv, vec<int>& inv, ReadPathVec& paths) :
    mHBV(hbv),
    mInv(inv),
    mPaths(paths)
    {};
    void find_informative_pairs();
    void compute_overlap_support();
    void analyse_complex_overlaps();
    std::vector<uint64_t> find_perfect_tips(uint16_t max_size);
    uint64_t collect_all_support(uint64_t vi, uint64_t e1, uint64_t e2);

private:
    HyperBasevector &mHBV;
    vec<int> &mInv;
    ReadPathVec &mPaths;
    std::vector<InformativePair> mInformativePairs;
    std::vector<std::vector<uint64_t>> mCross, mJump; //for each overlap


};


#endif //W2RAP_CONTIGGER_OVERLAPVALIDATOR_H
