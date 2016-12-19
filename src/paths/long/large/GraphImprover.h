//
// Created by Bernardo Clavijo (TGAC) on 14/12/2016.
//

#ifndef W2RAP_CONTIGGER_GRAPHIMPROVER_H
#define W2RAP_CONTIGGER_GRAPHIMPROVER_H

#include <paths/HyperBasevector.h>
#include <paths/long/ReadPath.h>
#include "OverlapValidator.h"

class GraphImprover {
public:
    GraphImprover( HyperBasevector& hbv, vec<int>& inv, ReadPathVec& paths, VecULongVec& invpaths) :
    mHBV(hbv),
    mInv(inv),
    mPaths(paths),
    mEdgeToPathIds(invpaths),
    mOverlapValidator(OverlapValidator(hbv,inv,paths)){};


    void improve_graph();

    //Graph modifying functions
    std::set<uint64_t> expand_cannonical_repeats(uint64_t min_support, uint64_t min_alternative_support);

    //Basic Graph Modifications

    std::map<uint64_t,std::vector<uint64_t>> separate_path(std::vector<uint64_t> p);
    void migrate_readpaths(std::map<uint64_t,std::vector<uint64_t>> edgemap);


private:
    HyperBasevector &mHBV;
    vec<int> &mInv;
    ReadPathVec &mPaths;
    OverlapValidator mOverlapValidator;
    VecULongVec& mEdgeToPathIds;
    vec<int> mToLeft;
    vec<int> mToRight;
};


#endif //W2RAP_CONTIGGER_GRAPHIMPROVER_H
