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
    GraphImprover( HyperBasevector& hbv, vec<int>& inv, ReadPathVec& paths) :
    mHBV(hbv),
    mInv(inv),
    mPaths(paths),
    mOverlapValidator(OverlapValidator(hbv,inv,paths)){};
    void improve_graph();
private:
    HyperBasevector &mHBV;
    vec<int> &mInv;
    ReadPathVec &mPaths;
    OverlapValidator mOverlapValidator;
};


#endif //W2RAP_CONTIGGER_GRAPHIMPROVER_H
