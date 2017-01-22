//
// Created by Bernardo Clavijo (TGAC) on 11/07/2016.
//

#ifndef W2RAP_CONTIGGER_PATHFINDER_H_TX
#define W2RAP_CONTIGGER_PATHFINDER_H_TX

#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "paths/long/large/GapToyTools.h"
#include "Vec.h"
#include "TenX/TenX_pather.h"
#include "paths/local/LocalPather.h"
#include "paths/PathFinder.h"

//class PathFinder_tx {
class PathFinder_tx: public PathFinder {
public:
    PathFinder_tx (TenXPather& txp, HyperBasevector& hbv, vec<int> inv, int min_reads, ReadPathVec& paths, VecULongVec& invPaths);
    void solve_region_using_TenX(uint64_t large_frontier_size, bool verbose_separation=false, float score_threshold=1.0);


//private:
    TenXPather& mTxp;
};

#endif //W2RAP_CONTIGGER_PATHFINDER_H_TX
