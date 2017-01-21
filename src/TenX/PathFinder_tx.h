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

class PathFinder_tx {
public:
    PathFinder_tx(TenXPather* txp, HyperBasevector* hbv, vec<int> inv, int min_reads = 5 );

    void untangle_complex_in_out_choices(uint64_t large_frontier_size, bool verbose_separation=false, float score_threshold=1.0);

    void init_prev_next_vectors();
    std::string path_str(std::vector<uint64_t> e);
    std::map<uint64_t,std::vector<uint64_t>> separate_path(std::vector<uint64_t> p, bool verbose_separation=false);
    std::array<std::vector<uint64_t>,2>  get_all_long_frontiers(uint64_t e,uint64_t large_frontier_size);
    void migrate_readpaths(std::map<uint64_t,std::vector<uint64_t>> edgemap);


private:
    TenXPather* mTxp;
    HyperBasevector* mHBV;
    vec<int> mInv;
    ReadPathVec& mPaths;
    VecULongVec& mEdgeToPathIds;
    vec<int> mToLeft;
    vec<int> mToRight;
    std::vector<std::vector<uint64_t>> next_edges,prev_edges;
    int mMinReads;

};

#endif //W2RAP_CONTIGGER_PATHFINDER_H_TX
