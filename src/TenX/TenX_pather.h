//
// Created by Gonzalo Garcia (TGAC) on 16/12/2016.
//

#ifndef W2RAP_CONTIGGER_TENX_PATHER_H
#define W2RAP_CONTIGGER_TENX_PATHER_H


#include "paths/local/LocalPather.h"

typedef  struct {
    int read_id;
//    int read_size;
    int edge_id;
    int inv_edge_id;
    int edge_offset;
    int read_offset;
    int kmer;
} tenXLink;

struct linkreg_less_than {
    inline bool operator() (const tenXLink& struct1, const tenXLink& struct2)
    {
      return (struct1.read_offset < struct2.read_offset);
    }
};



class TenXPather: public KMatch {
public:
    typedef std::uint16_t tagktype;
    std::map<uint64_t, std::map<tagktype, int>> kmerTagMap;

    TenXPather(std::vector<tenXRead>* aseqVector, HyperBasevector* ahbv);

    int createEmptyMap(HyperBasevector* hbv);
    int reads2kmerTagMap();
    float kmerTagDensity();

    typedef std::pair<std::vector<tenXLink>, std::vector<tenXLink>> pairLink; // A pairLinks are the links of r1,r2 of one tenxread
    typedef std::vector<pairLink> tagLink; // A tagLink are all the pairLink for a particular tag

    // Map reads to graph
    std::vector<TenXPather::tagktype> getSequenceTags(std::string seq);
    float edgeTagIntersection(std::string edgeFrom, std::string edgeTo, int roi);

//    int resolve_regions(int large_frontier_size=500);

private:
    // Reads and graph
    std::vector<tenXRead>* seqVector;
    vec<int> inv;
    HyperBasevector* hbv;
    tagktype kmerize_tag(std::string tag);

};

class LocalPaths_TX: public LocalPaths {
public:
    LocalPaths_TX(HyperBasevector* hbv, std::vector<std::vector<uint64_t>> pair_solutions, vec<int>& to_right, TenXPather* txp, std::vector<BaseVec>& edges)
        : LocalPaths(hbv, pair_solutions, to_right, edges), mTxp (txp) {};

    std::vector<uint64_t> choose_best_path(std::vector<std::vector<uint64_t>> *alternative_paths);

    TenXPather* mTxp;

};

#endif //W2RAP_CONTIGGER_TENX_PATHER_H
