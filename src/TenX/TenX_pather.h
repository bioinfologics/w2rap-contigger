//
// Created by Gonzalo Garcia (TGAC) on 16/12/2016.
//

#ifndef W2RAP_CONTIGGER_TENX_PATHER_H
#define W2RAP_CONTIGGER_TENX_PATHER_H


#include "paths/local/LocalPather.h"
#include "paths/PathFinder.h"

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


class TenXPather: public KMatch, public PathFinder {
public:
    // Data types
    typedef std::uint16_t tagktype;
    std::map<uint64_t, std::map<tagktype, int>> kmerTagMap;

    TenXPather(std::vector<tenXRead>& aseqVector, HyperBasevector& ahbv, vec<int>& ainv, int min_reads, std::vector<BaseVec>& edges, ReadPathVec& apaths, VecULongVec& ainvPaths);

    int createEmptyMap(HyperBasevector* hbv);
    int reads2kmerTagMap();
    float kmerTagDensity();

    typedef std::pair<std::vector<tenXLink>, std::vector<tenXLink>> pairLink; // A pairLinks are the links of r1,r2 of one tenxread
    typedef std::vector<pairLink> tagLink; // A tagLink are all the pairLink for a particular tag

    // Map reads to graph
    std::vector<TenXPather::tagktype> getSequenceTags(std::string seq);
    float edgeTagIntersection(std::string edgeFrom, std::string edgeTo, int roi);

    // Pathfinder
    void solve_region_using_TenX(uint64_t large_frontier_size, bool verbose_separation=false, float score_threshold=1.0);

    std::vector<uint64_t> choose_best_path(std::vector<std::vector<uint64_t>> *alternative_paths);

private:
    // Reads and graph
    std::vector<tenXRead>& seqVector;
    tagktype kmerize_tag(std::string tag);
    std::vector<BaseVec>& mEdges;
};

#endif //W2RAP_CONTIGGER_TENX_PATHER_H
