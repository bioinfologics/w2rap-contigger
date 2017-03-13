//
// Created by Gonzalo Garcia (TGAC) on 04/11/2016.
//

#ifndef W2RAP_CONTIGGER_PACBIO_PATHER_H
#define W2RAP_CONTIGGER_PACBIO_PATHER_H

#include "kmers/kmatch/KMatch.h"
#include "paths/long/ReadPath.h"
#include "paths/PathFinder.h"

typedef struct {
    int read_id;
    int read_size;
    int edge_id;
    int inv_edge_id;
    int edge_offset;
    int read_offset;
    int kmer;
} linkReg;



class LongReadPather: public KMatch, public PathFinder {
public:
    LongReadPather(const vecbvec &aseqVector, HyperBasevector &ahbv, vec<int> &ainv, const int min_reads,
                   std::vector<BaseVec> &edges, ReadPathVec &apaths, VecULongVec &ainvPaths);
    ReadPathVec mapReads();

    std::vector<uint64_t> choose_best_path(std::vector<std::vector<uint64_t>>* alternative_paths){};

    void solve_using_long_read(uint64_t large_frontier_size, bool verbose_separation);
private:

    struct less_than {
        inline bool operator() (const linkReg& struct1, const linkReg& struct2) const
        {
            return (struct1.read_offset < struct2.read_offset);
        }
    };

    const vecbvec& seqVector;
    std::vector<BaseVec>& mEdges;

    std::vector<std::vector<linkReg>> getReadsLinks(bool output_to_file=true);
    std::vector<linkReg> readOffsetFilter(const vector<linkReg> &data) const;
    std::vector<linkReg> minCoverageFilter(const vector<linkReg> &data) const;
};


#endif //W2RAP_CONTIGGER_PACBIO_PATHER_H
