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

struct linkreg_less_than_pb {
    inline bool operator() (const linkReg& struct1, const linkReg& struct2)
    {
      return (struct1.read_offset < struct2.read_offset);
    }
};

class PacbioPather: public KMatch, public PathFinder {
public:
    PacbioPather(vecbvec& aseqVector, HyperBasevector& ahbv, vec<int>& ainv, int min_reads, std::vector<BaseVec>& edges, ReadPathVec& apaths, VecULongVec& ainvPaths);
    ReadPathVec mapReads();

    std::vector<uint64_t> choose_best_path(std::vector<std::vector<uint64_t>>* alternative_paths){};

    void PacbioPather::solve_using_long_read(uint64_t large_frontier_size, bool verbose_separation);
private:
    vecbvec& seqVector;
    std::vector<BaseVec>& mEdges;

    std::vector<std::vector<linkReg>> getReadsLinks(bool output_to_file=true);
    std::vector<linkReg> readOffsetFilter(std::vector<linkReg> data);
    std::vector<linkReg> matchLengthFilter(std::vector<linkReg> data);
};


#endif //W2RAP_CONTIGGER_PACBIO_PATHER_H
