//
// Created by Gonzalo Garcia (TGAC) on 16/12/2016.
//

#ifndef W2RAP_CONTIGGER_TENX_PATHER_H
#define W2RAP_CONTIGGER_TENX_PATHER_H

#include "kmers/kmatch/KMatch.h"
#include "paths/long/ReadPath.h"
#include "paths/long/large/ExtractReads.h"

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

    int createEmptyMap(HyperBasevector* hbv);

    int reads2kmerTagMap();

    typedef std::pair<std::vector<tenXLink>, std::vector<tenXLink>> pairLink; // A pairLinks are the links of r1,r2 of one tenxread
    typedef std::vector<pairLink> tagLink; // A tagLink are all the pairLink for a particular tag

    TenXPather(std::vector<tenXRead>* aseqVector, HyperBasevector* ahbv);
//    ReadPathVec mapReads();

    // Qc the run
    std::vector<int> readsTagQc();

    //TODO: Map indexing records per tag and index
    //TODO: Filter the tags
    //TODO: Filter the reads

    // Map reads to graph
    std::vector<tagLink> getTagLinks(bool output_to_file=true);
    std::vector<TenXPather::tagktype> getSequenceTags(std::string seq);
    std::vector<TenXPather::tagktype> edgeTagIntersection(std::string edgeFrom, std::string edgeTo, int roi);

private:
    // Reads and graph
    std::vector<tenXRead>* seqVector;
    vec<int> inv;
    HyperBasevector* hbv;
    tagktype kmerize_tag(std::string tag);

    // Read processing data
//    std::vector<std::string> tagVector;
//    std::map<std::string, std::vector<int>> indexMap;
//    std::map<std::string, std::vector<int>> tagMap;

//    int makeIndexMap(bool to_disc=false);
//    int makeTagMap(bool to_disc=false);
//    int makeTagVector();
//    std::vector<tenXLink> processLinks(std::string read, int read_id);
};



#endif //W2RAP_CONTIGGER_TENX_PATHER_H
