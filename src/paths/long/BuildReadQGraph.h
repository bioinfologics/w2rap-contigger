///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * BuildReadQGraph.h
 *
 *  Created on: Jan 22, 2014
 *      Author: tsharpe
 */

#ifndef PATHS_LONG_BUILDREADQGRAPH_H_
#define PATHS_LONG_BUILDREADQGRAPH_H_

#include "dvString.h"
#include "feudal/PQVec.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "kmers/ReadPather.h"
#include "kmers/KMer.h"
#include "kmers/KMerContext.h"

typedef struct __attribute__((__packed__)) KMerNodeFreq_s {
    uint64_t kdata[2];
    uint8_t count;
    uint8_t kc;
    inline const bool operator==(KMerNodeFreq_s const & other) const {
        //return 0==memcmp(&kdata,&other.kdata,2*sizeof(uint64_t));
        return (kdata[0]==other.kdata[0] and kdata[1]==other.kdata[1]);
    }
    inline const operator<(KMerNodeFreq_s const & other) const{
        //return -1==memcmp(&kdata,&other.kdata,2*sizeof(uint64_t));
        if (kdata[0]<other.kdata[0]) return true;
        if (kdata[0]>other.kdata[0]) return false;
        return kdata[1]<other.kdata[1];
    }
    inline const operator>(KMerNodeFreq_s const & other) const{
        if (kdata[0]>other.kdata[0]) return true;
        if (kdata[0]<other.kdata[0]) return false;
        return kdata[1]>other.kdata[1];
    }
    inline void combine(KMerNodeFreq_s const & other){
        uint16_t newcount=count+other.count;
        count = (newcount > UINT8_MAX) ? UINT8_MAX : newcount;
        kc|=other.kc;
    }
    inline void merge(KMerNodeFreq_s const & other){
        uint16_t newcount=count+other.count;
        count = (newcount > UINT8_MAX) ? UINT8_MAX : newcount;
        kc|=other.kc;
    }
    /*inline const operator=(KMerNodeFreq_s const & other) const {
        memcpy(&this,&other,sizeof(this));
    }*/
    friend std::ostream& operator<<(std::ostream& os, const KMerNodeFreq_s& kmer);
    friend std::istream& operator>>(std::istream& is, KMerNodeFreq_s& kmer);

};

const unsigned K = 60;
typedef KMer<K> BRQ_Kmer;
typedef KMer<K-1> BRQ_SubKmer;
typedef KmerDictEntry<K> BRQ_Entry;
typedef KmerDict<K> BRQ_Dict;
//we can reduce this usage by storing blocks of kmers and using a difference with the previous. we would need slightly more complex heuristics.


class KMerNodeFreq: public KMer<K>{
public:
    KMerNodeFreq(){};
    template <class Itr>
    explicit KMerNodeFreq( Itr start )
    { assign(start,NopMapper()); }
    template <class Itr>
    explicit KMerNodeFreq( Itr start, bool dummy)
    { assign(start, CharToBaseMapper());}

    KMerNodeFreq (const KMerNodeFreq &other){
        *this=other;
        count=other.count;
        kc=other.kc;
    }
    KMerNodeFreq (const KMerNodeFreq_s &other) {
        memcpy (&this->mVal,&other.kdata,sizeof(KMer<K>));
        count=other.count;
        kc.mVal=other.kc;
    }

    KMerNodeFreq (const KMerNodeFreq &other, bool rc){
        *this=other;
        count=other.count;
        kc=other.kc;
        if (rc){
            this->rc();
            kc=kc.rc();
        }
    }
    void to_struct (KMerNodeFreq_s &other) const {
        memcpy (&other.kdata,&this->mVal,sizeof(KMer<K>));
        other.count=count;
        other.kc=kc.mVal;
    }
    unsigned char count;

    KMerContext kc;
};

class KmerList{
public:
    ~KmerList();
    void clear();
    void merge(KmerList & other);
    void sort();
    void uniq();
    void dump(std::string filename);
    void load(std::string filename);
    void resize(size_t new_size);
    KMerNodeFreq_s * kmers = nullptr;
    size_t size=0;
};

void create_read_lengths(std::vector<uint16_t> & rlen, VecPQVec const& quals, unsigned minQual);

void buildReadQGraph( vecbvec const & reads, VecPQVec const &quals, std::shared_ptr<KmerList> kmerlist,
                      bool doFillGaps, bool doJoinOverlaps,
                      unsigned minFreq, double minFreq2Fract, unsigned maxGapSize,  HyperBasevector* pHBV,
                      ReadPathVec* pPaths, int _K);

std::shared_ptr<KmerList> buildKMerCount( vecbvec const& reads,
                                                             std::vector<uint16_t> & rlen, unsigned minCount,
                                                             std::string workdir, std::string tmpdir,
                                                             unsigned char disk_batches, uint64_t count_batch_size );


#endif /* PATHS_LONG_BUILDREADQGRAPH_H_ */
