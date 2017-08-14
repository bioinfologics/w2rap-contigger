//
// Created by Luis Yanes (EI) on 27/07/2017.
//

#ifndef W2RAP_CONTIGGER_KMERFREQFACTORY_H
#define W2RAP_CONTIGGER_KMERFREQFACTORY_H

#define likely(x)       __builtin_expect((x),1)
#define unlikely(x)     __builtin_expect((x),0)

struct KMerParams {
    uint8_t k;
    unsigned int minQual;
};


template<typename FileRecord>
class KMerFreqFactory{
public:
    KMerFreqFactory(KMerParams params) : K(params.k), minQual(params.minQual),
                                         KMER_FIRSTOFFSET((uint64_t) (K - 1) * 2),
                                         KMER_MASK((((uint64_t) 1) << (K * 2)) - 1),
                                         p(0), bases(0) {
    }

    ~KMerFreqFactory() {
        std::cout << "Bases processed " << bases << "\n";
    }
    void setFileRecord(FileRecord& rec){
        std::swap(rec,currentRecord);
        nRecords = currentRecord.size()-K;
        p=0;
    }

    const bool next_element(KMerNodeFreq_s& mer){
/*        if (unlikely(K > currentRecord.seq.size()) or currentRecord.qual[p] < minQual) return false;  */
        if (unlikely(p >= currentRecord.size())) return false;
        if (unlikely(p == 0)) {
            kkk = KMerNodeFreq(currentRecord.begin());
            kkk.hash();
            kkk.kc = KMerContext::initialContext(*(currentRecord.begin()+K));
            kkk.count = 1;
            bases+=K;
            p = K;
            (kkk.isRev()) ? KMerNodeFreq(kkk,true).to_struct(mer) : kkk.to_struct(mer);
            return p < currentRecord.size();
        }

        auto itr = currentRecord.cbegin()+p;
        while (itr != currentRecord.cend()-1) {
            bases++;
            unsigned char pred = kkk.front();
            kkk.toSuccessor(*itr);
            ++itr;
            kkk.kc = KMerContext(pred, *itr);
            p++;
            (kkk.isRev()) ? KMerNodeFreq(kkk,true).to_struct(mer) : kkk.to_struct(mer);
            return true;
        }
        kkk.kc = KMerContext::finalContext(kkk.front());
        kkk.toSuccessor(*(currentRecord.cbegin() + (currentRecord.size()-1)));
        (kkk.isRev()) ? KMerNodeFreq(kkk,true).to_struct(mer) : kkk.to_struct(mer);
        p++;
        return true;
    }

private:
    KMerNodeFreq kkk;
    FileRecord currentRecord;
    uint64_t p;
    uint64_t bases;
    uint nRecords;
    unsigned int minQual;
    const uint8_t K;
    const uint64_t KMER_MASK;
    const uint64_t KMER_FIRSTOFFSET;
};


#endif //W2RAP_CONTIGGER_KMERFREQFACTORY_H
