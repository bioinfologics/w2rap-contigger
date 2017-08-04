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
        memset(b2f, 0, 255);
        b2f['a'] = b2f['A'] = 0;
        b2f['c'] = b2f['C'] = 1;
        b2f['g'] = b2f['G'] = 2;
        b2f['t'] = b2f['T'] = 3;
    }

    ~KMerFreqFactory() {
        std::cout << "Bases processed " << bases << "\n";
    }
    void setFileRecord(FileRecord& rec){
        std::swap(rec,currentRecord);
        p=0;
    }

    const bool next_element(KMerNodeFreq_s& mer){
/*        if (unlikely(K > currentRecord.seq.size()) or currentRecord.qual[p] < minQual) return false;  */
        if (unlikely(K > currentRecord.size())) return false;
        if (unlikely(p == 0)) {
            kkk = KMerNodeFreq(currentRecord.begin());

            kkk.hash();
/*            kkk.kc = KMerContext::initialContext(b2f[*(currentRecord.seq.begin()+K)]);*/
            kkk.kc = KMerContext::initialContext(*(currentRecord.begin()+K));
            kkk.count = 1;
            p = K;
            (kkk.isRev()) ? KMerNodeFreq(kkk,true).to_struct(mer) : kkk.to_struct(mer);
/*            for (int i = 0; i < K; ++i) if (currentRecord.qual[i] < minQual) return false;*/
            /*return p < currentRecord.seq.size();*/
            return p < currentRecord.size();
        }
        /*auto itr = currentRecord.seq.cbegin()+p;*/
        auto itr = currentRecord.cbegin()+p;
        /*while (itr != currentRecord.seq.cend() and currentRecord.qual[p] > minQual) {*/
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
        kkk.toSuccessor(*itr);
        (kkk.isRev()) ? KMerNodeFreq(kkk,true).to_struct(mer) : kkk.to_struct(mer);
        return true;
    }

private:
    KMerNodeFreq kkk;
    FileRecord currentRecord;
    uint64_t p;
    uint64_t bases;
    unsigned int minQual;
    const uint8_t K;
    const uint64_t KMER_MASK;
    const uint64_t KMER_FIRSTOFFSET;
    unsigned char b2f[256];
};


#endif //W2RAP_CONTIGGER_KMERFREQFACTORY_H
