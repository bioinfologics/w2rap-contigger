#include "KMatch.h"
#include <sys/time.h>
#include <thread>
#include <paths/HyperBasevector.h>

KMatch::KMatch(int kv){
    if (kv > 31){
        std::cout << "Kmer value is too big for this, using 31 instead... " << std::endl;
        this->K = 31;
    }
    this->K = kv;
}

std::vector<pKmer> KMatch::ProduceKmers(std::string seq){
    // get a sequence a produce the set of kmers ()
    std::vector<pKmer> kmer_vector;
    kmer_vector.reserve(seq.size()-K+1);
    int offset=0;
    const uint64_t KMER_MASK=( ((uint64_t)1)<<(this->K*2) )-1;

    if (seq.size()>this->K) {
        const char *s = seq.c_str();
        int64_t last_unknown = -1;
        uint64_t fkmer = 0;
        for (auto p=0; p < seq.size(); ++p) {
            //fkmer: grows from the right (LSB)
            switch (s[p]) {
                case 'A':
                case 'a':
                    fkmer = ((fkmer << 2) + 0) & KMER_MASK;
                    break;
                case 'C':
                case 'c':
                    fkmer = ((fkmer << 2) + 1) & KMER_MASK;
                    break;
                case 'G':
                case 'g':
                    fkmer = ((fkmer << 2) + 2) & KMER_MASK;
                    break;
                case 'T':
                case 't':
                    fkmer = ((fkmer << 2) + 3) & KMER_MASK;
                    break;
                default:
                    fkmer = ((fkmer << 2) + 0) & KMER_MASK;
                    last_unknown = p;
                    break;
            }
            //TODO: last unknown passed by?
            if (last_unknown + this->K <= p) {
                pKmer pair_temp;
                pair_temp.kmer = fkmer;
                pair_temp.offset = offset;
                kmer_vector.push_back(pair_temp);
                offset++;
            }
        }
    }
    return kmer_vector;
}

void KMatch::Hbv2Map(HyperBasevector & hbv){
    //  std::vector<kmer_position_t> karray;



    auto edges = hbv.Edges();
    uint64_t kc=0;
    for (auto &e:edges) kc+=e.size()+1-K;
    //edgeMap.reserve(kc);
    edgeKmerPositionV tmatch;
    std::vector<edgeKmerPositionV> nv;
    for (auto seqN=0; seqN<edges.size(); ++seqN) {
        auto seq = edges[seqN].ToString();
        auto kv = this->ProduceKmers(seq);

        for (auto a=0; a<kv.size(); ++a){
            tmatch.edge_id = seqN;
            tmatch.edge_offset = kv[a].offset;
            nv.clear();
            nv.push_back(tmatch);
            auto const &r=edgeMap.insert(KMAP_t::value_type(kv[a].kmer,nv));
            if (not r.second){
                r.first->second.push_back(tmatch);
            }
            /*if (this->edgeMap.find(kv[a].kmer) == this->edgeMap.end()){
                std::vector<edgeKmerPositionV> temp_vector;
                tmatch.edge_id = seqN;
                tmatch.edge_offset = kv[a].offset;
                temp_vector.push_back(tmatch);
                this->edgeMap[kv[a].kmer] = temp_vector;
            } else {
                tmatch.edge_id = seqN;
                tmatch.edge_offset = kv[a].offset;
                this->edgeMap[kv[a].kmer].push_back(tmatch);
            }*/
        }
        if ((seqN+1)%25000==0)std::cout<<Date()<<": "<<seqN+1<<"edges processed, "<<edgeMap.size()<<" "<< (int) K<<"-mers "<<std::endl;
    }
}

void KMatch::Hbv2Index(HyperBasevector &hbv){
    //  std::vector<kmer_position_t> karray;

    KmerIndex.clear();

    auto edges = hbv.Edges();
    uint64_t kc=0;
    for (auto &e:edges) kc+=e.size()+1-K;
    KmerIndex.reserve(kc);
    std::cout<<Date()<<": creating kmer positions" <<std::endl;
    //map
    for (uint64_t seqN=0; seqN<edges.size(); ++seqN) {
        auto seq = edges[seqN].ToString();
        auto kv = this->ProduceKmers(seq);
        edgeKmerPositionNR tmatch;
        for (auto i=0; i<kv.size(); ++i){
            tmatch.kmer=kv[i].kmer;
            tmatch.edge_id = seqN;
            tmatch.edge_offset = kv[i].offset;
            KmerIndex.push_back(tmatch);
        }
        if ((seqN+1)%100000==0)std::cout<<Date()<<": "<<seqN+1<<" edges processed, "<<KmerIndex.size()<<" "<<K<<"-mers "<<std::endl;
    }
    std::cout<<Date()<<": sorting"<<std::endl;
    std::sort(KmerIndex.begin(),KmerIndex.end());
}
std::vector<edgeKmerPosition> KMatch::lookupRead(std::string read){
    // produce kmers
    auto rkms = this->ProduceKmers(read);

    // look kmers in the dictionary
    std::vector<edgeKmerPosition> mapped_edges;
    int cont = 0; // Cont to hold the kmer offset in the read
    for (auto a: rkms){
        //std::cout<<" finding kmer "<<cont<<std::endl;
        edgeKmerPositionNR xs;
        xs.kmer = a.kmer;
        edgeKmerPosition x;
        x.kmer = a.kmer;
        auto kis = std::lower_bound(KmerIndex.begin(),KmerIndex.end(),xs);

        while (kis != KmerIndex.end() and kis->kmer==x.kmer){
            x.edge_id = kis->edge_id;
            x.edge_offset = kis->edge_offset;
            x.read_offset = cont;
            mapped_edges.push_back(x);
            kis++;
        }

        cont++;
    }
    return mapped_edges;
}

/*std::vector<edgeKmerPosition> KMatch::lookupRead(std::string read){
    // produce kmers
    auto rkms = this->ProduceKmers(read);

    // look kmers in the dictionary
    std::vector<edgeKmerPosition> mapped_edges;
    int cont = 0; // Cont to hold the kmer offset in the read
    for (auto a: rkms){
        auto tt = this->edgeMap.find(a.kmer);
        if (tt != this->edgeMap.end()){
            for (auto p: tt->second){
                edgeKmerPosition x;
                x.edge_id = p.edge_id;
                x.edge_offset = p.edge_offset;
                x.read_offset = cont;
                x.kmer = a.kmer;
                mapped_edges.push_back(x);
            }
        }
        cont++;
    }
    return mapped_edges;
}*/