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

void KMatch::Hbv2Map(HyperBasevector* hbv){
  //  std::vector<kmer_position_t> karray;

  std::map<uint64_t, std::vector<edgeKmerPosition>> edgeDict;
  uint32_t seq_index=0;

  auto edges = hbv->Edges();

  for (auto seqN=0; seqN<edges.size(); ++seqN) {
    auto seq = edges[seqN].ToString();
    auto kv = this->ProduceKmers(seq);

    for (auto &a: kv){
      if (this->edgeMap.find(a.kmer) == this->edgeMap.end()){
        std::vector<edgeKmerPosition> temp_vector;
        edgeKmerPosition tmatch;
        tmatch.edge_id = seq_index;
        tmatch.offset = a.offset;
        temp_vector.push_back(tmatch);
        this->edgeMap[a.kmer] = temp_vector;

      } else {
        auto temp_vector = this->edgeMap[a.kmer];
        edgeKmerPosition tmatch;
        tmatch.edge_id = seq_index;
        tmatch.offset = a.offset;
        temp_vector.push_back(tmatch);
        this->edgeMap[a.kmer] = temp_vector;
      }
    }
    seq_index++;
  }
}

std::vector<edgeKmerPosition> KMatch::lookupRead(std::string read){
  // produce kmers
  auto rkms = this->ProduceKmers(read);

  // look kmers in the dictionary
  std::vector<edgeKmerPosition> mapped_edges;
  for (auto a: rkms){
    std::map<uint64_t, std::vector<edgeKmerPosition>>::iterator tt = this->edgeMap.find(a.kmer);
    if (tt != this->edgeMap.end()){
      for (auto p: tt->second){
        edgeKmerPosition x;
        x.edge_id = p.edge_id;
        x.offset = p.offset;
        x.kmer = a.kmer;
        mapped_edges.push_back(x);
      }
    }
  }
  return mapped_edges;
}