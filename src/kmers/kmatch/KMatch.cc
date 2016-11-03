#include "KMatch.h"
#include <sys/time.h>
#include <thread>
#include <paths/HyperBasevector.h>
//#include <functional>

void timed_log(std::string s){
  struct timeval tp;
  gettimeofday(&tp,NULL);
  long int ms = tp.tv_sec * 1000 + tp.tv_usec / 1000;
  std::cout<<"TIME_LOG: "<<ms<<" - "<<s<<std::endl;
};

KMatch::KMatch(int kv){
  if (kv > 31){
    std::cout << "Kmer value is too big for this, using 31 instead... " << std::endl;
    this->K = 31;
  }
  this->K = kv;
}

std::vector<std::pair<uint_least64_t, int>> KMatch::ProduceKmers(std::string seq){
  // get a sequence a produce the set of kmers ()
  std::vector<std::pair<uint64_t, int>> kmer_vector;

  uint64_t kmer_index=0;
  const uint64_t KMER_MASK=( ((uint64_t)1)<<(this->K*2) )-1;
//  const uint64_t KMER_FIRSTOFFSET=(K-1)*2;

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
      //std::cout<<"c="<<s[p]<<" f="<<fkmer<<" r="<<rkmer<<std::endl;
      //TODO: last unknown passed by?
      if (last_unknown + this->K <= p) {
        auto pair_temp = std::make_pair(fkmer, kmer_index);
        kmer_vector.push_back(pair_temp);
        kmer_index++;
      }
    }
  } else {
    std::cout << "Sequence is too short to get kmers..." << std::endl;
  }

  return kmer_vector;
}

std::map<uint64_t, std::vector<std::pair<int, int>>> KMatch::Hbv2Map(HyperBasevector* hbv){
//  std::vector<kmer_position_t> karray;
  //read fasta and push_back(kmer,pos) (use pos as chr*CHR_CONST+offset)
  //open file

  std::map<uint64_t, std::vector<std::pair<int, int>>> edgeDict;
  uint32_t seq_index=0;

  auto edges = hbv->Edges();


  for (auto seqN=0; seqN<edges.size(); ++seqN) {
    auto seq = edges[seqN].ToString();
    auto kv = this->ProduceKmers(seq);

    for (auto &a: kv){
      if (edgeDict.find(a.first) == edgeDict.end()){
        std::vector<std::pair<int, int>> temp_vector;
        temp_vector.push_back(std::make_pair(seq_index, a.second));
        edgeDict[a.first] = temp_vector;
      } else {
        auto temp_vector = edgeDict[a.first];
        temp_vector.push_back(std::make_pair(seq_index, a.second));
        edgeDict[a.first] = temp_vector;
      }
    }
    seq_index++;
  }
  return edgeDict;
}

void KMatch::MapReads(vecbvec& seqVector){
  // get the reads and map them to the graph using the dictionary
}
