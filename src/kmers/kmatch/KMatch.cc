#include "KMatch.h"
#include <sys/time.h>
#include <thread>
#include <paths/HyperBasevector.h>

KMatch::KMatch(const int kv=31){
  if (kv > 31){
    std::cout << "Kmer value is too big for this, using 31 instead... " << std::endl;
    K = 31;
  }
  K = kv;
}

std::vector<pKmer> KMatch::ProduceKmers(const std::string &seq) const {
  // get a sequence a produce the set of kmers ()
  std::vector<pKmer> kmer_vector;

  int offset=0;
  const uint64_t KMER_MASK=( ((uint64_t)1)<<(K*2) )-1;

  if (seq.size()>K) {
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
      if (last_unknown + K <= p) {
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

void KMatch::Hbv2Map(const HyperBasevector *hbv){
  /* Takes the HBV and produces a map of kmer:edge
   * Each kmer is matched with a edgeKmerPostion object that stores extra information like:
   * - edge id
   * - position of the kmer in the edge
   * - position of the kmer in the read
   * - Kmer sequence
   * */

  //  std::vector<kmer_position_t> karray;

  std::map<uint64_t, std::vector<edgeKmerPosition>> edgeDict;
  uint32_t seq_index=0;

  auto edges = hbv->Edges();

  for (auto seqN=0; seqN<edges.size(); ++seqN) {
    const auto seq = edges[seqN].ToString();
    const auto kv (ProduceKmers(seq));

    for (auto a=0; a<kv.size(); ++a){
      if (edgeMap.find(kv[a].kmer) == edgeMap.end()){ // Not there, add the map entry (TODO: fix this, looks too complicated)
        std::vector<edgeKmerPosition> temp_vector;
        edgeKmerPosition tmatch;
        tmatch.edge_id = seq_index;
        tmatch.edge_offset = kv[a].offset;
        temp_vector.push_back(tmatch);
        edgeMap[kv[a].kmer] = temp_vector;

      } else {                                      // is there, first get the current content and then reattach (TODO: fix this, looks too complicated)
        auto temp_vector = edgeMap[kv[a].kmer];
        edgeKmerPosition tmatch;
        tmatch.edge_id = seq_index;
        tmatch.edge_offset = kv[a].offset;
        temp_vector.push_back(tmatch);
        edgeMap[kv[a].kmer] = temp_vector;
      }
    }
    seq_index++;
  }
}

std::vector<edgeKmerPosition> KMatch::lookupRead(const std::string &read){
  /* Find the kmers on the reads in the edges of the graph.
   * Each match produces a edgeKmerPosition object that stores the info of the match
   * returns all matches, one read can have multiple edges matching at this point
   * */

  // produce kmers
  const auto rkms (ProduceKmers(read));

  // look kmers in the dictionary
  std::vector<edgeKmerPosition> mapped_edges;
  int cont = 0; // Cont to hold the kmer offset in the read
  for (const auto &a: rkms){
    std::map<uint64_t, std::vector<edgeKmerPosition>>::iterator tt = edgeMap.find(a.kmer);
    if (tt != edgeMap.end()){
      for (const auto &p: tt->second){
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
}