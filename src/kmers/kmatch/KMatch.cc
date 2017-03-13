#include "KMatch.h"
#include <sys/time.h>
#include <thread>
#include <unordered_map>
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

void KMatch::Hbv2Map(const HyperBasevector &hbv){
  /* Takes the HBV and produces a map of kmer:edgeKmerPostion (is the position in the edge)
   * Each kmer is matched with a edgeKmerPostion object that stores extra information like:
   * - edge id
   * - position of the kmer in the edge
   * - position of the kmer in the read
   * - Kmer sequence
   * */

  //  std::vector<kmer_position_t> karray;

  std::unordered_map<uint64_t, std::vector<edgeKmerPosition>> edgeDict;
  uint32_t seq_index=0;

  const auto &edges = hbv.Edges();

  for (auto seqN=0; seqN<edges.size(); ++seqN) {
    const auto seq = edges[seqN].ToString();
    const auto kmer_vector (ProduceKmers(seq));

    for (const auto & kmer: kmer_vector){
      if (edgeMap.find(kmer.kmer) == edgeMap.cend()){
        std::vector<edgeKmerPosition> temp_vector;
        temp_vector.push_back(edgeKmerPosition (kmer.kmer, seq_index, kmer.offset));
        edgeMap[kmer.kmer] = temp_vector;
      } else {
        edgeMap[kmer.kmer].push_back(edgeKmerPosition (kmer.kmer, seq_index, kmer.offset));
      }
    }
    seq_index++;
  }
  std::cout << Date() << ": Done creating the hbv map. Keys: " << edgeMap.size() << std::endl;
}

std::vector<edgeKmerPosition> KMatch::lookupRead(const std::string &read) const{
  /* Find the kmers on the reads in the edges of the graph.
   * Each match produces a edgeKmerPosition object that stores the info of the match
   * returns all matches, one read can have multiple edges matching at this point
   * */

  // produce kmers
  const auto read_kmers (ProduceKmers(read));

  // look kmers in the dictionary
  std::vector<edgeKmerPosition> mapped_edges;
  int cont = 0; // Cont to hold the kmer offset in the read
  for (const auto &kmer: read_kmers){
    const auto kmer_in_graph = edgeMap.find(kmer.kmer);
    if (kmer_in_graph != edgeMap.end()){
      for (const auto &p: kmer_in_graph->second){
        mapped_edges.push_back(edgeKmerPosition (kmer.kmer, p.edge_id, p.edge_offset, cont));
      }
    }
    cont++;
  }
  return mapped_edges;
}