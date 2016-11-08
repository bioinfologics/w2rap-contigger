//
// Created by Gonzalo Garcia (TGAC) on 04/11/2016.
//

#include "pacbio_pather.h"

PacbioPather::PacbioPather() : KMatch(31) {
};

void PacbioPather::mapReads(vecbvec seqVector, HyperBasevector* hbv){
  // get the reads and map them to the graph using the dictionary
  // returns a vector of paths
  std::cout<< "Size of the dictionary: " << edgeMap.size() << std::endl;
  std::ofstream fout;
  fout.open("/Users/ggarcia/Documents/test_dataset/test/testlinks.txt");
  fout << "read edge_id inv_edge_id offset kmer" <<std::endl;
  auto edges = hbv->Edges();
  vec<int> inv;
  hbv->Involution(inv);

  int cont = 0;
  for (auto v=0; v<seqVector.size(); ++v){
    std::vector<int> s;
    auto g = this->lookupRead(seqVector[v].ToString());

    if (g.size()>10){
      for (auto a=0; a<g.size(); ++a){
        s.push_back(g[a].edge_id);
//        std::cout << "Read: " << cont << " mapped "<< this->lookupRead(seqVector[v].ToString()).size() << " places*kmers"  << std::endl;
//        std::cout << cont << " " << g[a].edge_id << " " << inv[g[a].edge_id] << " " << g[a].offset<< " " << g[a].kmer << std::endl;
        fout << cont << " " << g[a].edge_id << " " << inv[g[a].edge_id] << " " << g[a].offset<< " " << g[a].kmer << std::endl;
      }
//      std::set<int> ss(s.begin(), s.end());
//      std::cout << "Read mapped: " << seqVector[v] << std::endl;
//      if (ss.size()>1){
//        for (auto x: ss){
//          std::cout << "Edge index: " << x << " " << inv[x] << std::endl;
//          //        std::cout << edges[x] << std::endl;
//        }
//      }
    }
    cont ++;
  }
  std::cout << "mapped reads: " << cont << std::endl;
  fout.close();
}
