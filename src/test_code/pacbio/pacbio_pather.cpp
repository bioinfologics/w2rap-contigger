//
// Created by Gonzalo Garcia (TGAC) on 04/11/2016.
//

#include "pacbio_pather.h"

PacbioPather::PacbioPather(vecbvec* aseqVector, HyperBasevector* ahbv) : KMatch(31) {
  //
  vecbvec* seqVector = seqVector;
  HyperBasevector* hbv = hbv;
};

void PacbioPather::mapReads(){
  // get the reads and map them to the graph using the dictionary
  // returns a vector of paths

  std::vector<linkReg> links = this->getReadsLinks(true);
}

//std::vector<linkReg> PacbioPather::getReadsLinks(vecbvec* seqVector, HyperBasevector* hbv, bool output_to_file=true){
std::vector<linkReg> PacbioPather::getReadsLinks(bool output_to_file=true){
  // Get the reads
  std::vector<linkReg> links;

  std::cout<< "Size of the dictionary: " << edgeMap.size() << std::endl;
  std::ofstream fout;
  fout.open("/Users/ggarcia/Documents/test_dataset/test/testlinks.txt");
  fout << "read readlength edge_id inv_edge_id eoffset kmer roffset" <<std::endl;

  auto edges = hbv->Edges();
  vec<int> inv;
  hbv->Involution(inv);

  int cont = 0;
  for (auto v=0; v<this->seqVector->size(); ++v){
    std::string read = (*seqVector)[v].ToString();
    auto g = this->lookupRead(read);
    if (g.size()>0){
      for (size_t a=0; a<g.size(); ++a){
        linkReg temp_link;
        temp_link.read_id = cont;
        temp_link.read_size = read.size();
        temp_link.read_offset = g[a].read_offset;
        temp_link.edge_id = g[a].edge_id;
        temp_link.edge_offset = g[a].edge_offset;
        temp_link.inv_edge_id = inv[g[a].edge_id];
        temp_link.kmer = g[a].kmer;
        fout << cont << " " << read.size() << " " << g[a].edge_id << " " << inv[g[a].edge_id] << " " << g[a].edge_offset<< " " << g[a].kmer << " " << g[a].read_offset << std::endl;
        links.push_back(temp_link);
      }
    }
    cont ++;
  }
  std::cout << "mapped reads: " << cont << std::endl;
  fout.close();
  return links;
}

std::vector<linkReg> PacbioPather::readLinksFilter(std::vector<linkReg> data, int read_id){
  // filter the links that come from a specific read id
  std::vector<linkReg> filtered;
  for (auto l: data){
      if (l.read_id == read_id){
        filtered.push_back(l);
      }
  }
  return filtered;
}

std::vector<linkReg> PacbioPather::readOffsetFilter(std::vector<linkReg> data){
  // Count the number of edges_id that map to each kmer in the read, filter out the edges that map to a position that has more than one edge mapped to it
  // input is a vector with all the links for a specific read

  std::vector<linkReg> links;
  std::map<int, unsigned int> pos_count;

  for (auto l: data){
      // Aumentar el contador, valores son inicializados a 0!?
      pos_count[l.read_offset] += 1;
  }

  // Solo seleccionados lo que estan en offsets de la lectura con un solo edge mapeado
  for (auto l: data){
      if (pos_count[l.read_offset] == 1){
        links.push_back(l);
      }
  }
  return links;
}