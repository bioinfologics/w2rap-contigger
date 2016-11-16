//
// Created by Gonzalo Garcia (TGAC) on 04/11/2016.
//

#include "test_code/pacbio/pacbio_pather.h"


PacbioPather::PacbioPather(vecbvec* aseqVector, HyperBasevector* ahbv) : KMatch(31) {
  //
  seqVector = aseqVector;
  hbv = ahbv;
};

//std::vector<linkReg> PacbioPather::getReadsLinks(vecbvec* seqVector, HyperBasevector* hbv, bool output_to_file=true){
std::vector<linkReg> PacbioPather::getReadsLinks(bool output_to_file=true){
  // Get the reads
  std::vector<linkReg> links;

  std::cout<< "Size of the dictionary: " << edgeMap.size() << std::endl;
  std::ofstream fout;
  fout.open("/Users/ggarcia/Documents/test_dataset/test/testlinks.txt");
  fout << "read readlength edge_id inv_edge_id eoffset kmer roffset" <<std::endl;

  std::cout << "loading edges and involution." << std::endl;
  auto edges = hbv->Edges();
  std::cout << "edges loaded" << std::endl;
  vec<int> inv;
  hbv->Involution(inv);

  std::cout << "processing edges" << std::endl;
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



ReadPathVec PacbioPather::mapReads(){
  // get the reads and map them to the graph using the dictionary
  // returns a vector of paths
  std::cout << "Executing getReadsLines.." << std::endl;
//  std::vector<std::vector<int>> pb_paths;
  std::vector<linkReg> links = this->getReadsLinks(true);
  std::cout << "Done..." << std::endl;

  // To store the paths
  ReadPathVec pb_paths;

  // for each read
  for (std::uint32_t r=0; r<this->seqVector->size(); ++r){

    // filter the links for this reads
    auto read_links = this->readLinksFilter(links, r);

    // filter the shared roffsets
    auto offset_filter = this->readOffsetFilter(read_links);

    // filter the min match length TODO: this function

    // sort the vector
    std::sort(offset_filter.begin(), offset_filter.end(), linkreg_less_than());

    // Filter if the reverse complement (?)

    // Create vector of unique edge_ids
    std::vector<int> presentes;
    std::vector<linkReg> s_edges;
    for (auto a: offset_filter){
      if (std::find(presentes.begin(), presentes.end(), a.edge_id) == presentes.end() && std::find(presentes.begin(), presentes.end(), a.inv_edge_id) == presentes.end()){
        s_edges.push_back(a);
        presentes.push_back(a.edge_id);
      }
    }

    std::vector<int> temp_path;
    if (s_edges.size() >= 2) {
     std::cout << "-----------------------" << std::endl;
      std::cout << "read: " << r << std::endl;
      std::cout << "-----------------------" << std::endl;
      int poffset=s_edges[0].edge_offset;
      for (auto s: s_edges) {
        std::cout << s.edge_id << "," << s.inv_edge_id << std::endl;
        temp_path.push_back(s.edge_id);
      }
      // apend the path to the vector
      ReadPath tp(poffset, temp_path);
      pb_paths.push_back(tp);
    }
  }

  return pb_paths;
}


//void WriteReadPathVec(const ReadPathVec &rpv, const char * filename){
//  std::ofstream f(filename, std::ios::out | std::ios::trunc | std::ios::binary);
//  uint64_t pathcount=rpv.size();
//  f.write((const char *) &pathcount, sizeof(pathcount));
//  uint16_t ps;
//  int mOffset;
//  for (auto const &rp:rpv){
//    mOffset=rp.getOffset();
//    ps=rp.size();
//    f.write((const char *) &mOffset, sizeof(mOffset));
//    f.write((const char *) &ps, sizeof(ps));
//    f.write((const char *) rp.data(),ps*sizeof(int));
//  }
//  f.close();
//}