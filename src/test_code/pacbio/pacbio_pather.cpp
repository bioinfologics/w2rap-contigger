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
std::vector<std::vector<linkReg>> PacbioPather::getReadsLinks(bool output_to_file=true){
  // Get the reads
  std::vector<std::vector<linkReg>> links;

  std::cout<< Date() << ": Size of the dictionary: " << edgeMap.size() << std::endl;

  std::cout<< Date() << ": loading edges and involution." << std::endl;
  auto edges = hbv->Edges();
  std::cout << "edges loaded" << std::endl;
  vec<int> inv;
  hbv->Involution(inv);

  std::cout<< Date() << ": processing edges" << std::endl;
//  int cont = 0;
  links.resize(seqVector->size());
#pragma omp parallel for
  for (auto v=0; v<seqVector->size(); ++v){
    std::string read = (*seqVector)[v].ToString();
    auto g = lookupRead(read);
    if (g.size()>0){
      for (size_t a=0; a<g.size(); ++a){
        linkReg temp_link;
        temp_link.read_id = v;
        temp_link.read_size = read.size();
        temp_link.read_offset = g[a].read_offset;
        temp_link.edge_id = g[a].edge_id;
        temp_link.edge_offset = g[a].edge_offset;
        temp_link.inv_edge_id = inv[g[a].edge_id];
        temp_link.kmer = g[a].kmer;
#pragma omp critical
        links[v].push_back(temp_link);
      }
    }
  }

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

std::vector<linkReg> PacbioPather::matchLengthFilter(std::vector<linkReg> data){
  // Filter matches by length of the match

  std::map<std::string, int> index_map;
  for (auto l: data){
    // juntar la lectura y el eje en el mismo string concatenando asi pueod indexar por los dos
    std::string key = l.read_id + "-" +l.edge_id;
    index_map[key] += 1;
  }

  // contar los links por lectura y por edge
  std::vector<std::string> filtered_links;
  for (std::map<std::string, int>::iterator k=index_map.begin(); k != index_map.end(); ++k){
    //
    if (k->second > 10) {
      filtered_links.push_back(k->first);
    }
  }

  // Filtrar los que no pegan mas del limite
  std::vector<linkReg> good_links;
  for (auto l: data){
    std::string key = l.read_id + "-" +l.edge_id;
    if (std::find(filtered_links.begin(), filtered_links.end(), key) != filtered_links.end()){
      good_links.push_back(l);
    }
  }

  //devolver la lista de filtrados
  return good_links;
}


ReadPathVec PacbioPather::mapReads(){
  // get the reads and map them to the graph using the dictionary
  // returns a vector of paths
  std::cout << Date()<<": Executing getReadsLines..." << std::endl;
  auto links = getReadsLinks(true);
  std::cout << Date() << ": Done" << std::endl;

  // To store the paths
  ReadPath pb_paths_temp[seqVector->size()];
  ReadPathVec pb_paths;

  // for each read
  std::cout<<Date()<<": pathing "<<seqVector->size()<<" PacBio reads..."<<std::endl;
  std::atomic_uint_fast64_t pr(0),ppr(0);
  //#pragma omp parallel for
  for (std::uint32_t r=0; r<seqVector->size(); ++r){

    // filter the links for this reads
    auto read_links = links[r];

    // filter the shared roffsets
    auto offset_filter = readOffsetFilter(read_links);

    // filter the min match length TODO: this function

    // sort the vector
    std::sort(offset_filter.begin(), offset_filter.end(), linkreg_less_than());

    // Filter if the reverse complement (?)

    // Create vector of unique edge_ids
    std::vector<int> presentes;
    std::vector<linkReg> s_edges;
    for (auto a: offset_filter){
      if (std::find(presentes.begin(), presentes.end(), a.edge_id) == presentes.end()){
        s_edges.push_back(a);
        presentes.push_back(a.edge_id);
      }
    }



    std::vector<int> temp_path;

    if (s_edges.size() >= 2) {
      int poffset=s_edges[0].edge_offset;
      for (auto s: s_edges) {
        temp_path.push_back(s.edge_id);
      }
      pb_paths_temp[ppr++]=ReadPath(poffset, temp_path);
    }
    ++pr;
    //if (pr%500==0) std::cout<<Date()<<": "<<pr<<" reads processed, "<<ppr<<" pathed"<<std::endl;
  }
  std::cout<<Date()<<": "<<pr<<" reads processed, "<<ppr<<" pathed"<<std::endl;
  pb_paths.reserve(ppr);
  for (auto i=0;i<ppr;++i) pb_paths.push_back(pb_paths_temp[i]);

  return pb_paths;
}