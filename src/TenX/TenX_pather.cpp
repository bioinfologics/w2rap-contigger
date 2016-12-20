//
// Created by Gonzalo Garcia (TGAC) on 16/12/2016.
//

#include "TenX_pather.h"

TenXPather::TenXPather(std::vector<tenXRead>* aseqVector, HyperBasevector* ahbv): KMatch(31)  {
  //
  seqVector = aseqVector;
  hbv = ahbv;

//  makeIndexMap();
//  makeTagMap();
//  makeTagVector();
}

int TenXPather::createEmptyMap(HyperBasevector* hbv){
  // Create empty kmer map from the edges.
  auto edges = hbv->Edges();

  for (auto e = 0; e<edges.size(); ++e) {
    auto seq = edges[e].ToString();
    auto kv = ProduceKmers(seq);
    for (auto &k: kv)
      kmerTagMap[k.kmer]; // TODO: check if this is doing the correct thing
  }
  return 0;
}

int TenXPather::reads2kmerTagMap(){
  // for each read
  for (auto e = seqVector->begin(); e!=seqVector->end(); ++e){
    auto tag = e->tag.ToString();

    auto seq = e->r1.ToString();
    auto kv = ProduceKmers(seq);
    // for each kmer
    for (auto const& k: kv){
      if (kmerTagMap.find(k.kmer) != kmerTagMap.end()){
        kmerTagMap[k.kmer][tag]++;
      }
    }

    seq = e->r2.ToString();
    kv = ProduceKmers(seq);
    // for each kmer
    for (auto const& k: kv){
      if (kmerTagMap.find(k.kmer) != kmerTagMap.end()){
        kmerTagMap[k.kmer][tag]++;
      }
    }
  }
  return 0;
}

std::vector<int> TenXPather::readsTagQc(){
  // Calculate tags statistics
  //TODO: Should this do the same for indexes!? is that usefull!?
  std::vector<int> histogram (256, 0);
  for (const auto& k: tagMap){
    histogram[k.second.size()]++;
  }
  return histogram;
}

int TenXPather::makeTagVector(){
  // Create a vector of the tags to wrok with later (i can filter here as well)
  if (tagMap.size() == 0){
    // Is this corect or is better to create the dictionary with the funciton here is does not exist!?
    std::cout << "Empty tag map, create the map first to use this function"<< std::endl;
    return -1;
  }
  for (const auto& k: tagMap){
    tagVector.push_back(k.first);
  }
  return 0;
}

std::vector<tenXLink> TenXPather::processLinks(std::string read, int read_id){
  // TODO: avoid the copy returning the pointer to the object (LUIS)
//  std::cout << "Aligning read: " << read_id;
  auto kmer_matches = lookupRead(read);
  std::vector<tenXLink> temp_vector;
  if (kmer_matches.size()>0){
//    std::cout << kmer_matches.size()<<",";
    for (size_t a=0; a<kmer_matches.size(); ++a) {
      tenXLink temp_link;
      temp_link.read_id = read_id;
      temp_link.read_offset = kmer_matches[a].read_offset;
      temp_link.edge_id = kmer_matches[a].edge_id;
      temp_link.edge_offset = kmer_matches[a].edge_offset;
      temp_link.inv_edge_id = inv[kmer_matches[a].edge_id];
      temp_link.kmer = kmer_matches[a].kmer;
      temp_vector.push_back(temp_link);
      }
    }
//  std::cout << std::endl;
  return temp_vector;
}


std::vector<TenXPather::tagLink> TenXPather::getTagLinks(bool output_to_file=true){
  // This function takes the HBV, te tags and and the vector of tenx reads and aligns the tags to the
  // hbv using the reads.

  // Create the hbv map
  std::cout<< Date() << ": Size of the dictionary: " << edgeMap.size() << std::endl;
  std::cout<< Date() << ": loading edges and involution." << std::endl;

  auto edges = hbv->Edges();
  std::cout << "edges loaded: " << edges.size() << std::endl;
  std::cout << "edges loaded: " << inv.size() << std::endl;
//  vec<int> inv;
//  hbv->Involution(inv);

  std::vector<tagLink> tag_links;
  // Para cada tag en el vector de tag
  for (auto const &tag: tagVector){
//    std::cout << "Tag: " << tag;
    tagLink t_read_links;
    for (auto const &readid: tagMap[tag]){

      // Align R1
      std::string read1 = (*seqVector)[readid].r1.ToString();
      auto tv1 = processLinks(read1, readid);
      // Align R2
      std::string read2 = (*seqVector)[readid].r2.ToString();
      auto tv2 = processLinks(read2, readid);

//      std::cout << ","<< readid << " (" << tv1.size() << "," <<tv2.size()<< ") ";
      auto t_links_pair = std::make_pair(tv1, tv2);

      t_read_links.push_back(t_links_pair);
    }
    tag_links.push_back(t_read_links);

//    std::cout << std::endl;
  }
  return tag_links;
}

int TenXPather::makeIndexMap(bool to_disc=false){
  // Create a map with the index as key and a vector of indexes to the reads that have that index
  std::map<std::string, std::vector<int>> index_map;

  auto cont = 0;
  for (auto i=seqVector->begin(); i<seqVector->end(); ++i){
    auto idx = i->tag.ToString();
    index_map[idx].push_back(cont);
    cont++;
    }

  if (to_disc){
    for (auto const& i: index_map){
      std::cout << i.first << ":";
      for(auto const& v : i.second){
        std::cout << v << ",";
      }
      std::cout << std::endl;
    }
  }

  indexMap = index_map;
  return 0;
}

int TenXPather::makeTagMap(bool to_disc=false){
  // Create a map with the index as key and a vector of indexes to the reads that have that index
  std::map<std::string, std::vector<int>> tag_map;

  auto cont = 0;
  for (auto i=seqVector->begin(); i<seqVector->end(); ++i){
    auto idx = i->tag.ToString();
    tag_map[idx].push_back(cont);
    cont++;
  }

  if (to_disc){
    for (auto const& i: tag_map){
      std::cout << i.first << ":";
      for(auto const& v : i.second){
        std::cout << v << ",";
      }
      std::cout << std::endl;
    }
  }

  tagMap = tag_map;
  return 0;
}