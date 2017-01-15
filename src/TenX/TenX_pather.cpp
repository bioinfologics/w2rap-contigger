//
// Created by Gonzalo Garcia (TGAC) on 16/12/2016.
//

#include "TenX_pather.h"
#include "paths/PathFinder.h"

TenXPather::TenXPather(std::vector<tenXRead>* aseqVector, HyperBasevector* ahbv): KMatch(31)  {
  //
  seqVector = aseqVector;
  hbv = ahbv;

//  makeIndexMap();
//  makeTagMap();
//  makeTagVector();
}

TenXPather::tagktype TenXPather::kmerize_tag(std::string seq){

  int K=16;
  const tagktype KMER_MASK=( ((uint16_t)1)<<(K*2) )-1;

  const char *s = seq.c_str();
  tagktype fkmer = 0;
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
        return 0; // TODO: [GONZA] Check if this is correcttly done??
    }
  }
  return fkmer;
}


int TenXPather::createEmptyMap(HyperBasevector* hbv){
  // Create empty kmer map from the edges.
  int min_edge_length = 1500;
  auto & edges = hbv->Edges();
  auto kcount=0;

  for (auto &e:edges){
    if (e.size()>min_edge_length){
      kcount+=e.size()-31;//XXX: harcoded kmer size!!!!
    }
  }

  //store all kmers in an std::vector
  std::vector<std::pair<uint64_t, std::map<tagktype, int>>> all_kmers;
  all_kmers.reserve(kcount);
  std::cout << Date() << ": Number of kmers: " << kcount <<std::endl;


//  for (auto e = 0; e<edges.size(); ++e) {
#pragma omp parallel for
  for (auto e=0; e<edges.size(); ++e) { // TODO: [GONZA] fix this to run in a bigger machine, is like this for the map to fit in my laptop :/
    auto seq = edges[e].ToString();
    if (seq.length()>min_edge_length){
      auto kv = ProduceKmers(seq);
      for (auto k: kv) {
        std::map<tagktype, int> partest;
#pragma omp critical (vectorpush)
        all_kmers.push_back(std::make_pair(k.kmer, partest));
      }
    }
  }

//  std::sort(all_kmers.begin(),all_kmers.end(), kmer_pair_lessthan);
//  all_kmers.erase(std::unique(all_kmers.begin(), all_kmers.end()),all_kmers.end()); // This workds comparing keys!? do i need a comparator??
  kmerTagMap.insert(all_kmers.begin(), all_kmers.end());
  return 0;
}

int TenXPather::reads2kmerTagMap(){
  // for each read

#pragma omp parallel for
  for (auto e = 0; e < seqVector->size(); ++e){
    auto tag = kmerize_tag((*seqVector)[e].tag.ToString());
//    std::cout << "Tga: " << tag<<std::endl;

    auto seq = (*seqVector)[e].r1.ToString();
    auto kv = ProduceKmers(seq);
    // for each kmer
    for (auto const& k: kv){
      if (kmerTagMap.find(k.kmer) != kmerTagMap.end()){
#pragma omp critical (taginsert)
        kmerTagMap[k.kmer][tag]++;
      }
    }

    seq = (*seqVector)[e].r2.ToString();
    kv = ProduceKmers(seq);
    // for each kmer
    for (auto const& k: kv){
      if (kmerTagMap.find(k.kmer) != kmerTagMap.end()){
#pragma omp critical (taginsert)
        kmerTagMap[k.kmer][tag]++;
      }
    }
  }
  return 0;
}

std::vector<TenXPather::tagktype> TenXPather::getSequenceTags(std::string seq){
  // Given a sequence returns the tags asociated with those edge kmers in the kmerTagMap

  std::vector<TenXPather::tagktype> tags;
  auto kmers = ProduceKmers(seq);
  for (auto k: kmers){
    if (kmerTagMap.find(k.kmer) != kmerTagMap.end()){
      for (auto tag: kmerTagMap[k.kmer]){
        tags.push_back(tag.first);
      }
    }
  }
  return tags;
}

float TenXPather::edgeTagIntersection(std::string edgeFrom, std::string edgeTo, int roi) {
  // Given 2 edges as strings will return the set of tags that are common to both edges
  // Directional, edgeFrom (tail roi), edgeTo(head roi), takes the end of the first edge and the tail of the second

  std::vector<TenXPather::tagktype> intersection_tagset;

  // Set rois in both edges
  std::string edgeFrom_roi;
  if (edgeFrom.size()>roi){
    edgeFrom_roi = edgeFrom.substr(edgeFrom.size()-roi);
  } else {
    edgeFrom_roi = edgeFrom;
  }
  std::string edgeTo_roi;
  if (edgeTo.size()>roi){
    edgeTo_roi = edgeTo.substr(0, roi);
  } else {
    edgeTo_roi = edgeTo;
  }

  // Get tags from the map
  auto tagsFrom = getSequenceTags(edgeFrom_roi);
  auto tagsTo = getSequenceTags(edgeTo_roi);

  // [GONZA] TODO: check if this sort is really necesary with the set intersection!
  std::sort(tagsFrom.begin(), tagsFrom.end());
  std::sort(tagsTo.begin(), tagsTo.end());

  std::set_intersection(tagsFrom.begin(), tagsFrom.end(), tagsTo.begin(), tagsTo.end(), std::back_inserter(intersection_tagset));

  // Calculate the intersection score as the density of tags/kmer
  float intersection_score = (float)intersection_tagset.size() / (float)(edgeFrom_roi.size() + edgeTo_roi.size());
//  float intersection_score = (float)intersection_tagset.size();
  return intersection_score;
}

float TenXPather::kmerTagDensity(){
  // Get the tag kmer map and create a histogram of tag count, then choose a treshold
  std::vector<int> histogram (255, 0);
  for (auto kt: kmerTagMap){
    int count = 0;
    for (auto tagcount: kt.second){
      count += tagcount.second;
    }
    if (count > 255){
      count = 255;
    }
    histogram[count]++;
  }

  std::cout << "Tag count histogram begining" << std::endl;
  for (auto h=0; h<histogram.size(); ++h){
    std::cout << h << "," << histogram[h] << std::endl;
  }
  std::cout << "Histogram end" << std::endl;

  return 1.0;
}

//int TenXPather::resolve_regions(int large_frontier_size=500){
//  // get the regions and
//
//  // Create pathfinder instance
//  PathFinder pf(hbv, inv, ReadPathVec& paths, VecULongVec& invPaths, int min_reads = 5 );
//
//  for (int e = 0; e < hbv->EdgeObjectCount(); ++e) {
//    if (e < inv[e] && hbv->EdgeObject(e).size() < large_frontier_size) {
//      auto f=pf.get_all_long_frontiers(e, large_frontier_size);
//      for (auto x: f){
//        std::cout << f << std::endl;
//      }
//    }
//  }
//}

//std::vector<int> TenXPather::readsTagQc(){
//  // Calculate tags statistics
//  //TODO: Should this do the same for indexes!? is that usefull!?
//  std::vector<int> histogram (256, 0);
//  for (const auto& k: tagMap){
//    histogram[k.second.size()]++;
//  }
//  return histogram;
//}

//int TenXPather::makeTagVector(){
//  // Create a vector of the tags to wrok with later (i can filter here as well)
//  if (tagMap.size() == 0){
//    // Is this corect or is better to create the dictionary with the funciton here is does not exist!?
//    std::cout << "Empty tag map, create the map first to use this function"<< std::endl;
//    return -1;
//  }
//  for (const auto& k: tagMap){
//    tagVector.push_back(k.first);
//  }
//  return 0;
//}

//std::vector<tenXLink> TenXPather::processLinks(std::string read, int read_id){
//  // TODO: avoid the copy returning the pointer to the object (LUIS)
////  std::cout << "Aligning read: " << read_id;
//  auto kmer_matches = lookupRead(read);
//  std::vector<tenXLink> temp_vector;
//  if (kmer_matches.size()>0){
////    std::cout << kmer_matches.size()<<",";
//    for (size_t a=0; a<kmer_matches.size(); ++a) {
//      tenXLink temp_link;
//      temp_link.read_id = read_id;
//      temp_link.read_offset = kmer_matches[a].read_offset;
//      temp_link.edge_id = kmer_matches[a].edge_id;
//      temp_link.edge_offset = kmer_matches[a].edge_offset;
//      temp_link.inv_edge_id = inv[kmer_matches[a].edge_id];
//      temp_link.kmer = kmer_matches[a].kmer;
//      temp_vector.push_back(temp_link);
//      }
//    }
////  std::cout << std::endl;
//  return temp_vector;
//}


//std::vector<TenXPather::tagLink> TenXPather::getTagLinks(bool output_to_file=true){
//  // This function takes the HBV, te tags and and the vector of tenx reads and aligns the tags to the
//  // hbv using the reads.
//
//  // Create the hbv map
//  std::cout<< Date() << ": Size of the dictionary: " << edgeMap.size() << std::endl;
//  std::cout<< Date() << ": loading edges and involution." << std::endl;
//
//  auto edges = hbv->Edges();
//  std::cout << "edges loaded: " << edges.size() << std::endl;
//  std::cout << "edges loaded: " << inv.size() << std::endl;
////  vec<int> inv;
////  hbv->Involution(inv);
//
//  std::vector<tagLink> tag_links;
//  // Para cada tag en el vector de tag
//  for (auto const &tag: tagVector){
////    std::cout << "Tag: " << tag;
//    tagLink t_read_links;
//    for (auto const &readid: tagMap[tag]){
//
//      // Align R1
//      std::string read1 = (*seqVector)[readid].r1.ToString();
//      auto tv1 = processLinks(read1, readid);
//      // Align R2
//      std::string read2 = (*seqVector)[readid].r2.ToString();
//      auto tv2 = processLinks(read2, readid);
//
////      std::cout << ","<< readid << " (" << tv1.size() << "," <<tv2.size()<< ") ";
//      auto t_links_pair = std::make_pair(tv1, tv2);
//
//      t_read_links.push_back(t_links_pair);
//    }
//    tag_links.push_back(t_read_links);
//
////    std::cout << std::endl;
//  }
//  return tag_links;
//}
//
//int TenXPather::makeIndexMap(bool to_disc=false){
//  // Create a map with the index as key and a vector of indexes to the reads that have that index
//  std::map<std::string, std::vector<int>> index_map;
//
//  auto cont = 0;
//  for (auto i=seqVector->begin(); i<seqVector->end(); ++i){
//    auto idx = i->tag.ToString();
//    index_map[idx].push_back(cont);
//    cont++;
//    }
//
//  if (to_disc){
//    for (auto const& i: index_map){
//      std::cout << i.first << ":";
//      for(auto const& v : i.second){
//        std::cout << v << ",";
//      }
//      std::cout << std::endl;
//    }
//  }
//
//  indexMap = index_map;
//  return 0;
//}
//
//int TenXPather::makeTagMap(bool to_disc=false){
//  // Create a map with the index as key and a vector of indexes to the reads that have that index
//  std::map<std::string, std::vector<int>> tag_map;
//
//  auto cont = 0;
//  for (auto i=seqVector->begin(); i<seqVector->end(); ++i){
//    auto idx = i->tag.ToString();
//    tag_map[idx].push_back(cont);
//    cont++;
//  }
//
//  if (to_disc){
//    for (auto const& i: tag_map){
//      std::cout << i.first << ":";
//      for(auto const& v : i.second){
//        std::cout << v << ",";
//      }
//      std::cout << std::endl;
//    }
//  }
//
//  tagMap = tag_map;
//  return 0;
//}