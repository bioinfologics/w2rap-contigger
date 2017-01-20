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
  for (auto e=0; e<10000; ++e) { // TODO: [GONZA] fix this to run in a bigger machine, is like this for the map to fit in my laptop :/
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

  kmerTagMap.insert(all_kmers.begin(), all_kmers.end());
  return 0;
}

int TenXPather::reads2kmerTagMap(){
  // Load the reads into the map
  std::cout << Date() << ": Filling the map with reads" << std::endl;
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
  std::cout << Date() << ": Done filling the maps" << std::endl;
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
  float intersection_score = (float)intersection_tagset.size()*2 / (float)(edgeFrom_roi.size() + edgeTo_roi.size());
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


//LocalPaths_TX::LocalPaths_TX(HyperBasevector* hbv, std::vector<std::vector<uint64_t>> pair_solutions, vec<int>& to_right, TenXPather* txp, std::vector<BaseVec>& edges)
//    : LocalPaths(hbv, pair_solutions, to_right, edges){
//  mTxp = txp;
//}

std::vector<uint64_t> LocalPaths_TX::choose_best_path(std::vector<std::vector<uint64_t>> *alternative_paths) {
  // Choose the best path for the combination from the pairs list

  // If there is only one posible path return that path and finish
  if (alternative_paths->size() == 0){
    std::cout <<"This path still CEROOOO" <<std::endl;
    std::vector<uint64_t> nopaths;
    return nopaths;
  }
  if (alternative_paths->size() == 1){
    std::cout << "Only one patha available: " << std::endl;
    return (*alternative_paths)[0];
  } else {
    // If there is more than one path vote for the best (the criteria here is most tag density (presentTags/totalKmers)
    float best_path = 0.0;
    float best_path_score = 0.0;
    for (auto path_index = 0; path_index < alternative_paths->size(); ++path_index) {
      float cpath_score = 0;
      for (auto ei = 0; ei < (*alternative_paths)[path_index].size() - 1; ++ei) {
        auto from_edge_string = (*mEdges)[(*alternative_paths)[path_index][ei]].ToString();
        auto to_edge_string = (*mEdges)[(*alternative_paths)[path_index][ei + 1]].ToString();
        cpath_score += mTxp->edgeTagIntersection(from_edge_string, to_edge_string, 1500);
      }
      if (cpath_score > best_path_score) {
        best_path = path_index;
        best_path_score = cpath_score;
      }
    }
    std::cout << "Best path selected: " << best_path << ", score: " << best_path_score << std::endl;
    for (auto p=0; p<(*alternative_paths)[best_path].size(); ++p){
      std::cout << (*alternative_paths)[best_path][p] << ",";
    }
    std::cout << std::endl;

    return (*alternative_paths)[best_path];
  }
}