//
// Created by Gonzalo Garcia (TGAC) on 16/12/2016.
//

#include <kmers/kmatch/KMatch.h>
#include "TenX_pather.h"
#include "paths/PathFinder.h"


TenXPather::TenXPather(std::vector<tenXRead>& aseqVector, HyperBasevector& ahbv, vec<int>& ainv, int min_reads, std::vector<BaseVec>& edges, ReadPathVec& apaths, VecULongVec& ainvPaths)
    : KMatch(31), PathFinder(ahbv, ainv, apaths, ainvPaths, 5), seqVector(aseqVector), mEdges(edges) {}

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
  std::cout << Date() << ": Create empty map Begins" << std::endl;
  int min_edge_length = 1500;
  auto & edges = hbv->Edges();
  auto kcount=0;

  for (auto &e:edges){
    if (e.size()>min_edge_length){
      kcount+=e.size()-31;//XXX: harcoded kmer size!!!!
    }
  }

  //store all kmers in an std::vector
  std::vector<std::pair<uint64_t, tbb::concurrent_unordered_multiset<tagktype>>> all_kmers;
  all_kmers.reserve(kcount);
  std::cout << Date() << ": Number of kmers: " << kcount <<std::endl;


  std::cout  << Date() << ": Creating the kmer:set pair to insert" << std::endl;
//  for (auto e = 0; e<edges.size(); ++e) {
#pragma omp parallel for
  for (auto e=0; e<10000; ++e) { // TODO: [GONZA] fix this to run in a bigger machine, is like this for the map to fit in my laptop :/
    auto seq = edges[e].ToString();
    if (seq.length()>min_edge_length){
      auto kv = ProduceKmers(seq);
      for (auto k: kv) {
        tbb::concurrent_unordered_multiset<tagktype> empty_multiset;
#pragma omp critical (vectorpush)
        all_kmers.push_back(std::make_pair(k.kmer, empty_multiset));
      }
    }
  }

  std::cout  << Date() << ": Inserting pairs in the map" << std::endl;
  kmerTagMap.insert(all_kmers.begin(), all_kmers.end());
  std::cout  << Date() << ": Done, number of keys in the vector: "<< all_kmers.size() << ", Keys in th map: " << kmerTagMap.size() << std::endl;
  return 0;
}

int TenXPather::reads2kmerTagMap(){
  // Load the reads into the map
  std::cout << Date() << ": Filling the map with reads" << std::endl;
#pragma omp parallel for
  for (auto e = 0; e < seqVector.size(); ++e){
    auto tag = kmerize_tag(seqVector[e].tag.ToString());
//    std::cout << "Tag: " << tag<<std::endl;

    auto seq = seqVector[e].r1.ToString();
    insert_kmertags_in_edgemap(tag, seq);

    seq = seqVector[e].r2.ToString();
    insert_kmertags_in_edgemap(tag, seq);

  }
  std::cout << Date() << ": Done filling the maps" << std::endl;
  return 0;
}

void TenXPather::insert_kmertags_in_edgemap(const TenXPather::tagktype tag, const String &seq) {
  auto kv = ProduceKmers(seq);
  // for each kmer
  for (auto const& k: kv){
      const auto klookup(kmerTagMap.find(k.kmer));

      // Reduce the number of lookups by storing the found iterator and inserting there rather than looking it up again with "at"
      if (klookup != kmerTagMap.cend()){
        klookup->second.insert(tag);
      }
    }
}

std::vector<TenXPather::tagktype> TenXPather::getSequenceTags(const std::string seq){
  // Given a sequence returns the tags associated with those edge kmers in the kmerTagMap

  std::vector<TenXPather::tagktype> tags;
  const auto kmers = ProduceKmers(seq);
  for (const auto &k: kmers){
    const auto klookup(kmerTagMap.find(k.kmer));
    if (klookup != kmerTagMap.cend()){
      for (const auto tag: klookup->second){
        tags.push_back(tag);
      }
    }
  }
  return tags;
}

float TenXPather::edgeTagIntersection(const string edgeFrom, const string edgeTo, const int roi) {
  // Given 2 edges as strings will return the set of tags that are common to both edges
  // Directional, edgeFrom (tail roi), edgeTo(head roi), takes the end of the first edge and the tail of the second

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

  // [GONZA] TODO: check if this sort is really necessary with the set intersection!
  std::sort(tagsFrom.begin(), tagsFrom.end());
  std::sort(tagsTo.begin(), tagsTo.end());

  std::vector<TenXPather::tagktype> intersection_tagset;
  std::set_intersection(tagsFrom.begin(), tagsFrom.end(), tagsTo.begin(), tagsTo.end(), std::back_inserter(intersection_tagset));

  // Calculate the intersection score as the density of tags/kmer
  float intersection_score = (float)intersection_tagset.size()*2 / (float)(edgeFrom_roi.size() + edgeTo_roi.size());
  return intersection_score;
}

//TODO: Esto no anda...
float TenXPather::kmerTagDensity(){
  // Get the tag kmer map and create a histogram of tag count, then choose a threshold
  std::vector<int> histogram (255, 0);
  for (auto kt: kmerTagMap){
    int count (kt.second.size());
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

std::vector<uint64_t> TenXPather::choose_best_path(std::vector<std::vector<uint64_t>> *alternative_paths) {
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
        auto from_edge_string = mEdges[(*alternative_paths)[path_index][ei]].ToString();
        auto to_edge_string = mEdges[(*alternative_paths)[path_index][ei + 1]].ToString();
        cpath_score += edgeTagIntersection(from_edge_string, to_edge_string, 1500);
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


void TenXPather::solve_region_using_TenX(const uint64_t large_frontier_size, const bool verbose_separation,
                                         const float score_threshold) {
  //find a complex path
  // [GOnza] TODO: separate this function into 2 differets, one to create the map an the other to test the permutations
  const auto edges = mHBV.Edges();

  uint64_t qsf=0,qsf_paths=0;
  uint64_t msf=0,msf_paths=0;
  init_prev_next_vectors();
  std::cout<<"vectors initialised"<<std::endl;
  std::set<std::array<std::vector<uint64_t>,2>> seen_frontiers,solved_frontiers;
  std::vector<std::vector<uint64_t>> paths_to_separate;
  int solved_regions = 0;
  int unsolved_regions = 0;
  //TODO: EdgeObjectCount can overflow
  for (int e = 0; e < mHBV.EdgeObjectCount(); ++e) {
    if (e < mInv[e] && mHBV.EdgeObject(e).size() < large_frontier_size) {
      // Get the frontiers fo the edge
      // [GONZA] TODO: check the return details to document
      auto f=get_all_long_frontiers(e, large_frontier_size);
      if (f[0].size()>1 and f[1].size()>1 and f[0].size() == f[1].size() and seen_frontiers.count(f)==0){
        seen_frontiers.insert(f);
        bool single_dir=true;

        std::unordered_map<std::string, float> shared_paths;
        for (auto in_e:f[0]) for (auto out_e:f[1]) if (in_e==out_e) {single_dir=false;break;}

        if (single_dir) {
          // If there is a region to resolve within the frontiers
          std::cout<<" Single direction frontiers for complex region on edge "<<e<<" IN:"<<path_str(f[0])<<" OUT: "<<path_str(f[1])<<std::endl;
          std::vector<int> in_used(f[0].size(),0);
          std::vector<int> out_used(f[1].size(),0);
          std::vector<std::vector<uint64_t>> first_full_paths;
          bool reversed=false;

          // intersect all pairs of ins and outs to score, save the scores in a map to score the combinations later
          for (auto in_i=0;in_i<f[0].size();++in_i) {
            auto in_e=f[0][in_i];
            for (auto out_i=0;out_i<f[1].size();++out_i) {
              auto out_e=f[1][out_i];

              int edges_in_path;
              const auto in_e_seq = edges[in_e].ToString();
              const auto out_e_seq = edges[out_e].ToString();

              // Intersect the tags for the edges
              auto intersection_score = edgeTagIntersection(in_e_seq, out_e_seq, 1500);
              if (intersection_score>score_threshold){
                // if the edges overlap in the tagspace thay are added to the map and the combination is markes in the used edges
                const std::string pid(std::to_string(in_e) + "-" + std::to_string(out_e));
                shared_paths[pid] += intersection_score; // This should score the link based in the number of tags that tha pair shares
//                std::cout << " from-to edges: " << pid << "= " << intersection_score << std::endl;
                out_used[out_i]++;
                in_used[in_i]++;
              }
            }
          } // When this is done i get a map of all combinations of ins and outs to score the permutations in the next step

          // Here all combinations are counted, now i need to get the best configuration between nodes
          auto in_frontiers = f[0];
          auto out_frontiers = f[1];
          float max_score = -9999.0;
          int perm_number = 0;

          std::vector<int> max_score_permutation;
          do {
            // Vectors to count seen edges (to check that all edges are included in th solution)
            std::vector<int> seen_in(in_frontiers.size(), 0);
            std::vector<int> seen_out(in_frontiers.size(), 0);

            //std::cout << "-----------------------------Testing permutaiton ------------------" << std::endl;
            float current_score = 0;
            for (auto pi=0; pi<in_frontiers.size(); ++pi){
              std::string index = std::to_string(in_frontiers[pi])+"-"+std::to_string(out_frontiers[pi]);
              if ( shared_paths.find(index) != shared_paths.end() ){
                // Mark the pair as seen in this iteration
                seen_in[pi]++;
                seen_out[pi]++;
                current_score += shared_paths[index];
              }
            }
            // Check that all boundaries are used in the permutation
            bool all_used = true;
            for (auto a=0; a<seen_in.size(); ++a){
              if (0==seen_in[a] or 0==seen_out[a]){
                all_used = false;
              }
            }

            if (current_score>max_score and all_used){
              max_score = current_score;
              max_score_permutation.clear();
              for (auto aix: out_frontiers){
                max_score_permutation.push_back(aix);
              }
            }
            perm_number++;
          } while (std::next_permutation(out_frontiers.begin(), out_frontiers.end()));

          // Get the solution
          if (max_score>score_threshold){
            std::cout << " Found solution to region: " <<std::endl;
            solved_regions++;

            std::vector<std::vector<uint64_t>> wining_permutation;
            for (auto ri=0; ri<max_score_permutation.size(); ++ri){
              std::cout << in_frontiers[ri] << "(" << mInv[in_frontiers[ri]] << ") --> " << max_score_permutation[ri] << "("<< mInv[max_score_permutation[ri]] <<"), Score: "<<  max_score << std::endl;
              std::vector<uint64_t> tp = {in_frontiers[ri], max_score_permutation[ri]};
              wining_permutation.push_back(tp);
            }
            std::cout << "--------------------" << std::endl;

            // Get all paths between the solved pair of edges and get the best one for each pair
            auto all_paths_found = find_all_solution_paths(wining_permutation);
            std::cout << "All paths done" << std::endl;
            for (auto spi = 0; spi<all_paths_found.size(); ++spi){
              auto sv = all_paths_found[spi];
              paths_to_separate.push_back(sv);
            }
            std::cout << "--------------------" << std::endl;
          } else {
//                        std::cout << "Region not resolved, not enough links or bad combinations" << max_score << std::endl;
            unsolved_regions++;
          }

        }
      }
    }
  }
  std::cout << "========================" << std::endl;
  std::cout << "Solved regions: " << solved_regions << ", Unsolved regions: " << unsolved_regions << " --> " << (float)solved_regions / (float)(solved_regions + unsolved_regions) * 100 << " %" << std::endl;
  std::cout << "========================" << std::endl;

  std::cout << "List of paths to separate" << std::endl;
  for (auto ss: paths_to_separate){
    for (auto sx: ss){
      std::cout << sx << ",";
    }
    std::cout << std::endl;
  }
  std::cout << "Paths to separate: " << paths_to_separate.size() << std::endl;

  uint64_t sep=0;
  std::map<uint64_t,std::vector<uint64_t>> old_edges_to_new;
  for (auto p: paths_to_separate){
    if (old_edges_to_new.count(p.front()) > 0 or old_edges_to_new.count(p.back()) > 0) {
      std::cout<<"WARNING: path starts or ends in an already modified edge, skipping"<<std::endl;
      continue;
    }

    auto oen=separate_path(p, verbose_separation);
    std::cout << "End separate paths, eon size: " << oen.size() << std::endl;
    if (oen.size()>0) {
      for (auto et:oen){
        if (old_edges_to_new.count(et.first)==0) old_edges_to_new[et.first]={};
        for (auto ne:et.second) old_edges_to_new[et.first].push_back(ne);
      }
      sep++;
      std::cout << "separated_paths: " << sep << std::endl;
    }
  }
//    if (old_edges_to_new.size()>0) {
//        migrate_readpaths(old_edges_to_new);
//    }
  std::cout<<" "<<sep<<" paths separated!"<<std::endl;
}