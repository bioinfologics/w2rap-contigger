//
// Created by Bernardo Clavijo (TGAC) on 11/07/2016.
//

#include "TenX/PathFinder_tx.h"

PathFinder_tx::PathFinder_tx (TenXPather* txp, HyperBasevector* hbv, vec<int> inv, int min_reads = 5) {
    mHBV = hbv;
    mInv = inv;
    mMinReads = min_reads;
    mTxp = txp;

    hbv->ToLeft(mToLeft);
    hbv->ToRight(mToRight);
}

void PathFinder_tx::init_prev_next_vectors(){
    //TODO: this is stupid duplication of the digraph class, but it's so weird!!!
    prev_edges.resize(mToLeft.size());
    next_edges.resize(mToRight.size());
    for (auto e=0;e<mToLeft.size();++e){

        uint64_t prev_node=mToLeft[e];

        prev_edges[e].resize(mHBV->ToSize(prev_node));
        for (int i=0;i<mHBV->ToSize(prev_node);++i){
            prev_edges[e][i]=mHBV->EdgeObjectIndexByIndexTo(prev_node,i);
        }

        uint64_t next_node=mToRight[e];

        next_edges[e].resize(mHBV->FromSize(next_node));
        for (int i=0;i<mHBV->FromSize(next_node);++i){
            next_edges[e][i]=mHBV->EdgeObjectIndexByIndexFrom(next_node,i);
        }
    }
}


std::string PathFinder_tx::path_str(std::vector<uint64_t> path) {
    std::string s="[";
    for (auto p:path){
        s+=std::to_string(p)+":"+std::to_string(mInv[p])+" ";//+" ("+std::to_string(mHBV->EdgeObject(p).size())+"bp "+std::to_string(paths_per_kbp(p))+"ppk)  ";
    }
    s+="]";
    return s;
}


void PathFinder_tx::untangle_complex_in_out_choices(uint64_t large_frontier_size, bool verbose_separation, float score_threshold) {
    //find a complex path
    // [GOnza] TODO: separate this function into 2 differets, one to create the map an the other to test the permutations
    auto edges = mHBV->Edges();

    uint64_t qsf=0,qsf_paths=0;
    uint64_t msf=0,msf_paths=0;
    init_prev_next_vectors();
    std::cout<<"vectors initialised"<<std::endl;
    std::set<std::array<std::vector<uint64_t>,2>> seen_frontiers,solved_frontiers;
    std::vector<std::vector<uint64_t>> paths_to_separate;
    int solved_regions = 0;
    int unsolved_regions = 0;
    for (int e = 0; e < mHBV->EdgeObjectCount(); ++e) {
        if (e < mInv[e] && mHBV->EdgeObject(e).size() < large_frontier_size) {
            // Get the frontiers fo the edge
            // [GONZA] TODO: check the return details to document
            auto f=get_all_long_frontiers(e, large_frontier_size);
            if (f[0].size()>1 and f[1].size()>1 and f[0].size() == f[1].size() and seen_frontiers.count(f)==0){
                seen_frontiers.insert(f);
                bool single_dir=true;
                std::map<std::string, float> shared_paths;
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
                            std::string pid;
                            auto in_e_seq = edges[in_e].ToString();
                            auto out_e_seq = edges[out_e].ToString();

                            // Intersect the tags for the edges
                            auto intersection_score = mTxp->edgeTagIntersection(in_e_seq, out_e_seq, 1500);
                            if (intersection_score>score_threshold){
                                // if the edges overlap in the tagspace thay are added to the map and the combination is markes in the used edges
                                pid = std::to_string(in_e) + "-" + std::to_string(out_e);
                                shared_paths[pid] += intersection_score; // This should score the link based in the number of tags that tha pair shares
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
                        // Vectors to count seen edges (check that all edges are included in th solution)
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

                        // Fill intermediate nodes
                        LocalPaths_TX lp (mHBV, wining_permutation, mToRight, mTxp, edges);
                        lp.find_all_solution_paths();
                        std::cout << "All paths done" << std::endl;
                        for (auto spi = 0; spi<lp.all_paths.size(); ++spi){
                            auto sv = lp.all_paths[spi];
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

//    std::cout<<"Complex Regions solved by paths: "<<solved_frontiers.size() <<"/"<<seen_frontiers.size()<<" comprising "<<paths_to_separate.size()<<" paths to separate"<< std::endl;
//    //std::cout<<"Complex Regions quasi-solved by paths (not acted on): "<< qsf <<"/"<<seen_frontiers.size()<<" comprising "<<qsf_paths<<" paths to separate"<< std::endl;
//    //std::cout<<"Multiple Solution Regions (not acted on): "<< msf <<"/"<<seen_frontiers.size()<<" comprising "<<msf_paths<<" paths to separate"<< std::endl;
    std::cout << "Paths to separate" << std::endl;

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

std::array<std::vector<uint64_t>,2> PathFinder_tx::get_all_long_frontiers(uint64_t e, uint64_t large_frontier_size){
    //TODO: return all components in the region
    std::set<uint64_t> seen_edges, to_explore={e}, in_frontiers, out_frontiers;

    while (to_explore.size()>0){
        std::set<uint64_t> next_to_explore;
        for (auto x:to_explore){ //to_explore: replace rather and "update" (use next_to_explore)

            if (!seen_edges.count(x)){

                //What about reverse complements and paths that include loops that "reverse the flow"?
                if (seen_edges.count(mInv[x])) return std::array<std::vector<uint64_t>,2>(); //just cancel for now

                for (auto p:prev_edges[x]) {
                    if (mHBV->EdgeObject(p).size() >= large_frontier_size )  {
                        //What about frontiers on both sides?
                        in_frontiers.insert(p);
                        for (auto other_n:next_edges[p]){
                            if (!seen_edges.count(other_n)) {
                                if (mHBV->EdgeObject(other_n).size() >= large_frontier_size) {
                                    out_frontiers.insert(other_n);
                                    seen_edges.insert(other_n);
                                }
                                else next_to_explore.insert(other_n);
                            }
                        }
                    }
                    else if (!seen_edges.count(p)) next_to_explore.insert(p);
                }

                for (auto n:next_edges[x]) {
                    if (mHBV->EdgeObject(n).size() >= large_frontier_size) {
                        //What about frontiers on both sides?
                        out_frontiers.insert(n);
                        for (auto other_p:prev_edges[n]){
                            if (!seen_edges.count(other_p)) {
                                if (mHBV->EdgeObject(other_p).size() >= large_frontier_size) {
                                    in_frontiers.insert(other_p);
                                    seen_edges.insert(other_p);
                                }
                                else next_to_explore.insert(other_p);
                            }
                        }
                    }
                    else if (!seen_edges.count(n)) next_to_explore.insert(n);
                }
                seen_edges.insert(x);
            }
            if (seen_edges.size()>50) {
                return std::array<std::vector<uint64_t>,2>();
            }

        }
        to_explore=next_to_explore;
    }

    if (to_explore.size()>0) return std::array<std::vector<uint64_t>,2>();
    //the "canonical" representation is the one that has the smalled edge on the first vector, and bot ordered

    if (in_frontiers.size()>0 and out_frontiers.size()>0) {
        uint64_t min_in=*in_frontiers.begin();
        for (auto i:in_frontiers){
            if (i<min_in) min_in=i;
            if (mInv[i]<min_in) min_in=mInv[i];
        }
        uint64_t min_out=*out_frontiers.begin();
        for (auto i:out_frontiers){
            if (i<min_out) min_out=i;
            if (mInv[i]<min_out) min_out=mInv[i];
        }
        if (min_out<min_in){
            std::set<uint64_t> new_in, new_out;
            for (auto e:in_frontiers) new_out.insert(mInv[e]);
            for (auto e:out_frontiers) new_in.insert(mInv[e]);
            in_frontiers=new_in;
            out_frontiers=new_out;
        }
    }

    //std::sort(in_frontiers.begin(),in_frontiers.end());
    //std::sort(out_frontiers.begin(),out_frontiers.end());
    std::array<std::vector<uint64_t>,2> frontiers={std::vector<uint64_t>(in_frontiers.begin(),in_frontiers.end()),std::vector<uint64_t>(out_frontiers.begin(),out_frontiers.end())};



    return frontiers;
}


std::map<uint64_t,std::vector<uint64_t>> PathFinder_tx::separate_path(std::vector<uint64_t> p, bool verbose_separation){

    //TODO XXX: proposed version 1 (never implemented)
    //Creates new edges for the "repeaty" parts of the path (either those shared with other edges or those appearing multiple times in this path).
    //moves paths across to the new reapeat instances as needed
    //changes neighbourhood (i.e. creates new vertices and moves the to and from for the implicated edges).

    //creates a copy of each node but the first and the last, connects only linearly to the previous copy,
    //std::cout<<std::endl<<"Separating path"<<std::endl;
    std::set<uint64_t> edges_fw;
    std::set<uint64_t> edges_rev;
    for (auto e:p){//TODO: this is far too astringent...
        edges_fw.insert(e);
        edges_rev.insert(mInv[e]);

        if (edges_fw.count(mInv[e]) ||edges_rev.count(e) ){ //std::cout<<"PALINDROME edge detected, aborting!!!!"<<std::endl;
            return {};}
    }
    //create two new vertices (for the FW and BW path)
    uint64_t current_vertex_fw=mHBV->N(),current_vertex_rev=mHBV->N()+1;
    mHBV->AddVertices(2);
    //migrate connections (dangerous!!!)
    if (verbose_separation) std::cout<<"Migrating edge "<<p[0]<<" To node old: "<<mToRight[p[0]]<<" new: "<<current_vertex_fw<<std::endl;
    mHBV->GiveEdgeNewToVx(p[0],mToRight[p[0]],current_vertex_fw);
    if (verbose_separation) std::cout<<"Migrating edge "<<mInv[p[0]]<<" From node old: "<<mToLeft[mInv[p[0]]]<<" new: "<<current_vertex_rev<<std::endl;
    mHBV->GiveEdgeNewFromVx(mInv[p[0]],mToLeft[mInv[p[0]]],current_vertex_rev);
    std::map<uint64_t,std::vector<uint64_t>> old_edges_to_new;

    for (auto ei=1;ei<p.size()-1;++ei){
        //add a new vertex for each of FW and BW paths
        uint64_t prev_vertex_fw=current_vertex_fw,prev_vertex_rev=current_vertex_rev;
        //create two new vertices (for the FW and BW path)
        current_vertex_fw=mHBV->N();
        current_vertex_rev=mHBV->N()+1;
        mHBV->AddVertices(2);

        //now, duplicate next edge for the FW and reverse path
        auto nef=mHBV->AddEdge(prev_vertex_fw,current_vertex_fw,mHBV->EdgeObject(p[ei]));
        if (verbose_separation)  std::cout<<"Edge "<<nef<<": copy of "<<p[ei]<<": "<<prev_vertex_fw<<" - "<<current_vertex_fw<<std::endl;
        mToLeft.push_back(prev_vertex_fw);
        mToRight.push_back(current_vertex_fw);
        if (! old_edges_to_new.count(p[ei]))  old_edges_to_new[p[ei]]={};
        old_edges_to_new[p[ei]].push_back(nef);

        auto ner=mHBV->AddEdge(current_vertex_rev,prev_vertex_rev,mHBV->EdgeObject(mInv[p[ei]]));
        if (verbose_separation) std::cout<<"Edge "<<ner<<": copy of "<<mInv[p[ei]]<<": "<<current_vertex_rev<<" - "<<prev_vertex_rev<<std::endl;
        mToLeft.push_back(current_vertex_rev);
        mToRight.push_back(prev_vertex_rev);
        if (! old_edges_to_new.count(mInv[p[ei]]))  old_edges_to_new[mInv[p[ei]]]={};
        old_edges_to_new[mInv[p[ei]]].push_back(ner);
        std::cout << "before pushing" << std::endl;
        mInv.push_back(ner);
        mInv.push_back(nef);
//        mEdgeToPathIds.resize(mEdgeToPathIds.size()+2); // [GONZA] TODO: this fals for some reason commented now

    }
    if (verbose_separation) std::cout<<"Migrating edge "<<p[p.size()-1]<<" From node old: "<<mToLeft[p[p.size()-1]]<<" new: "<<current_vertex_fw<<std::endl;
    mHBV->GiveEdgeNewFromVx(p[p.size()-1],mToLeft[p[p.size()-1]],current_vertex_fw);
    if (verbose_separation) std::cout<<"Migrating edge "<<mInv[p[p.size()-1]]<<" To node old: "<<mToRight[mInv[p[p.size()-1]]]<<" new: "<<current_vertex_rev<<std::endl;
    mHBV->GiveEdgeNewToVx(mInv[p[p.size()-1]],mToRight[mInv[p[p.size()-1]]],current_vertex_rev);

    //TODO: cleanup new isolated elements and leading-nowhere paths.
    //for (auto ei=1;ei<p.size()-1;++ei) mHBV->DeleteEdges({p[ei]});
    return old_edges_to_new;

}

//TODO: this should probably be called just once
void PathFinder_tx::migrate_readpaths(std::map<uint64_t,std::vector<uint64_t>> edgemap){
    //Migrate readpaths: this changes the readpaths from old edges to new edges
    //if an old edge has more than one new edge it tries all combinations until it gets the paths to map
    //if more than one combination is valid, this chooses at random among them (could be done better? should the path be duplicated?)
    mHBV->ToLeft(mToLeft);
    mHBV->ToRight(mToRight);
    for (auto &p:mPaths){
        std::vector<std::vector<uint64_t>> possible_new_edges;
        bool translated=false,ambiguous=false;
        for (auto i=0;i<p.size();++i){
            if (edgemap.count(p[i])) {
                possible_new_edges.push_back(edgemap[p[i]]);
                if (not translated) translated=true;
                if (possible_new_edges.back().size()>1) ambiguous=true;
            }
            else possible_new_edges.push_back({p[i]});
        }
        if (translated){
            if (not ambiguous){ //just straigh forward translation
                for (auto i=0;i<p.size();++i) p[i]=possible_new_edges[i][0];
            }
            else {
                //ok, this is the complicated case, we first generate all possible combinations
                std::vector<std::vector<uint64_t>> possible_paths={{}};
                for (auto i=0;i<p.size();++i) {//for each position
                    std::vector<std::vector<uint64_t>> new_possible_paths;
                    for (auto pp:possible_paths) { //take every possible one
                        for (auto e:possible_new_edges[i]) {
                            //if i>0 check there is a real connection to the previous edge
                          if (i == 0 and pp[i] == e){
                            std::cout << "first edge of path does not need migrating" << std::endl;
                          }
                            if (i == 0 or (mToRight[pp.back()]==mToLeft[e])) {
                                new_possible_paths.push_back(pp);
                                new_possible_paths.back().push_back(e);
                            }
                        }
                    }
                    possible_paths=new_possible_paths;
                    if (possible_paths.size()==0) break;
                }
                if (possible_paths.size()==0){
                    std::cout<<"Warning, a path could not be updated, truncating it to its first element!!!!"<<std::endl;
                    p.resize(1);
                }
                else{
                    std::srand (std::time(NULL));
                    //randomly choose a path
                    int r=std::rand()%possible_paths.size();
                    for (auto i=0;i<p.size();++i) p[i]=possible_paths[r][i];
                }
            }
        }
    }

}
