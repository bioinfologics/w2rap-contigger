//
// Created by Bernardo Clavijo (TGAC) on 11/07/2016.
//

#include "PathFinder.h"


std::string PathFinder::edge_pstr(uint64_t e){
    return "e"+std::to_string(e)+"("+std::to_string(mHBV.EdgeObject(e).size())+"bp "+std::to_string(paths_per_kbp(e))+"ppk)";
};

void PathFinder::init_prev_next_vectors(){
    //TODO: this is stupid duplication of the digraph class, but it's so weird!!!
    prev_edges.resize(mToLeft.size());
    next_edges.resize(mToRight.size());
    for (auto e=0;e<mToLeft.size();++e){

        uint64_t prev_node=mToLeft[e];

        prev_edges[e].resize(mHBV.ToSize(prev_node));
        for (int i=0;i<mHBV.ToSize(prev_node);++i){
            prev_edges[e][i]=mHBV.EdgeObjectIndexByIndexTo(prev_node,i);
        }

        uint64_t next_node=mToRight[e];

        next_edges[e].resize(mHBV.FromSize(next_node));
        for (int i=0;i<mHBV.FromSize(next_node);++i){
            next_edges[e][i]=mHBV.EdgeObjectIndexByIndexFrom(next_node,i);
        }
    }


}

std::array<uint64_t,3> PathFinder::transition_votes(uint64_t left_e,uint64_t right_e){
    std::array<uint64_t,3> tv={0,0,0};
    std::cout<<"Scoring transition from "<<left_e<<"to"<<right_e<<std::endl;

    std::cout<<"Paths on right edge "<<mEdgeToPathIds[right_e].size()<<" inv "<<mEdgeToPathIds[mInv[right_e]].size()<<std::endl;
    return tv;
};

std::array<uint64_t,3> PathFinder::path_votes(std::vector<uint64_t> path){
    //Returns a vote vector: { FOR, PARTIAL (reads end halfway, but validate at least one transition), AGAINST }
    //does this on forward and reverse paths, just in case
    std::vector<ReadPath> vfor,vpartial,vagainst;
    //TODO: needs to be done in both directions? (votes for only in one direction?)

    //std::cout<<std::endl<<"scoring path!!!!"<<std::endl;
    //first detect paths going out of first edge, also add them to open_paths
    std::list<std::pair<ReadPath,uint16_t>> initial_paths, open_paths;
    //std::cout<<"starting at edge "<<path[0]<<std::endl;
    for (auto pi:mEdgeToPathIds[path[0]]){

        auto p=mPaths[pi];
        //std::cout<<p<<std::endl;
        if (p.size()>1){
            uint16_t i=0;
            while (p[i]!=path[0]) ++i;
            if (i<p.size()-1) {
                open_paths.insert(open_paths.end(),std::make_pair(p,i));
                //std::cout<<" inserting into open paths"<<std::endl;
            }
        }
    }
    initial_paths.insert(initial_paths.begin(),open_paths.begin(),open_paths.end());
    //std::cout<<"initial paths generated from edge"<<path[0]<<" "<<initial_paths.size() <<std::endl;
    // basically every path in the mEdgeToPathIds[e] is either on the openPaths, starts here o
    for (auto ei=1;ei<path.size();++ei) {

        auto e=path[ei];
        //std::cout<<" going through edge "<<e<<" open list size: "<<open_paths.size()<<" paths on this edge: "<<mEdgeToPathIds[e].size()<<std::endl;
        //First go through the open list
        for (auto o=open_paths.begin();o!=open_paths.end();) {
            //if goes somewhere else, vote against and remove
            if (o->first[o->second+1]!=e){
                //std::cout<<"AGAINST: previous path goes to "<<o->first[o->second+1]<<std::endl;
                vagainst.push_back(o->first);
                o=open_paths.erase(o);
            } else { //else, advance
                ++(o->second);
                ++o;
            }
        }

        //std::cout<<" open list processed, votes" << pv[0]<<":"<<pv[1]<<":"<<pv[2]<<std::endl;

        std::list<std::pair<ReadPath,uint16_t>> new_paths;

        //check paths coming here from somewhere else and vote against
        for (auto ip:mEdgeToPathIds[e]) {
            //path is of size 1: irrelevant
            auto p=mPaths[ip];
            //std::cout<<"Considering path "<<ip<<":"<<p<<std::endl;
            if (p.size()==1) continue;

            //path starts_here, add to the open list later TODO: no need on the last edge!
            if (p[0]==e) {
                new_paths.insert(new_paths.end(),std::make_pair(p,0));
                //std::cout<<"adding new path to the open list"<<std::endl;
                continue;
            }

            //path in the open list?
            auto lp=open_paths.begin();
            for (;lp!=open_paths.end();++lp){
                if (p==lp->first) break;
            }

            if (lp!=open_paths.end()){
                //path is in the open list

                //last edge?
                if (ei==path.size()-1) {
                    //check if path in initial paths
                    auto ipp=initial_paths.begin();
                    for (;ipp!=initial_paths.end();++ipp){
                        if (p==ipp->first) break;
                    }
                    //path in initial_paths?
                    if (ipp!=initial_paths.end()){
                        //++pv[0];
                        vfor.push_back(p);
                        //std::cout<<"FOR: last edge and this path was in the initial list"<<std::endl;
                    }
                    else{
                        //++pv[1];
                        vpartial.push_back(p);
                        //std::cout<<"PARTIAL: last edge and this path was NOT in the initial list"<<std::endl;
                    }
                } else if (lp->first.size()-1==lp->second){
                    //++pv[1];
                    vpartial.push_back(p);
                    open_paths.erase(lp);
                    //std::cout<<"PARTIAL: path finished before the last edge"<<std::endl;
                }
            } else {
                //path comes from somewhere else, vote against
                //++pv[2];
                vagainst.push_back(p);
                //std::cout<<"AGAINST: path comes from somewhere else"<<std::endl;
            }

        }

        open_paths.insert(open_paths.end(),new_paths.begin(),new_paths.end());
    }

    //auto tv=transition_votes(path[i],path[i+1])
    //std::cout<<"Paths on left edge "<<mEdgeToPathIds[left_e].size()<<" inv "<<mEdgeToPathIds[mInv[left_e]].size()<<std::endl;
    std::array<uint64_t,3> pv={0,0,0};
    //collect votes:
    //on favour are on favour:
    //pv[0]=vfor.size();
    std::vector<ReadPath> votes_used;

    for (auto vf:vfor){
        bool u=false;
        for (auto vu:votes_used) if (vu.same_read(vf)) {u=true;break;}
        if (!u) {
            votes_used.insert(votes_used.end(), vf);
            pv[0]++;
        }
    }
    for (auto vp:vpartial){
        bool u=false;
        for (auto vu:votes_used) if (vu.same_read(vp)) {u=true;break;}
        if (!u) {
            votes_used.insert(votes_used.end(), vp);
            pv[1]++;
        }
    }
    for (auto va:vagainst){
        bool u=false;
        for (auto vu:votes_used) if (vu.same_read(va)) {u=true;break;}
        if (!u) {
            votes_used.insert(votes_used.end(), va);
            pv[2]++;
        }
    }
    return pv;
}
std::string PathFinder::path_str(std::vector<uint64_t> path) {
    std::string s="[";
    for (auto p:path){
        s+=std::to_string(p)+":"+std::to_string(mInv[p])+" ";//+" ("+std::to_string(mHBV.EdgeObject(p).size())+"bp "+std::to_string(paths_per_kbp(p))+"ppk)  ";
    }
    s+="]";
    return s;
}
std::array<uint64_t,3> PathFinder::multi_path_votes(std::vector<std::vector<uint64_t>> paths){
    //Returns a vote vector: { FOR, PARTIAL (reads end halfway, but validate at least one transition), AGAINST }
    //does this on forward and reverse paths, just in case
    std::vector<ReadPath> vfor,vpartial,vagainst;
    //TODO: needs to be done in both directions? (votes for only in one direction?)
    for (auto path:paths) {
        //std::cout<<std::endl<<"scoring path!!!!"<<std::endl;
        //first detect paths going out of first edge, also add them to open_paths
        std::list<std::pair<ReadPath, uint16_t>> initial_paths, open_paths;
        //std::cout<<"starting at edge "<<path[0]<<std::endl;
        for (auto pi:mEdgeToPathIds[path[0]]) {

            auto p = mPaths[pi];
            //std::cout<<p<<std::endl;
            if (p.size() > 1) {
                uint16_t i = 0;
                while (p[i] != path[0]) ++i;
                if (i < p.size() - 1) {
                    open_paths.insert(open_paths.end(), std::make_pair(p, i));
                    //std::cout<<" inserting into open paths"<<std::endl;
                }
            }
        }
        initial_paths.insert(initial_paths.begin(), open_paths.begin(), open_paths.end());
        //std::cout<<"initial paths generated from edge"<<path[0]<<" "<<initial_paths.size() <<std::endl;
        // basically every path in the mEdgeToPathIds[e] is either on the openPaths, starts here o
        for (auto ei = 1; ei < path.size(); ++ei) {

            auto e = path[ei];
            //std::cout<<" going through edge "<<e<<" open list size: "<<open_paths.size()<<" paths on this edge: "<<mEdgeToPathIds[e].size()<<std::endl;
            //First go through the open list
            for (auto o = open_paths.begin(); o != open_paths.end();) {
                //if goes somewhere else, vote against and remove
                if (o->first[o->second + 1] != e) {
                    //std::cout<<"AGAINST: previous path goes to "<<o->first[o->second+1]<<std::endl;
                    vagainst.push_back(o->first);
                    o = open_paths.erase(o);
                } else { //else, advance
                    ++(o->second);
                    ++o;
                }
            }

            //std::cout<<" open list processed, votes" << pv[0]<<":"<<pv[1]<<":"<<pv[2]<<std::endl;

            std::list<std::pair<ReadPath, uint16_t>> new_paths;

            //check paths coming here from somewhere else and vote against
            for (auto ip:mEdgeToPathIds[e]) {
                //path is of size 1: irrelevant
                auto p = mPaths[ip];
                //std::cout<<"Considering path "<<ip<<":"<<p<<std::endl;
                if (p.size() == 1) continue;

                //path starts_here, add to the open list later TODO: no need on the last edge!
                if (p[0] == e) {
                    new_paths.insert(new_paths.end(), std::make_pair(p, 0));
                    //std::cout<<"adding new path to the open list"<<std::endl;
                    continue;
                }

                //path in the open list?
                auto lp = open_paths.begin();
                for (; lp != open_paths.end(); ++lp) {
                    if (p == lp->first) break;
                }

                if (lp != open_paths.end()) {
                    //path is in the open list

                    //last edge?
                    if (ei == path.size() - 1) {
                        //check if path in initial paths
                        auto ipp = initial_paths.begin();
                        for (; ipp != initial_paths.end(); ++ipp) {
                            if (p == ipp->first) break;
                        }
                        //path in initial_paths?
                        if (ipp != initial_paths.end()) {
                            //++pv[0];
                            vfor.push_back(p);
                            //std::cout<<"FOR: last edge and this path was in the initial list"<<std::endl;
                        }
                        else {
                            //++pv[1];
                            vpartial.push_back(p);
                            //std::cout<<"PARTIAL: last edge and this path was NOT in the initial list"<<std::endl;
                        }
                    } else if (lp->first.size() - 1 == lp->second) {
                        //++pv[1];
                        vpartial.push_back(p);
                        open_paths.erase(lp);
                        //std::cout<<"PARTIAL: path finished before the last edge"<<std::endl;
                    }
                } else {
                    //path comes from somewhere else, vote against
                    //++pv[2];
                    vagainst.push_back(p);
                    //std::cout<<"AGAINST: path comes from somewhere else"<<std::endl;
                }

            }

            open_paths.insert(open_paths.end(), new_paths.begin(), new_paths.end());
        }
    }
    //auto tv=transition_votes(path[i],path[i+1])
    //std::cout<<"Paths on left edge "<<mEdgeToPathIds[left_e].size()<<" inv "<<mEdgeToPathIds[mInv[left_e]].size()<<std::endl;
    std::array<uint64_t,3> pv={0,0,0};
    //collect votes:
    //on favour are on favour:
    //pv[0]=vfor.size();
    std::vector<ReadPath> votes_used;

    for (auto vf:vfor){
        bool u=false;
        for (auto vu:votes_used) if (vu.same_read(vf)) {u=true;break;}
        if (!u) {
            votes_used.insert(votes_used.end(), vf);
            pv[0]++;
        }
    }
    for (auto vp:vpartial){
        bool u=false;
        for (auto vu:votes_used) if (vu.same_read(vp)) {u=true;break;}
        if (!u) {
            votes_used.insert(votes_used.end(), vp);
            pv[1]++;
        }
    }
    for (auto va:vagainst){
        bool u=false;
        for (auto vu:votes_used) if (vu.same_read(va)) {u=true;break;}
        if (!u) {
            votes_used.insert(votes_used.end(), va);
            pv[2]++;
        }
    }
    return pv;
}



void PathFinder::classify_forks(){
    int64_t nothing_fw=0, line_fw=0,join_fw=0,split_fw=0,join_split_fw=0;
    uint64_t nothing_fw_size=0, line_fw_size=0,join_fw_size=0,split_fw_size=0,join_split_fw_size=0;
    for ( int i = 0; i < mHBV.EdgeObjectCount(); ++i ) {
        //checks for both an edge and its complement, because it only looks forward
        uint64_t out_node=mToRight[i];
        if (mHBV.FromSize(out_node)==0) {
            nothing_fw++;
            nothing_fw_size+=mHBV.EdgeObject(i).size();
        } else if (mHBV.FromSize(out_node)==1) {
            if (mHBV.ToSize(out_node)==1){
                line_fw++;
                line_fw_size+=mHBV.EdgeObject(i).size();
            } else {
                split_fw++;
                split_fw_size+=mHBV.EdgeObject(i).size();
            }
        } else if (mHBV.ToSize(out_node)==1){
            join_fw++;
            join_fw_size+=mHBV.EdgeObject(i).size();
        } else {
            join_split_fw++;
            join_split_fw_size+=mHBV.EdgeObject(i).size();
        }
    }
    std::cout<<"Forward Node Edge Classification: "<<std::endl
    <<"nothing_fw: "<<nothing_fw<<" ( "<<nothing_fw_size<<" kmers )"<<std::endl
    <<"line_fw: "<<line_fw<<" ( "<<line_fw_size<<" kmers )"<<std::endl
    <<"join_fw: "<<join_fw<<" ( "<<join_fw_size<<" kmers )"<<std::endl
    <<"split_fw: "<<split_fw<<" ( "<<split_fw_size<<" kmers )"<<std::endl
    <<"join_split_fw: "<<join_split_fw<<" ( "<<join_split_fw_size<<" kmers )"<<std::endl;

}

void PathFinder::unroll_loops(uint64_t min_side_sizes) {
    //find nodes where in >1 or out>1 and in>0 and out>0
    uint64_t uloop=0,ursize=0;
    std::cout<<"Starting loop finding"<<std::endl;
    init_prev_next_vectors();
    std::cout<<"Prev and Next vectors initialised"<<std::endl;
    //score all possible transitions, discards all decidible and

        // is there any score>0 transition that is not incompatible with any other transitions?
    std::vector<std::vector<uint64_t>> new_paths; //these are solved paths, they will be materialised later
    for ( int e = 0; e < mHBV.EdgeObjectCount(); ++e ) {
        if (e<mInv[e]) {
            auto urs=is_unrollable_loop(e,min_side_sizes);

            auto iurs=is_unrollable_loop(mInv[e],min_side_sizes);
            if (urs.size()>0 && iurs.size()>0) {
                //std::cout<<"unrolling loop on edge"<<e<<std::endl;
                new_paths.push_back(urs[0]);
            }
        }

    }
    //std::cout<<"Unrollable loops: "<<uloop<<" ("<<ursize<<"bp)"<<std::endl;

    std::cout<<"Loop finding finished, "<<new_paths.size()<< " loops to unroll" <<std::endl;
    uint64_t sep=0;
    std::map<uint64_t,std::vector<uint64_t>> old_edges_to_new;
    for (auto p:new_paths){
        auto oen=separate_path(p);
        if (oen.size()>0) {
            for (auto et:oen){
                if (old_edges_to_new.count(et.first)==0) old_edges_to_new[et.first]={};
                for (auto ne:et.second) old_edges_to_new[et.first].push_back(ne);
            }
            sep++;
        }
    }
    if (old_edges_to_new.size()>0) {
        migrate_readpaths(old_edges_to_new);
    }
    std::cout<<sep<<" loops unrolled, re-initing the prev and next vectors, just in case :D"<<std::endl;
    init_prev_next_vectors();
    std::cout<<"Prev and Next vectors initialised"<<std::endl;
}
/*
void PathFinder::untangle_complex_in_out_choices() {
    //find a complex path
    init_prev_next_vectors();
    for (int e = 0; e < mHBV.EdgeObjectCount(); ++e) {
        if (e < mInv[e] && mHBV.EdgeObject(e).size()>1000) {
            //Ok, so why do we stop?
            //is next edge a join? how long till it splits again? can we choose the split?
            if (next_edges[e].size()==1 and prev_edges[next_edges[e][0]].size()>1){
                //std::cout<<"next edge from "<<e<<" is a join!"<<std::endl;
                if (mHBV.EdgeObject(next_edges[e][0]).size()<500 and next_edges[next_edges[e][0]].size()>1){

                    if (prev_edges[next_edges[e][0]].size()==prev_edges[next_edges[e][0]].size()){
                        std::cout<<"next edge from "<<e<<" is a small join! with "<<prev_edges[next_edges[e][0]].size()<<"in-outs"<<std::endl;
                        auto join_edge=next_edges[e][0];
                        std::vector<std::vector<uint64_t>> p0011={
                                {prev_edges[join_edge][0],join_edge,next_edges[join_edge][0]},
                                {prev_edges[join_edge][1],join_edge,next_edges[join_edge][1]}
                        };
                        std::vector<std::vector<uint64_t>> p1001={
                                {prev_edges[join_edge][1],join_edge,next_edges[join_edge][0]},
                                {prev_edges[join_edge][0],join_edge,next_edges[join_edge][1]}
                        };
                        auto v0011=multi_path_votes(p0011);
                        auto v1001=multi_path_votes(p1001);
                        std::cout<<"votes for: "<<path_str(p0011[0])<<" / "<<path_str(p0011[1])<<":  "<<v0011[0]<<":"<<v0011[1]<<":"<<v0011[2]<<std::endl;
                        std::cout<<"votes for: "<<path_str(p1001[0])<<" / "<<path_str(p1001[1])<<":  "<<v1001[0]<<":"<<v1001[1]<<":"<<v1001[2]<<std::endl;
                    }
                }
            }

        }
    }
}*/

void PathFinder::untangle_pins() {

    init_prev_next_vectors();
    uint64_t pins=0;
    for (int e = 0; e < mHBV.EdgeObjectCount(); ++e) {
        if (mToLeft[e]==mToLeft[mInv[e]] and next_edges[e].size()==1 ) {
            std::cout<<" Edge "<<e<<" forms a pinhole!!!"<<std::endl;
            if (next_edges[next_edges[e][0]].size()==2) {
                std::vector<uint64_t> pfw = {mInv[next_edges[next_edges[e][0]][0]],mInv[next_edges[e][0]],e,next_edges[e][0],next_edges[next_edges[e][0]][1]};
                std::vector<uint64_t> pbw = {mInv[next_edges[next_edges[e][0]][1]],mInv[next_edges[e][0]],e,next_edges[e][0],next_edges[next_edges[e][0]][0]};
                auto vpfw=multi_path_votes({pfw});
                auto vpbw=multi_path_votes({pbw});
                std::cout<<"votes FW: "<<vpfw[0]<<":"<<vpfw[1]<<":"<<vpfw[2]<<"     BW: "<<vpbw[0]<<":"<<vpbw[1]<<":"<<vpbw[2]<<std::endl;
            }
            ++pins;
        }
    }
    std::cout<<"Total number of pinholes: "<<pins;
}

void PathFinder::untangle_complex_in_out_choices(uint64_t large_frontier_size, bool verbose_separation) {
    //find a complex path
    uint64_t qsf=0,qsf_paths=0;
    uint64_t msf=0,msf_paths=0;
    init_prev_next_vectors();
    std::cout<<"vectors initialised"<<std::endl;
    std::set<std::array<std::vector<uint64_t>,2>> seen_frontiers,solved_frontiers;
    std::vector<std::vector<uint64_t>> paths_to_separate;
    for (int e = 0; e < mHBV.EdgeObjectCount(); ++e) {
        if (e < mInv[e] && mHBV.EdgeObject(e).size() < large_frontier_size) {
            auto f=get_all_long_frontiers(e, large_frontier_size);
            if (f[0].size()>1 and f[1].size()>1 and seen_frontiers.count(f)==0){
                seen_frontiers.insert(f);
                bool single_dir=true;
                for (auto in_e:f[0]) for (auto out_e:f[1]) if (in_e==out_e) {single_dir=false;break;}
                if (single_dir) {
                    std::cout<<" Single direction frontiers for complex region on edge "<<e<<" IN:"<<path_str(f[0])<<" OUT: "<<path_str(f[1])<<std::endl;
                    std::vector<int> in_used(f[0].size(),0);
                    std::vector<int> out_used(f[1].size(),0);
                    std::vector<std::vector<uint64_t>> first_full_paths;
                    bool reversed=false;
                    for (auto in_i=0;in_i<f[0].size();++in_i) {
                        auto in_e=f[0][in_i];
                        for (auto out_i=0;out_i<f[1].size();++out_i) {
                            auto out_e=f[1][out_i];
                            auto shared_paths = 0;

                            for (auto inp:mEdgeToPathIds[in_e])
                                for (auto outp:mEdgeToPathIds[out_e])
                                    if (inp == outp) {

                                        shared_paths++;
                                        if (shared_paths==1){//not the best solution, but should work-ish
                                            std::vector<uint64_t > pv;
                                            for (auto e:mPaths[inp]) pv.push_back(e);
                                            std::cout<<"found first path from "<<in_e<<" to "<< out_e << path_str(pv)<< std::endl;
                                            first_full_paths.push_back({});
                                            int16_t ei=0;
                                            while (mPaths[inp][ei]!=in_e) ei++;

                                            while (mPaths[inp][ei]!=out_e && ei<mPaths[inp].size()) first_full_paths.back().push_back(mPaths[inp][ei++]);
                                            if (ei>=mPaths[inp].size()) {
                                                std::cout<<"reversed path detected!"<<std::endl;
                                                reversed=true;
                                            }
                                            first_full_paths.back().push_back(out_e);
                                            //std::cout<<"added!"<<std::endl;
                                        }
                                    }
                            //check for reverse paths too
                            for (auto inp:mEdgeToPathIds[mInv[out_e]])
                                for (auto outp:mEdgeToPathIds[mInv[in_e]])
                                    if (inp == outp) {

                                        shared_paths++;
                                        if (shared_paths==1){//not the best solution, but should work-ish
                                            std::vector<uint64_t > pv;
                                            for (auto e=mPaths[inp].rbegin();e!=mPaths[inp].rend();++e) pv.push_back(mInv[*e]);
                                            std::cout<<"found first path from "<<in_e<<" to "<< out_e << path_str(pv)<< std::endl;
                                            first_full_paths.push_back({});
                                            int16_t ei=0;
                                            while (pv[ei]!=in_e) ei++;

                                            while (pv[ei]!=out_e && ei<pv.size()) first_full_paths.back().push_back(pv[ei++]);
                                            if (ei>=pv.size()) {
                                                std::cout<<"reversed path detected!"<<std::endl;
                                                reversed=true;
                                            }
                                            first_full_paths.back().push_back(out_e);
                                            //std::cout<<"added!"<<std::endl;
                                        }
                                    }
                            if (shared_paths) {
                                out_used[out_i]++;
                                in_used[in_i]++;
                                //std::cout << "  Shared paths " << in_e << " --> " << out_e << ": " << shared_paths << std::endl;

                            }
                        }
                    }
                    if ((not reversed) and std::count(in_used.begin(),in_used.end(),1) == in_used.size() and
                            std::count(out_used.begin(),out_used.end(),1) == out_used.size()){
                        std::cout<<" REGION COMPLETELY SOLVED BY PATHS!!!"<<std::endl;
                        solved_frontiers.insert(f);
                        for (auto p:first_full_paths) paths_to_separate.push_back(p);
                    } /*else if (std::count(in_used.begin(),in_used.end(),1) == in_used.size()-1 and
                            std::count(in_used.begin(),in_used.end(),0) == 1 and
                            std::count(out_used.begin(),out_used.end(),1) == out_used.size()-1 and
                            std::count(out_used.begin(),out_used.end(),0) == 1){
                        //std::cout<<" REGION SOLVED BY PATHS and a jump (not acted on!!!)"<<std::endl;
                        //solved_frontiers.insert(f);
                        qsf++;
                        qsf_paths+=in_used.size();
                        unsigned int in_index=0;
                        while (in_used[in_index]!=0) in_index++;
                        unsigned int out_index=0;
                        while (out_used[out_index]!=0) out_index++;
                        std::cout<<"Trying to solve region by reducing last choice to unique path between "<<f[0][in_index]<< " and "<< f[1][out_index]<<std::endl;

                        auto all_paths=AllPathsFromTo({f[0][in_index]},{f[1][out_index]},10);
                        if (all_paths.size()==1){
                            solved_frontiers.insert(f);
                            for (auto p:first_full_paths) paths_to_separate.push_back(p);
                            paths_to_separate.push_back(all_paths[0]);

                            std::cout<<"Solved!!!"<<std::endl;
                        }
                        else std::cout<<"Not solved, "<<all_paths.size()<<" different paths :("<<std::endl;

                    } else if (std::count(in_used.begin(),in_used.end(),0) == 0 and
                        std::count(out_used.begin(),out_used.end(),0) == 0){
                        msf++;
                        msf_paths+=in_used.size();
                    }*/

                }

            }
        }
    }
    std::cout<<"Complex Regions solved by paths: "<<solved_frontiers.size() <<"/"<<seen_frontiers.size()<<" comprising "<<paths_to_separate.size()<<" paths to separate"<< std::endl;
    //std::cout<<"Complex Regions quasi-solved by paths (not acted on): "<< qsf <<"/"<<seen_frontiers.size()<<" comprising "<<qsf_paths<<" paths to separate"<< std::endl;
    //std::cout<<"Multiple Solution Regions (not acted on): "<< msf <<"/"<<seen_frontiers.size()<<" comprising "<<msf_paths<<" paths to separate"<< std::endl;

    uint64_t sep=0;
    std::map<uint64_t,std::vector<uint64_t>> old_edges_to_new;
    for (auto p:paths_to_separate){

        if (old_edges_to_new.count(p.front()) > 0 or old_edges_to_new.count(p.back()) > 0) {
            std::cout<<"WARNING: path starts or ends in an already modified edge, skipping"<<std::endl;
            continue;
        }

        auto oen=separate_path(p, verbose_separation);
        if (oen.size()>0) {
            for (auto et:oen){
                if (old_edges_to_new.count(et.first)==0) old_edges_to_new[et.first]={};
                for (auto ne:et.second) old_edges_to_new[et.first].push_back(ne);
            }
            sep++;
        }
    }
    if (old_edges_to_new.size()>0) {
        migrate_readpaths(old_edges_to_new);
    }
    std::cout<<" "<<sep<<" paths separated!"<<std::endl;
}

std::vector<std::vector<uint64_t>> PathFinder::AllPathsFromTo(std::vector<uint64_t> in_edges, std::vector<uint64_t> out_edges, uint64_t max_length) {
    std::vector<std::vector<uint64_t>> current_paths;
    std::vector<std::vector<uint64_t>> paths;
    //First, start all paths from in_edges
    for (auto ine:in_edges) current_paths.push_back({ine});
    while (max_length--) {
        auto old_paths=current_paths;
        current_paths.clear();
        for (auto op:old_paths) {
            //grow each path, adding variations if needed
            for (auto ne:next_edges[op.back()]){
                //if new edge on out_edges, add to paths
                op.push_back(ne);
                if (std::count(out_edges.begin(),out_edges.end(),ne)) paths.push_back(op);
                else current_paths.push_back(op);
            }
        }
    }
    return paths;

}

std::array<std::vector<uint64_t>,2> PathFinder::get_all_long_frontiers(uint64_t e, uint64_t large_frontier_size){
    //TODO: return all components in the region
    std::set<uint64_t> seen_edges, to_explore={e}, in_frontiers, out_frontiers;

    while (to_explore.size()>0){
        std::set<uint64_t> next_to_explore;
        for (auto x:to_explore){ //to_explore: replace rather and "update" (use next_to_explore)

            if (!seen_edges.count(x)){

                //What about reverse complements and paths that include loops that "reverse the flow"?
                if (seen_edges.count(mInv[x])) return std::array<std::vector<uint64_t>,2>(); //just cancel for now

                for (auto p:prev_edges[x]) {
                    if (mHBV.EdgeObject(p).size() >= large_frontier_size )  {
                        //What about frontiers on both sides?
                        in_frontiers.insert(p);
                        for (auto other_n:next_edges[p]){
                            if (!seen_edges.count(other_n)) next_to_explore.insert(other_n);
                        }
                    }
                    else if (!seen_edges.count(p)) next_to_explore.insert(p);
                }
                for (auto n:next_edges[x]) {
                    if (mHBV.EdgeObject(n).size() >= large_frontier_size) {
                        //What about frontiers on both sides?
                        out_frontiers.insert(n);
                        for (auto other_p:prev_edges[n]){
                            if (!seen_edges.count(other_p)) next_to_explore.insert(other_p);
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

std::vector<std::vector<uint64_t>> PathFinder::is_unrollable_loop(uint64_t loop_e, uint64_t min_size_sizes){
    //Checks if loop can be unrolled
    //Input: looping edge (i.e. prev_e--->R---->loop_e---->R---->next_e)
    uint64_t prev_e, repeat_e, next_e;
    //Conditions for unrollable loop:

    //1) only one neighbour on each direction, and the same one (repeat_e).
    if (prev_edges[loop_e].size()!=1 or
        next_edges[loop_e].size()!=1 or
        prev_edges[loop_e][0]!=next_edges[loop_e][0]) return {};
    repeat_e=prev_edges[loop_e][0];


    //2) the repeat edge has only one other neighbour on each direction, and it is a different one;
    if (prev_edges[repeat_e].size()!=2 or
        next_edges[repeat_e].size()!=2) return {};

    prev_e=(prev_edges[repeat_e][0]==loop_e ? prev_edges[repeat_e][1]:prev_edges[repeat_e][0]);

    next_e=(next_edges[repeat_e][0]==loop_e ? next_edges[repeat_e][1]:next_edges[repeat_e][0]);

    if (prev_e==next_e or prev_e==mInv[next_e]) return {};

    //3) size constraints: prev_e and next_e must be at least 1Kbp
    if (mHBV.EdgeObject(prev_e).size()<min_size_sizes or mHBV.EdgeObject(next_e).size()<min_size_sizes) return {};


    //std::cout<<" LOOP: "<< edge_pstr(prev_e)<<" ---> "<<edge_pstr(repeat_e)<<" -> "<<edge_pstr(loop_e)<<" -> "<<edge_pstr(repeat_e)<<" ---> "<<edge_pstr(next_e)<<std::endl;
    auto pvlin=path_votes({prev_e,repeat_e,loop_e,repeat_e,next_e});
    auto pvloop=path_votes({prev_e,repeat_e,loop_e,repeat_e,loop_e,repeat_e,next_e});
    //std::cout<<" Votes to p->r->l->r->n: "<<pvlin[0]<<":"<<pvlin[1]<<":"<<pvlin[2]<<std::endl;
    //std::cout<<" Votes to p->r->l->r->l->r->n: "<<pvloop[0]<<":"<<pvloop[1]<<":"<<pvloop[2]<<std::endl;

    auto pvcircleline=multi_path_votes({{loop_e,repeat_e,loop_e},{prev_e,repeat_e,next_e}});
    //std::cout<<" Multi-votes to circle+line "<<pvcircleline[0]<<":"<<pvcircleline[1]<<":"<<pvcircleline[2]<<std::endl;

    if (pvcircleline[0]>0 or pvloop[2]>0 or (pvcircleline[2]==0 and pvcircleline[1]>pvlin[1] and pvcircleline[1]>pvloop[1])) {
        //std::cout<<"   CAN'T be reliably unrolled!"<<std::endl;
        return {};
    }
    if (pvloop[3]==0 and pvloop[0]>pvlin[0]){
        //std::cout<<"   LOOP should be traversed as loop at least once!"<<std::endl;
        return {};
    }
    if (pvlin==pvcircleline){
        //std::cout<<"   UNDECIDABLE as path support problem, looking at coverages"<<std::endl;
        float prev_cov=paths_per_kbp(prev_e);
        float repeat_cov=paths_per_kbp(repeat_e);
        float loop_cov=paths_per_kbp(loop_e);
        float next_cov=paths_per_kbp(next_e);
        auto sc_min=prev_cov*.8;
        auto sc_max=prev_cov*1.2;
        auto dc_min=prev_cov*1.8;
        auto dc_max=prev_cov*2.2;
        if (repeat_cov<dc_min or repeat_cov>dc_max or
                loop_cov<sc_min or loop_cov>sc_max or
                next_cov<sc_min or next_cov>sc_max ){
            //std::cout<<"    Coverage FAIL conditions, CAN'T be reliably unrolled!"<<std::endl;
            return {};
        }
        //std::cout<<"    Coverage OK"<<std::endl;
    }
    //std::cout<<"   UNROLLABLE"<<std::endl;
    return {{prev_e,repeat_e,loop_e,repeat_e,next_e}};
    //return mHBV.EdgeObject(prev_e).size()+mHBV.EdgeObject(repeat_e).size()+
            //mHBV.EdgeObject(loop_e).size()+mHBV.EdgeObject(repeat_e).size()+mHBV.EdgeObject(next_e).size()-4*mHBV.K();
}

uint64_t PathFinder::paths_per_kbp(uint64_t e){
    return 1000 * mEdgeToPathIds[e].size()/mHBV.EdgeObject(e).size();
};

std::map<uint64_t,std::vector<uint64_t>> PathFinder::separate_path(std::vector<uint64_t> p, bool verbose_separation){

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
    uint64_t current_vertex_fw=mHBV.N(),current_vertex_rev=mHBV.N()+1;
    mHBV.AddVertices(2);
    //migrate connections (dangerous!!!)
    if (verbose_separation) std::cout<<"Migrating edge "<<p[0]<<" To node old: "<<mToRight[p[0]]<<" new: "<<current_vertex_fw<<std::endl;
    mHBV.GiveEdgeNewToVx(p[0],mToRight[p[0]],current_vertex_fw);
    if (verbose_separation) std::cout<<"Migrating edge "<<mInv[p[0]]<<" From node old: "<<mToLeft[mInv[p[0]]]<<" new: "<<current_vertex_rev<<std::endl;
    mHBV.GiveEdgeNewFromVx(mInv[p[0]],mToLeft[mInv[p[0]]],current_vertex_rev);
    std::map<uint64_t,std::vector<uint64_t>> old_edges_to_new;

    for (auto ei=1;ei<p.size()-1;++ei){
        //add a new vertex for each of FW and BW paths
        uint64_t prev_vertex_fw=current_vertex_fw,prev_vertex_rev=current_vertex_rev;
        //create two new vertices (for the FW and BW path)
        current_vertex_fw=mHBV.N();
        current_vertex_rev=mHBV.N()+1;
        mHBV.AddVertices(2);

        //now, duplicate next edge for the FW and reverse path
        auto nef=mHBV.AddEdge(prev_vertex_fw,current_vertex_fw,mHBV.EdgeObject(p[ei]));
        if (verbose_separation)  std::cout<<"Edge "<<nef<<": copy of "<<p[ei]<<": "<<prev_vertex_fw<<" - "<<current_vertex_fw<<std::endl;
        mToLeft.push_back(prev_vertex_fw);
        mToRight.push_back(current_vertex_fw);
        if (! old_edges_to_new.count(p[ei]))  old_edges_to_new[p[ei]]={};
        old_edges_to_new[p[ei]].push_back(nef);

        auto ner=mHBV.AddEdge(current_vertex_rev,prev_vertex_rev,mHBV.EdgeObject(mInv[p[ei]]));
        if (verbose_separation) std::cout<<"Edge "<<ner<<": copy of "<<mInv[p[ei]]<<": "<<current_vertex_rev<<" - "<<prev_vertex_rev<<std::endl;
        mToLeft.push_back(current_vertex_rev);
        mToRight.push_back(prev_vertex_rev);
        if (! old_edges_to_new.count(mInv[p[ei]]))  old_edges_to_new[mInv[p[ei]]]={};
        old_edges_to_new[mInv[p[ei]]].push_back(ner);

        mInv.push_back(ner);
        mInv.push_back(nef);
        mEdgeToPathIds.resize(mEdgeToPathIds.size()+2);
    }
    if (verbose_separation) std::cout<<"Migrating edge "<<p[p.size()-1]<<" From node old: "<<mToLeft[p[p.size()-1]]<<" new: "<<current_vertex_fw<<std::endl;
    mHBV.GiveEdgeNewFromVx(p[p.size()-1],mToLeft[p[p.size()-1]],current_vertex_fw);
    if (verbose_separation) std::cout<<"Migrating edge "<<mInv[p[p.size()-1]]<<" To node old: "<<mToRight[mInv[p[p.size()-1]]]<<" new: "<<current_vertex_rev<<std::endl;
    mHBV.GiveEdgeNewToVx(mInv[p[p.size()-1]],mToRight[mInv[p[p.size()-1]]],current_vertex_rev);

    //TODO: cleanup new isolated elements and leading-nowhere paths.
    //for (auto ei=1;ei<p.size()-1;++ei) mHBV.DeleteEdges({p[ei]});
    return old_edges_to_new;

}

//TODO: this should probably be called just once
void PathFinder::migrate_readpaths(std::map<uint64_t,std::vector<uint64_t>> edgemap){
    //Migrate readpaths: this changes the readpaths from old edges to new edges
    //if an old edge has more than one new edge it tries all combinations until it gets the paths to map
    //if more than one combination is valid, this chooses at random among them (could be done better? should the path be duplicated?)
    mHBV.ToLeft(mToLeft);
    mHBV.ToRight(mToRight);
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

bool PathFinder::join_edges_in_path(std::vector<uint64_t> p){
    //if a path is scrictly linear, it joins it
}