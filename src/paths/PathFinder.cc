//
// Created by Bernardo Clavijo (TGAC) on 11/07/2016.
//

#include <util/OutputLog.h>
#include "PathFinder.h"

void PathFinder::clip_tips(int max_length, int max_reads) {
    update_prev_next();
    vec<int> to_delete;
    for (uint64_t e=0;e<next_edges.size();++e){
        if (next_edges[e].empty() and prev_edges[e].size()==1 and mHBV.Edges()[e].size()<=max_length){
            //std::cout<<"Edge "<<e<<" ("<<mInv[e]<<") leads nowhere!!!"<<std::endl;
            //now check for support
            //std::cout<<" Support: "<< mEdgeToPathIds[e].size()<<" "<<mEdgeToPathIds[e].size()<<std::endl;
            if (mEdgeToPathIds[e].size()<=max_reads and mEdgeToPathIds[e].size()<=max_reads) {
                to_delete.push_back(e);
                to_delete.push_back(mInv[e]);
            }
        }
    }
    mHBV.DeleteEdges(to_delete);
    Cleanup( mHBV, mInv, mPaths );
    invert(mPaths,mEdgeToPathIds,mHBV.EdgeObjectCount());
    update_prev_next();
}

void PathFinder::extend_until_repeat(int max_collapsed_size) {
    update_prev_next();
    vec<int> to_delete;
    for (uint64_t e=0;e<next_edges.size();++e){
        if (next_edges[e].size()==1 and prev_edges[e].size()>1){
            if (mHBV.EdgeObject(e).size()<=max_collapsed_size) {
                std::cout << "Edge " << e << " (" << mInv[e] << ") can be back-filled into its inputs" << std::endl;
            }
        }
    }

}

std::string PathFinder::edge_pstr(uint64_t e){
    return "e"+std::to_string(e)+"("+std::to_string(mHBV.EdgeObject(e).size())+"bp "+std::to_string(paths_per_kbp(e))+"ppk)";
};

void PathFinder::update_prev_next(){
    //TODO: this is stupid duplication of the digraph class, but it's so weird!!!
    mHBV.ToLeft(mToLeft);
    mHBV.ToRight(mToRight);
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

void PathFinder::extend_bridging_paths() {
    //takes each pair of paths and if they connect consecutive edges it just extend the forward and reverse paths to include the other part
    if (mVerbose) std::cout<<Date()<<": pathfinder extending concurrent paths"<<std::endl;
    std::atomic_uint_fast64_t pairs_same(0),pairs_both_unplaced(0),pairs_single(0),pairs_overlap_contained(0),pairs_overlap_extendable(0),pairs_neighbour(0), pairs_complex(0);
    for (auto i=0;i<mPaths.size(); i+=2){
        auto r1path=mPaths[i];
        ReadPath r2rpath;
        for (int ei=mPaths[i+1].size()-1;ei>=0;--ei) r2rpath.push_back(mInv[mPaths[i+1][ei]]);
        if (r1path.size()==0 and r2rpath.size()==0){
            pairs_both_unplaced++;
        }
        else if (r1path.size()==0 or r2rpath.size()==0) {
            pairs_single++;
        }
        else if (r1path==r2rpath) {
            pairs_same++;
        }
        else {
            //first hipothesis:paths overlap
            auto r1last = r1path[r1path.size() - 1];
            int overlap=-1;
            int r2i=0;
            for (auto &e:r2rpath){
                if (overlap==-1 and r1last==e){
                    overlap=r2i;
                    break;
                }
                r2i++;
            }
            if (overlap!=-1){
                //if (overlap+1<r2rpath.size() or overlap<r1path.size()-1) {
                if (overlap+1<r2rpath.size() and overlap+1<r1path.size()) {
                    pairs_overlap_extendable++;
                    /*std::cout << "path overlap: ";
                    for (auto &e:r1path) std::cout << e << "|" << mInv[e] << "(" << mHBV.EdgeObject(e).size() - mHBV.K() + 1 << ") ";
                    std::cout << "  -> <-  ";
                    for (auto &e:r2rpath) std::cout << e << "|" << mInv[e] << "(" << mHBV.EdgeObject(e).size() - mHBV.K() + 1 << ") ";
                    std::cout << std::endl;*/
                    for (auto r2ai=overlap; r2ai<r2rpath.size(); ++r2ai) r1path.push_back(r2rpath[r2ai]);
                    mPaths[i]=r1path;
                    r2rpath.clear();
                    for (int ei=r1path.size()-1;ei>=0;--ei) mPaths[i+1].push_back(mInv[r1path[ei]]);

                }
                else pairs_overlap_contained++;
                //}
                //else pairs_complex++;
            }
            else if (mToRight[r1last]==r2rpath[0]) pairs_neighbour++; //second hipothesis, r1last and r2rfirst are perfect neighbours
            else pairs_complex++;

        }
    }
    if (mVerbose) std::cout<<Date()<<": patfinder pairs: unplaced: " <<
                           pairs_both_unplaced << "  single: " << pairs_single << "  same-edge: " << pairs_same <<
                           "  overlap(cont): " << pairs_overlap_contained << "  overlap(ext): " <<
                           pairs_overlap_extendable << "  neighbour: " << pairs_neighbour << "  others: " <<
                           pairs_complex <<std::endl;

}
std::array<uint64_t,3> PathFinder::transition_votes(uint64_t left_e,uint64_t right_e){
    std::array<uint64_t,3> tv={0,0,0};
    if (mVerbose) std::cout<<"Scoring transition from "<<left_e<<"to"<<right_e<<std::endl;

    if (mVerbose) std::cout<<"Paths on right edge "<<mEdgeToPathIds[right_e].size()<<" inv "<<mEdgeToPathIds[mInv[right_e]].size()<<std::endl;
    return tv;
};

std::array<uint64_t,3> PathFinder::path_votes(std::vector<uint64_t> path){
    //Returns a vote vector: { FOR, PARTIAL (reads end halfway, but validate at least one transition), AGAINST }
    //does this on forward and reverse paths, just in case
    std::vector<ReadPath> vfor,vpartial,vagainst;
    //TODO: needs to be done in both directions? (votes for only in one direction?)

    //if (mVerbose) std::cout<<std::endl<<"scoring path!!!!"<<std::endl;
    //first detect paths going out of first edge, also add them to open_paths
    std::list<std::pair<ReadPath,uint16_t>> initial_paths, open_paths;
    //if (mVerbose) std::cout<<"starting at edge "<<path[0]<<std::endl;
    for (auto pi:mEdgeToPathIds[path[0]]){

        auto p=mPaths[pi];
        //if (mVerbose) std::cout<<p<<std::endl;
        if (p.size()>1){
            uint16_t i=0;
            while (p[i]!=path[0]) ++i;
            if (i<p.size()-1) {
                open_paths.insert(open_paths.end(),std::make_pair(p,i));
                //if (mVerbose) std::cout<<" inserting into open paths"<<std::endl;
            }
        }
    }
    initial_paths.insert(initial_paths.begin(),open_paths.begin(),open_paths.end());
    //if (mVerbose) std::cout<<"initial paths generated from edge"<<path[0]<<" "<<initial_paths.size() <<std::endl;
    // basically every path in the mEdgeToPathIds[e] is either on the openPaths, starts here o
    for (auto ei=1;ei<path.size();++ei) {

        auto e=path[ei];
        //if (mVerbose) std::cout<<" going through edge "<<e<<" open list size: "<<open_paths.size()<<" paths on this edge: "<<mEdgeToPathIds[e].size()<<std::endl;
        //First go through the open list
        for (auto o=open_paths.begin();o!=open_paths.end();) {
            //if goes somewhere else, vote against and remove
            if (o->first[o->second+1]!=e){
                //if (mVerbose) std::cout<<"AGAINST: previous path goes to "<<o->first[o->second+1]<<std::endl;
                vagainst.push_back(o->first);
                o=open_paths.erase(o);
            } else { //else, advance
                ++(o->second);
                ++o;
            }
        }

        //if (mVerbose) std::cout<<" open list processed, votes" << pv[0]<<":"<<pv[1]<<":"<<pv[2]<<std::endl;

        std::list<std::pair<ReadPath,uint16_t>> new_paths;

        //check paths coming here from somewhere else and vote against
        for (auto ip:mEdgeToPathIds[e]) {
            //path is of size 1: irrelevant
            auto p=mPaths[ip];
            //if (mVerbose) std::cout<<"Considering path "<<ip<<":"<<p<<std::endl;
            if (p.size()==1) continue;

            //path starts_here, add to the open list later TODO: no need on the last edge!
            if (p[0]==e) {
                new_paths.insert(new_paths.end(),std::make_pair(p,0));
                //if (mVerbose) std::cout<<"adding new path to the open list"<<std::endl;
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
                        //if (mVerbose) std::cout<<"FOR: last edge and this path was in the initial list"<<std::endl;
                    }
                    else{
                        //++pv[1];
                        vpartial.push_back(p);
                        //if (mVerbose) std::cout<<"PARTIAL: last edge and this path was NOT in the initial list"<<std::endl;
                    }
                } else if (lp->first.size()-1==lp->second){
                    //++pv[1];
                    vpartial.push_back(p);
                    open_paths.erase(lp);
                    //if (mVerbose) std::cout<<"PARTIAL: path finished before the last edge"<<std::endl;
                }
            } else {
                //path comes from somewhere else, vote against
                //++pv[2];
                vagainst.push_back(p);
                //if (mVerbose) std::cout<<"AGAINST: path comes from somewhere else"<<std::endl;
            }

        }

        open_paths.insert(open_paths.end(),new_paths.begin(),new_paths.end());
    }

    //auto tv=transition_votes(path[i],path[i+1])
    //if (mVerbose) std::cout<<"Paths on left edge "<<mEdgeToPathIds[left_e].size()<<" inv "<<mEdgeToPathIds[mInv[left_e]].size()<<std::endl;
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
        //if (mVerbose) std::cout<<std::endl<<"scoring path!!!!"<<std::endl;
        //first detect paths going out of first edge, also add them to open_paths
        std::list<std::pair<ReadPath, uint16_t>> initial_paths, open_paths;
        //if (mVerbose) std::cout<<"starting at edge "<<path[0]<<std::endl;
        for (auto pi:mEdgeToPathIds[path[0]]) {

            auto p = mPaths[pi];
            //if (mVerbose) std::cout<<p<<std::endl;
            if (p.size() > 1) {
                uint16_t i = 0;
                while (p[i] != path[0]) ++i;
                if (i < p.size() - 1) {
                    open_paths.insert(open_paths.end(), std::make_pair(p, i));
                    //if (mVerbose) std::cout<<" inserting into open paths"<<std::endl;
                }
            }
        }
        initial_paths.insert(initial_paths.begin(), open_paths.begin(), open_paths.end());
        //if (mVerbose) std::cout<<"initial paths generated from edge"<<path[0]<<" "<<initial_paths.size() <<std::endl;
        // basically every path in the mEdgeToPathIds[e] is either on the openPaths, starts here o
        for (auto ei = 1; ei < path.size(); ++ei) {

            auto e = path[ei];
            //if (mVerbose) std::cout<<" going through edge "<<e<<" open list size: "<<open_paths.size()<<" paths on this edge: "<<mEdgeToPathIds[e].size()<<std::endl;
            //First go through the open list
            for (auto o = open_paths.begin(); o != open_paths.end();) {
                //if goes somewhere else, vote against and remove
                if (o->first[o->second + 1] != e) {
                    //if (mVerbose) std::cout<<"AGAINST: previous path goes to "<<o->first[o->second+1]<<std::endl;
                    vagainst.push_back(o->first);
                    o = open_paths.erase(o);
                } else { //else, advance
                    ++(o->second);
                    ++o;
                }
            }

            //if (mVerbose) std::cout<<" open list processed, votes" << pv[0]<<":"<<pv[1]<<":"<<pv[2]<<std::endl;

            std::list<std::pair<ReadPath, uint16_t>> new_paths;

            //check paths coming here from somewhere else and vote against
            for (auto ip:mEdgeToPathIds[e]) {
                //path is of size 1: irrelevant
                auto p = mPaths[ip];
                //if (mVerbose) std::cout<<"Considering path "<<ip<<":"<<p<<std::endl;
                if (p.size() == 1) continue;

                //path starts_here, add to the open list later TODO: no need on the last edge!
                if (p[0] == e) {
                    new_paths.insert(new_paths.end(), std::make_pair(p, 0));
                    //if (mVerbose) std::cout<<"adding new path to the open list"<<std::endl;
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
                            //if (mVerbose) std::cout<<"FOR: last edge and this path was in the initial list"<<std::endl;
                        }
                        else {
                            //++pv[1];
                            vpartial.push_back(p);
                            //if (mVerbose) std::cout<<"PARTIAL: last edge and this path was NOT in the initial list"<<std::endl;
                        }
                    } else if (lp->first.size() - 1 == lp->second) {
                        //++pv[1];
                        vpartial.push_back(p);
                        open_paths.erase(lp);
                        //if (mVerbose) std::cout<<"PARTIAL: path finished before the last edge"<<std::endl;
                    }
                } else {
                    //path comes from somewhere else, vote against
                    //++pv[2];
                    vagainst.push_back(p);
                    //if (mVerbose) std::cout<<"AGAINST: path comes from somewhere else"<<std::endl;
                }

            }

            open_paths.insert(open_paths.end(), new_paths.begin(), new_paths.end());
        }
    }
    //auto tv=transition_votes(path[i],path[i+1])
    //if (mVerbose) std::cout<<"Paths on left edge "<<mEdgeToPathIds[left_e].size()<<" inv "<<mEdgeToPathIds[mInv[left_e]].size()<<std::endl;
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





void PathFinder::unroll_loops(uint64_t min_side_sizes) {
    //extend_bridging_paths();
    //find nodes where in >1 or out>1 and in>0 and out>0
    uint64_t uloop=0,ursize=0;
    if (mVerbose) std::cout<<"Starting loop finding"<<std::endl;
    update_prev_next();
    if (mVerbose) std::cout<<"Prev and Next vectors initialised"<<std::endl;
    //score all possible transitions, discards all decidible and

        // is there any score>0 transition that is not incompatible with any other transitions?
    std::vector<std::vector<uint64_t>> new_paths; //these are solved paths, they will be materialised later
    for ( int e = 0; e < mHBV.EdgeObjectCount(); ++e ) {
        if (e<mInv[e]) {
            auto urs=is_unrollable_loop(e,min_side_sizes);

            auto iurs=is_unrollable_loop(mInv[e],min_side_sizes);
            if (urs.size()>0 && iurs.size()>0) {
                //if (mVerbose) std::cout<<"unrolling loop on edge"<<e<<std::endl;
                new_paths.push_back(urs[0]);
            }
        }

    }
    //if (mVerbose) std::cout<<"Unrollable loops: "<<uloop<<" ("<<ursize<<"bp)"<<std::endl;

    if (mVerbose) std::cout<<"Loop finding finished, "<<new_paths.size()<< " loops to unroll" <<std::endl;
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
    if (mVerbose) std::cout<<sep<<" loops unrolled, re-initing the prev and next vectors, just in case :D"<<std::endl;
    update_prev_next();
    if (mVerbose) std::cout<<"Prev and Next vectors updated"<<std::endl;
}


void PathFinder::untangle_pins() {

    update_prev_next();
    uint64_t pins=0;
    for (int e = 0; e < mHBV.EdgeObjectCount(); ++e) {
        if (mToLeft[e]==mToLeft[mInv[e]] and next_edges[e].size()==1 ) {
            if (mVerbose) std::cout<<" Edge "<<e<<" forms a pinhole!!!"<<std::endl;
            if (next_edges[next_edges[e][0]].size()==2) {
                std::vector<uint64_t> pfw = {mInv[next_edges[next_edges[e][0]][0]],mInv[next_edges[e][0]],e,next_edges[e][0],next_edges[next_edges[e][0]][1]};
                std::vector<uint64_t> pbw = {mInv[next_edges[next_edges[e][0]][1]],mInv[next_edges[e][0]],e,next_edges[e][0],next_edges[next_edges[e][0]][0]};
                auto vpfw=multi_path_votes({pfw});
                auto vpbw=multi_path_votes({pbw});
                if (mVerbose) std::cout<<"votes FW: "<<vpfw[0]<<":"<<vpfw[1]<<":"<<vpfw[2]<<"     BW: "<<vpbw[0]<<":"<<vpbw[1]<<":"<<vpbw[2]<<std::endl;
            }
            ++pins;
        }
    }
    if (mVerbose) std::cout<<"Total number of pinholes: "<<pins;
}

void PathFinder::untangle_complex_in_out_choices(uint64_t large_frontier_size, bool verbose_separation) {
    //extend_bridging_paths();
    //find a complex path
    uint64_t qsf=0,qsf_paths=0;
    uint64_t msf=0,msf_paths=0;
    update_prev_next();
    if (mVerbose) std::cout<<"vectors initialised"<<std::endl;
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
                    if (mVerbose) std::cout<<" Single direction frontiers for complex region on edge "<<e<<" IN:"<<path_str(f[0])<<" OUT: "<<path_str(f[1])<<std::endl;
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
                                            if (mVerbose) std::cout<<"found first path from "<<in_e<<" to "<< out_e << path_str(pv)<< std::endl;
                                            first_full_paths.push_back({});
                                            int16_t ei=0;
                                            while (mPaths[inp][ei]!=in_e) ei++;

                                            while (mPaths[inp][ei]!=out_e && ei<mPaths[inp].size()) first_full_paths.back().push_back(mPaths[inp][ei++]);
                                            if (ei>=mPaths[inp].size()) {
                                                if (mVerbose) std::cout<<"reversed path detected!"<<std::endl;
                                                reversed=true;
                                            }
                                            first_full_paths.back().push_back(out_e);
                                            //if (mVerbose) std::cout<<"added!"<<std::endl;
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
                                            if (mVerbose) std::cout<<"found first path from "<<in_e<<" to "<< out_e << path_str(pv)<< std::endl;
                                            first_full_paths.push_back({});
                                            int16_t ei=0;
                                            while (pv[ei]!=in_e) ei++;

                                            while (pv[ei]!=out_e && ei<pv.size()) first_full_paths.back().push_back(pv[ei++]);
                                            if (ei>=pv.size()) {
                                                if (mVerbose) std::cout<<"reversed path detected!"<<std::endl;
                                                reversed=true;
                                            }
                                            first_full_paths.back().push_back(out_e);
                                            //if (mVerbose) std::cout<<"added!"<<std::endl;
                                        }
                                    }
                            if (shared_paths) {
                                out_used[out_i]++;
                                in_used[in_i]++;
                                //if (mVerbose) std::cout << "  Shared paths " << in_e << " --> " << out_e << ": " << shared_paths << std::endl;

                            }
                        }
                    }
                    if ((not reversed) and std::count(in_used.begin(),in_used.end(),1) == in_used.size() and
                            std::count(out_used.begin(),out_used.end(),1) == out_used.size()){
                        if (mVerbose) std::cout<<" REGION COMPLETELY SOLVED BY PATHS!!!"<<std::endl;
                        solved_frontiers.insert(f);
                        for (auto p:first_full_paths) paths_to_separate.push_back(p);
                    } /*else if (std::count(in_used.begin(),in_used.end(),1) == in_used.size()-1 and
                            std::count(in_used.begin(),in_used.end(),0) == 1 and
                            std::count(out_used.begin(),out_used.end(),1) == out_used.size()-1 and
                            std::count(out_used.begin(),out_used.end(),0) == 1){
                        //if (mVerbose) std::cout<<" REGION SOLVED BY PATHS and a jump (not acted on!!!)"<<std::endl;
                        //solved_frontiers.insert(f);
                        qsf++;
                        qsf_paths+=in_used.size();
                        unsigned int in_index=0;
                        while (in_used[in_index]!=0) in_index++;
                        unsigned int out_index=0;
                        while (out_used[out_index]!=0) out_index++;
                        if (mVerbose) std::cout<<"Trying to solve region by reducing last choice to unique path between "<<f[0][in_index]<< " and "<< f[1][out_index]<<std::endl;

                        auto all_paths=AllPathsFromTo({f[0][in_index]},{f[1][out_index]},10);
                        if (all_paths.size()==1){
                            solved_frontiers.insert(f);
                            for (auto p:first_full_paths) paths_to_separate.push_back(p);
                            paths_to_separate.push_back(all_paths[0]);

                            if (mVerbose) std::cout<<"Solved!!!"<<std::endl;
                        }
                        else if (mVerbose) std::cout<<"Not solved, "<<all_paths.size()<<" different paths :("<<std::endl;

                    } else if (std::count(in_used.begin(),in_used.end(),0) == 0 and
                        std::count(out_used.begin(),out_used.end(),0) == 0){
                        msf++;
                        msf_paths+=in_used.size();
                    }*/

                }

            }
        }
    }
    if (mVerbose) std::cout<<"Complex Regions solved by paths: "<<solved_frontiers.size() <<"/"<<seen_frontiers.size()<<" comprising "<<paths_to_separate.size()<<" paths to separate"<< std::endl;
    //if (mVerbose) std::cout<<"Complex Regions quasi-solved by paths (not acted on): "<< qsf <<"/"<<seen_frontiers.size()<<" comprising "<<qsf_paths<<" paths to separate"<< std::endl;
    //if (mVerbose) std::cout<<"Multiple Solution Regions (not acted on): "<< msf <<"/"<<seen_frontiers.size()<<" comprising "<<msf_paths<<" paths to separate"<< std::endl;

    uint64_t sep=0;
    std::map<uint64_t,std::vector<uint64_t>> old_edges_to_new;
    for (auto p:paths_to_separate){

        if (old_edges_to_new.count(p.front()) > 0 or old_edges_to_new.count(p.back()) > 0) {
            if (mVerbose) std::cout<<"WARNING: path starts or ends in an already modified edge, skipping"<<std::endl;
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
    if (mVerbose) std::cout<<" "<<sep<<" paths separated!"<<std::endl;

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
                            if (!seen_edges.count(other_n)) {
                                if (mHBV.EdgeObject(other_n).size() >= large_frontier_size) {
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
                    if (mHBV.EdgeObject(n).size() >= large_frontier_size) {
                        //What about frontiers on both sides?
                        out_frontiers.insert(n);
                        for (auto other_p:prev_edges[n]){
                            if (!seen_edges.count(other_p)) {
                                if (mHBV.EdgeObject(other_p).size() >= large_frontier_size) {
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


    //if (mVerbose) std::cout<<" LOOP: "<< edge_pstr(prev_e)<<" ---> "<<edge_pstr(repeat_e)<<" -> "<<edge_pstr(loop_e)<<" -> "<<edge_pstr(repeat_e)<<" ---> "<<edge_pstr(next_e)<<std::endl;
    auto pvlin=path_votes({prev_e,repeat_e,loop_e,repeat_e,next_e});
    auto pvloop=path_votes({prev_e,repeat_e,loop_e,repeat_e,loop_e,repeat_e,next_e});
    //if (mVerbose) std::cout<<" Votes to p->r->l->r->n: "<<pvlin[0]<<":"<<pvlin[1]<<":"<<pvlin[2]<<std::endl;
    //if (mVerbose) std::cout<<" Votes to p->r->l->r->l->r->n: "<<pvloop[0]<<":"<<pvloop[1]<<":"<<pvloop[2]<<std::endl;

    auto pvcircleline=multi_path_votes({{loop_e,repeat_e,loop_e},{prev_e,repeat_e,next_e}});
    //if (mVerbose) std::cout<<" Multi-votes to circle+line "<<pvcircleline[0]<<":"<<pvcircleline[1]<<":"<<pvcircleline[2]<<std::endl;

    if (pvcircleline[0]>0 or pvloop[2]>0 or (pvcircleline[2]==0 and pvcircleline[1]>pvlin[1] and pvcircleline[1]>pvloop[1])) {
        //if (mVerbose) std::cout<<"   CAN'T be reliably unrolled!"<<std::endl;
        return {};
    }
    if (pvloop[3]==0 and pvloop[0]>pvlin[0]){
        //if (mVerbose) std::cout<<"   LOOP should be traversed as loop at least once!"<<std::endl;
        return {};
    }
    if (pvlin==pvcircleline){
        //if (mVerbose) std::cout<<"   UNDECIDABLE as path support problem, looking at coverages"<<std::endl;
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
            //if (mVerbose) std::cout<<"    Coverage FAIL conditions, CAN'T be reliably unrolled!"<<std::endl;
            return {};
        }
        //if (mVerbose) std::cout<<"    Coverage OK"<<std::endl;
    }
    //if (mVerbose) std::cout<<"   UNROLLABLE"<<std::endl;
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
    //if (mVerbose) std::cout<<std::endl<<"Separating path"<<std::endl;
    std::set<uint64_t> edges_fw;
    std::set<uint64_t> edges_rev;
    for (auto e:p){//TODO: this is far too astringent...
        edges_fw.insert(e);
        edges_rev.insert(mInv[e]);

        if (edges_fw.count(mInv[e]) ||edges_rev.count(e) ){ //if (mVerbose) std::cout<<"PALINDROME edge detected, aborting!!!!"<<std::endl;
            return {};}
    }
    //create two new vertices (for the FW and BW path)
    uint64_t current_vertex_fw=mHBV.N(),current_vertex_rev=mHBV.N()+1;
    mHBV.AddVertices(2);
    //migrate connections (dangerous!!!)
    if (verbose_separation) if (mVerbose) std::cout<<"Migrating edge "<<p[0]<<" To node old: "<<mToRight[p[0]]<<" new: "<<current_vertex_fw<<std::endl;
    mHBV.GiveEdgeNewToVx(p[0],mToRight[p[0]],current_vertex_fw);
    if (verbose_separation) if (mVerbose) std::cout<<"Migrating edge "<<mInv[p[0]]<<" From node old: "<<mToLeft[mInv[p[0]]]<<" new: "<<current_vertex_rev<<std::endl;
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
        if (verbose_separation)  if (mVerbose) std::cout<<"Edge "<<nef<<": copy of "<<p[ei]<<": "<<prev_vertex_fw<<" - "<<current_vertex_fw<<std::endl;
        mToLeft.push_back(prev_vertex_fw);
        mToRight.push_back(current_vertex_fw);
        if (! old_edges_to_new.count(p[ei]))  old_edges_to_new[p[ei]]={};
        old_edges_to_new[p[ei]].push_back(nef);

        auto ner=mHBV.AddEdge(current_vertex_rev,prev_vertex_rev,mHBV.EdgeObject(mInv[p[ei]]));
        if (verbose_separation) if (mVerbose) std::cout<<"Edge "<<ner<<": copy of "<<mInv[p[ei]]<<": "<<current_vertex_rev<<" - "<<prev_vertex_rev<<std::endl;
        mToLeft.push_back(current_vertex_rev);
        mToRight.push_back(prev_vertex_rev);
        if (! old_edges_to_new.count(mInv[p[ei]]))  old_edges_to_new[mInv[p[ei]]]={};
        old_edges_to_new[mInv[p[ei]]].push_back(ner);

        mInv.push_back(ner);
        mInv.push_back(nef);
        mEdgeToPathIds.resize(mEdgeToPathIds.size()+2);
    }
    if (verbose_separation) if (mVerbose) std::cout<<"Migrating edge "<<p[p.size()-1]<<" From node old: "<<mToLeft[p[p.size()-1]]<<" new: "<<current_vertex_fw<<std::endl;
    mHBV.GiveEdgeNewFromVx(p[p.size()-1],mToLeft[p[p.size()-1]],current_vertex_fw);
    if (verbose_separation) if (mVerbose) std::cout<<"Migrating edge "<<mInv[p[p.size()-1]]<<" To node old: "<<mToRight[mInv[p[p.size()-1]]]<<" new: "<<current_vertex_rev<<std::endl;
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
                    if (mVerbose) std::cout<<"Warning, a path could not be updated, truncating it to its first element!!!!"<<std::endl;
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

std::vector<int64_t> PathFinder::get_full_paired_path(int r1id, bool collapse_overlaps) {
    std::vector<int64_t> full_path;
    std::vector<int64_t> r2_path;
    for (auto e:mPaths[r1id]) full_path.emplace_back(e);
    if (not mPaths[r1id+1].empty()) {
        if (not full_path.empty()) {
            for (auto ei = mPaths[r1id + 1].rbegin(); ei != mPaths[r1id + 1].rend(); ++ei)
                r2_path.emplace_back(mInv[*ei]);
            //If this overlaps, the first element on the second path needs to be on the first path too
            size_t r1_ovlstart = 0;
            size_t ovl_size;
            while (r1_ovlstart < full_path.size()) {
                if (full_path[r1_ovlstart] == r2_path[0]) {
                    ovl_size = 0;
                    while (r1_ovlstart + ovl_size < full_path.size()
                           and ovl_size < r2_path.size()
                           and full_path[r1_ovlstart + ovl_size] == r2_path[ovl_size])
                        ++ovl_size;
                    if (r1_ovlstart + ovl_size == full_path.size())break;
                }
                ++r1_ovlstart;
            }
            if (r1_ovlstart < full_path.size()) {
                for (auto r2ei = r2_path.begin() + ovl_size; r2ei < r2_path.end(); ++r2ei)
                    full_path.emplace_back(*r2ei);
            } else {
                full_path.emplace_back(-1);
                for (auto &r2e:r2_path)
                    full_path.emplace_back(r2e);
            }
        }
        else full_path=r2_path;
    }

    return full_path;

}

std::vector<std::pair<uint64_t,std::vector<int64_t>>> PathFinder::collect_paths_for_edge(int e,bool collect_inverse, bool collapse_overlaps) {
    std::vector<std::pair<uint64_t,std::vector<int64_t>>> counted_paths;
    std::vector<std::vector<int64_t>> paths;
    std::cout<<"Collecting paths for edge "<<e<<" (inv="<<mInv[e]<<")"<<std::endl;
    std::set<uint64_t> used_paths;
    //process paths for this edge first
    for (auto &pid:mEdgeToPathIds[e]){
        if (used_paths.count(pid)) continue;
        auto ppid=(pid%2 ? pid-1:pid+1);
        used_paths.insert(pid);
        used_paths.insert(ppid);
        if (pid<ppid){
            paths.emplace_back(get_full_paired_path(pid,collapse_overlaps));
        }
        else{
            paths.emplace_back();
            auto fpp=get_full_paired_path(ppid,collapse_overlaps);
            for (auto ei=fpp.rbegin();ei<fpp.rend();++ei) paths.back().emplace_back((*ei==-1?-1:mInv[*ei]));
        }
    }
    if (collect_inverse) {
        for (auto &pid:mEdgeToPathIds[mInv[e]]){
            if (used_paths.count(pid)) continue;
            auto ppid=(pid%2 ? pid-1:pid+1);
            used_paths.insert(pid);
            used_paths.insert(ppid);
            if (pid<ppid){
                paths.emplace_back();
                auto fpp=get_full_paired_path(pid);
                for (auto ei=fpp.rbegin();ei<fpp.rend();++ei) paths.back().emplace_back((*ei==-1?-1:mInv[*ei]));
            }
            else{
                paths.emplace_back(get_full_paired_path(ppid));
            }
        }
    }
    if (paths.empty()) return {};
    //TODO: sort, count, add to answer vector
    std::sort(paths.begin(),paths.end());
    //uniq count

    for (auto &p:paths){
        if (counted_paths.empty() or counted_paths.back().second!=p){
            counted_paths.emplace_back(std::make_pair(1,p));
        }
        else {
            counted_paths.back().first++;
        }
    }

    return counted_paths;
}

void PathFinder::unroll_loops(int repeat_size, int loop_size, int side_size) {
    for (auto loop_e=1;loop_e<next_edges.size();++loop_e){
        //if (mHBV.EdgeObject(loop_e).size()>loop_size) continue;
        if (prev_edges[loop_e].size()!=1 or
            next_edges[loop_e].size()!=1 or
            prev_edges[loop_e][0]!=next_edges[loop_e][0]) continue;
        auto repeat_e=prev_edges[loop_e][0];
        //if (mHBV.EdgeObject(repeat_e).size()>repeat_size) continue;
        //2) the repeat edge has only one other neighbour on each direction, and it is a different one;
        if (prev_edges[repeat_e].size()!=2 or
            next_edges[repeat_e].size()!=2) continue;

        auto prev_e=(prev_edges[repeat_e][0]==loop_e ? prev_edges[repeat_e][1]:prev_edges[repeat_e][0]);

        auto next_e=(next_edges[repeat_e][0]==loop_e ? next_edges[repeat_e][1]:next_edges[repeat_e][0]);

        if (prev_e==next_e or prev_e==mInv[next_e]) break;

        std::cout<<"loop detected in edge"<<loop_e<<",edge"<<repeat_e<<std::endl;

        //3) size constraints: prev_e and next_e must be at least 1Kbp
        //if (mHBV.EdgeObject(prev_e).size()<min_size_sizes or mHBV.EdgeObject(next_e).size()<min_size_sizes) return {};
        //if (loop_e==1170) std::cout<<"fourth check passed!"<<std::endl;
    }
}

void PathFinder::solve_perfect_repeats(int max_size) {
    /* First part -> find repeats
     * repeat_e is the minimum of the inversion pair
     * size(repeat_e) <= max_size
     * repeat_e has same number of ins and outs >1
     * all ins and outs only connect to repeat_e
     * all ins and outs are different edges (including inversion)
     * every in is connected to an out only
     * every out is connected to only an in
     * WARNING: this may exclude palindromic edges!
     */
    std::vector<std::vector<uint64_t >> sol_paths;
    for (auto repeat_e=1;repeat_e<next_edges.size();++repeat_e){

        if (repeat_e>mInv[repeat_e]) continue;

        if (mHBV.EdgeObject(repeat_e).size()>max_size) continue;

        if (prev_edges[repeat_e].size()<=1 or prev_edges[repeat_e].size()!=next_edges[repeat_e].size()) continue;

        bool too_complex=false;
        for (auto pe:prev_edges[repeat_e]) if (next_edges[pe].size()>1) {too_complex=true;}
        for (auto ne:next_edges[repeat_e]) if (prev_edges[ne].size()>1) {too_complex=true;}
        if (too_complex) continue;

        std::set<uint64_t> seen_edges;
        seen_edges.insert(repeat_e);
        seen_edges.insert(mInv[repeat_e]);
        for (auto pe:prev_edges[repeat_e]) {
            seen_edges.insert(pe);
            seen_edges.insert(mInv[pe]);
        }
        for (auto ne:next_edges[repeat_e]) {
            seen_edges.insert(ne);
            seen_edges.insert(mInv[ne]);
        }
        if (seen_edges.size()!=4*next_edges[repeat_e].size()+2) continue;

        std::cout << std::endl << "Perfect repeat x" << prev_edges[repeat_e].size() << " ,edge" << repeat_e << std::endl;
//        for (auto e:prev_edges[repeat_e]) std::cout<<" "<<e;
//        std::cout<<" -- ";
//        for (auto e:next_edges[repeat_e]) std::cout<<" "<<e;
//        std::cout<<std::endl;

        //TODO: make this generic, rather than hardcoded for CN=2
        if (prev_edges[repeat_e].size()==2){
            auto p1=prev_edges[repeat_e][0];
            auto p2=prev_edges[repeat_e][1];
            auto n1=next_edges[repeat_e][0];
            auto n2=next_edges[repeat_e][1];

            uint64_t p1n1=0,p1n2=0,p2n1=0,p2n2=0;
            for (auto p:collect_paths_for_edge(p1)) {
                if (std::find(p.second.begin(), p.second.end(), n1) != p.second.end()) p1n1+=p.first;
                if (std::find(p.second.begin(), p.second.end(), n2) != p.second.end()) p1n2+=p.first;
            }
            for (auto p:collect_paths_for_edge(p2)) {
                if (std::find(p.second.begin(), p.second.end(), n1) != p.second.end()) p2n1+=p.first;
                if (std::find(p.second.begin(), p.second.end(), n2) != p.second.end()) p2n2+=p.first;
            }
            std::cout<<"   "<<p1n1<<" "<<p1n2<<" "<<p2n1<<" "<<p2n2<<std::endl;
            if (p1n1>0 and p2n2>0 and p1n2==0 and p2n1==0) {
                std::cout<<" SOLVED: "<<p1<<" -> "<<repeat_e<<" -> "<<n1<<" / "<<p2<<" -> "<<repeat_e<<" -> "<<n2<<std::endl;
                sol_paths.push_back({p1,(uint64_t)repeat_e,n1});
            }
            if (p1n2>0 and p2n1>0 and p1n1==0 and p2n2==0) {
                std::cout<<" SOLVED: "<<p1<<" -> "<<repeat_e<<" -> "<<n2<<" / "<<p2<<" -> "<<repeat_e<<" -> "<<n1<<std::endl;
                sol_paths.push_back({p1,(uint64_t)repeat_e,n2});
            }
            //Just separate one, the inverse will be separated the same
        }


        //3) size constraints: prev_e and next_e must be at least 1Kbp
        //if (mHBV.EdgeObject(prev_e).size()<min_size_sizes or mHBV.EdgeObject(next_e).size()<min_size_sizes) return {};
        //if (loop_e==1170) std::cout<<"fourth check passed!"<<std::endl;
    }
    //Separation can be done as in
    std::map<uint64_t,std::vector<uint64_t>> old_edges_to_new;
    for (auto p:sol_paths){
        auto oen=separate_path(p);
        if (oen.size()>0) {
            for (auto et:oen){
                if (old_edges_to_new.count(et.first)==0) old_edges_to_new[et.first]={};
                for (auto ne:et.second) old_edges_to_new[et.first].push_back(ne);
            }
        }
    }
    if (old_edges_to_new.size()>0) {
        migrate_readpaths(old_edges_to_new);
    }
    Cleanup( mHBV, mInv, mPaths );
    invert(mPaths,mEdgeToPathIds,mHBV.EdgeObjectCount());
    update_prev_next();
}

void PathFinder::improve_paths() {
    ReroutePaths(mHBV,mInv,mPaths,mBases,mQuals);
    invert(mPaths,mEdgeToPathIds,mHBV.EdgeObjectCount());
}

void simplifyWithPathFinder( HyperBasevector& hbv, vec<int>& inv, ReadPathVec& paths, VecULongVec& invPaths, const vecbasevector& bases, const VecPQVec& quals, int min_reads, bool verbose, bool dump_intermediate_gfas){

    /*Basically, simplifying in discovar was:
     * - Improve read placements and delete funky pairs ( ReroutePaths / DeleteFunkyPathPairs )
     * - Remove "unsupported edges"
     * - Tip clipping
     * - "Analyze branches"
     * - PopBubbles
     * - Degloop
     * - UnwindThreeEdgePlasmids(hb, inv, paths)
     * - remove small components
     * - RemoveUnneededVerticesGeneralizedLoops
     */

    /* In w2rap, the pathfinder added more loop resolution and Complex Repeat Assembly Paths separation (i.e. directional regions)
     * This worked only to a certain extent. The loop resolution had problems, and complex repeats were only sometimes solved.
     *
     * The new ideas behind this new implementation are:
     * - simplify as in discovar by expanding edges that create lines that are then simplified by RemoveUnnecessaryVertices
     * - Start by improving paths (reasonably) then make an effort not to discard or truncate paths lazily (DONE using DV functions)
     * - Use PE all the way through.
     * - Tip-clipping -> repeats -> shared-path-separation (extend_to_repeat) -> repeats -> loops
     * - It would be important to have clean functions that perform operations clearly and update the inversion of the graph while doing so
     */

    OutputLog(2)<<"======= Simplifying with pathfinder!!!!! ======="<<std::endl;
    auto pf1=PathFinder(hbv,inv,paths,invPaths,bases,quals,min_reads,verbose);
    //First clip stupid tips

    pf1.improve_paths();
    pf1.clip_tips(400,20);

    pf1.improve_paths();
    pf1.clip_tips(800,2);

    pf1.improve_paths();
    pf1.solve_perfect_repeats(1000);
    //RemoveUnneededVertices2(hbv, inv, paths);
    pf1.improve_paths();
    //then unroll loops
    // new heuristic: is the loop in the range of reads going full-through? if yes, check first for full through without collapsing.
    //  if not, collapse and check for looping end? if R is too big, it may be impossible, probably best to leave in peace

    //maybe add a function to find full-path-next and pe-next for an edge?


    //then expand perfectly-solved repeats

    //then extend-to-repeat



    //Now extend-to-repeats
    //pf1.extend_until_repeat();
    //Now unroll loops
    //pf1.mVerbose=true;
    //pf1.unroll_loops(0,0,0);
    /*for (auto cp:pf1.collect_paths_for_edge(1170)){
        std::cout<<cp.first<<"x ";
        for (auto &e:cp.second) std::cout<<e<<" ";
        std::cout<<std::endl;
    }
    for (auto cp:pf1.collect_paths_for_edge(1460)){
        std::cout<<cp.first<<"x ";
        for (auto &e:cp.second) std::cout<<e<<" ";
        std::cout<<std::endl;
    }
    for (auto cp:pf1.collect_paths_for_edge(374)){
        std::cout<<cp.first<<"x ";
        for (auto &e:cp.second) std::cout<<e<<" ";
        std::cout<<std::endl;
    }*/
    ;
    //Now find perfect repeats ?

    //
};