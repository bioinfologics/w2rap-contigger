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

void PathFinder::untangle_single_choices() {
    //find nodes where in >1 or out>1 and in>0 and out>0
    const int MIN_START_SEQ_SIZE=1000;
    const int MAX_NEXT_SEQ_SIZE=1000;
    uint64_t uloop=0,ursize=0;
    std::cout<<"Starting path finding"<<std::endl;
    init_prev_next_vectors();
    std::cout<<"Prev and Next vectors initialised"<<std::endl;
    //score all possible transitions, discards all decidible and

        // is there any score>0 transition that is not incompatible with any other transitions?

    for ( int e = 0; e < mHBV.EdgeObjectCount(); ++e ) {
        if (e<mInv[e]) {
            auto urs=is_unrollable_loop(e,1000);
            if (urs>0) {
                ++uloop;
                ursize+=urs;
            }
        }

    }
    std::cout<<"Unrollable loops: "<<uloop<<" ("<<ursize<<"bp)"<<std::endl;

    std::cout<<"Path finding finished"<<std::endl;
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

void PathFinder::untangle_complex_in_out_choices() {
    //find a complex path
    init_prev_next_vectors();
    for (int e = 0; e < mHBV.EdgeObjectCount(); ++e) {
        if (e < mInv[e] && mHBV.EdgeObject(e).size() < 1000) {
            std::cout<<" Large frontiers for small edge "<<e<<" "<<path_str(get_all_long_frontiers(e))<<std::endl;
        }
    }
}
std::vector<uint64_t> PathFinder::get_all_long_frontiers(uint64_t e){
    std::set<uint64_t> seen_edges, long_frontiers={e};

    int horizon=10;
    while (horizon--){
        for (auto f:long_frontiers){
            if (mHBV.EdgeObject(f).size()<2000 ) {
                if (!seen_edges.count(f)){
                    for (auto p:prev_edges[f]) if (!seen_edges.count(p)) long_frontiers.insert(p);
                    for (auto n:next_edges[f]) if (!seen_edges.count(n)) long_frontiers.insert(n);
                }
                seen_edges.insert(f);
                long_frontiers.erase(f);
            }
        }
    }

    return std::vector<uint64_t> (long_frontiers.begin(),long_frontiers.end());

}

uint64_t PathFinder::is_unrollable_loop(uint64_t loop_e, uint64_t min_size_sizes){
    //Checks if loop can be unrolled
    //Input: looping edge (i.e. prev_e--->R---->loop_e---->R---->next_e)
    uint64_t prev_e, repeat_e, next_e;
    //Conditions for unrollable loop:

    //1) only one neighbour on each direction, and the same one (repeat_e).
    if (prev_edges[loop_e].size()!=1 or
        next_edges[loop_e].size()!=1 or
        prev_edges[loop_e][0]!=next_edges[loop_e][0]) return 0;
    repeat_e=prev_edges[loop_e][0];


    //2) the repeat edge has only one other neighbour on each direction, and it is a different one;
    if (prev_edges[repeat_e].size()!=2 or
        next_edges[repeat_e].size()!=2) return 0;

    prev_e=(prev_edges[repeat_e][0]==loop_e ? prev_edges[repeat_e][1]:prev_edges[repeat_e][0]);

    next_e=(next_edges[repeat_e][0]==loop_e ? next_edges[repeat_e][1]:next_edges[repeat_e][0]);

    if (prev_e==next_e or prev_e==mInv[next_e]) return 0;

    //3) size constraints: prev_e and next_e must be at least 1Kbp
    if (mHBV.EdgeObject(prev_e).size()<min_size_sizes or mHBV.EdgeObject(next_e).size()<min_size_sizes) return 0;


    //std::cout<<" LOOP: "<< edge_pstr(prev_e)<<" ---> "<<edge_pstr(repeat_e)<<" -> "<<edge_pstr(loop_e)<<" -> "<<edge_pstr(repeat_e)<<" ---> "<<edge_pstr(next_e)<<std::endl;
    auto pvlin=path_votes({prev_e,repeat_e,loop_e,repeat_e,next_e});
    auto pvloop=path_votes({prev_e,repeat_e,loop_e,repeat_e,loop_e,repeat_e,next_e});
    //std::cout<<" Votes to p->r->l->r->n: "<<pvlin[0]<<":"<<pvlin[1]<<":"<<pvlin[2]<<std::endl;
    //std::cout<<" Votes to p->r->l->r->l->r->n: "<<pvloop[0]<<":"<<pvloop[1]<<":"<<pvloop[2]<<std::endl;

    auto pvcircleline=multi_path_votes({{loop_e,repeat_e,loop_e},{prev_e,repeat_e,next_e}});
    //std::cout<<" Multi-votes to circle+line "<<pvcircleline[0]<<":"<<pvcircleline[1]<<":"<<pvcircleline[2]<<std::endl;

    if (pvcircleline[0]>0 or pvloop[2]>0 or (pvcircleline[2]==0 and pvcircleline[1]>pvlin[1] and pvcircleline[1]>pvloop[1])) {
        //std::cout<<"   CAN'T be reliably unrolled!"<<std::endl;
        return 0;
    }
    if (pvloop[3]==0 and pvloop[0]>pvlin[0]){
        //std::cout<<"   LOOP should be traversed as loop at least once!"<<std::endl;
        return 0;
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
            return 0;
        }
        //std::cout<<"    Coverage OK"<<std::endl;
    }
    //std::cout<<"   UNROLLABLE"<<std::endl;
    return mHBV.EdgeObject(prev_e).size()+mHBV.EdgeObject(repeat_e).size()+
            mHBV.EdgeObject(loop_e).size()+mHBV.EdgeObject(repeat_e).size()+mHBV.EdgeObject(next_e).size()-4*mHBV.K();
}

uint64_t PathFinder::paths_per_kbp(uint64_t e){
    return 1000 * mEdgeToPathIds[e].size()/mHBV.EdgeObject(e).size();
};