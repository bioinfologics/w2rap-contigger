//
// Created by Bernardo Clavijo (TGAC) on 14/12/2016.
//


#include "GraphImprover.h"
#include "Simplify.h"
#include "GapToyTools.h"
#include <GFADump.h>

void GraphImprover::improve_graph() {
    //first tip clipping

    path_status(mPaths);
    OverlapValidator oval(mHBV,mInv,mPaths);
    oval.compute_overlap_support();

    std::vector<int> paint;
    paint.reserve(mHBV.EdgeObjectCount());
    for (int e:oval.find_perfect_tips(1000,5)) paint.push_back(e);
    std::cout<<Date()<<": "<<paint.size()<<" perfect tips found"<<std::endl;
    GFADumpDetail("ovlpval_perfect_tips_detail",mHBV,mInv,paint);
    mHBV.DeleteEdges(paint);
    Cleanup(mHBV,mInv,mPaths);
    graph_status(mHBV);
    path_status(mPaths);

    //todo: expand canonnical repeats

}

std::set<uint64_t> GraphImprover::expand_cannonical_repeats(uint64_t min_support, uint64_t min_alternative_support) {
    //mark all repetitions that could be expanded (if solved!)
    mHBV.ToRight(mToRight);
    mHBV.ToLeft(mToLeft);
    uint64_t solved=0;
    OverlapValidator oval(mHBV,mInv,mPaths);
    oval.compute_overlap_support();
    auto ssvp=oval.shared_support_vertex_pairs();
    std::cout<<Date()<<": "<<ssvp.size()<<" vertex pairs with shared support found"<<std::endl;
    std::vector<int> paintcr;
    paintcr.reserve(mHBV.EdgeObjectCount());
    for (auto ssp:ssvp) paintcr.push_back(mHBV.EdgeObjectIndexByIndexFrom(ssp.first,0));
    GFADumpDetail("consecutive_support",mHBV,mInv,paintcr);

    paintcr.clear();

    std::vector<std::vector<uint64_t>> canonical_separations;
    for (auto ssp:ssvp){
        //How many "ins into v1"?
        auto ins=mHBV.ToSize(ssp.first);
        //How many "outs" into v2?
        auto outs=mHBV.FromSize(ssp.second);
        //std::cout<<"checking v1="<<ssp.first<<"(in: "<<ins<<") -> v2="<<ssp.second<<"(out:"<<outs<<")"<<std::endl;
        if (ins!=outs) continue;
        //All different?
        uint64_t in_edges[ins];
        uint64_t out_edges[outs];
        for (auto i=0;i<ins;++i) {
            in_edges[i]=mHBV.EdgeObjectIndexByIndexTo(ssp.first,i);
            //std::cout<<"IN: "<<in_edges[i]<<std::endl;
        }
        bool clash=false;
        for (auto i=0;i<outs;++i) {
            out_edges[i]=mHBV.EdgeObjectIndexByIndexFrom(ssp.second,i);
            //std::cout<<"OUT: "<<out_edges[i]<<std::endl;
            for (auto j=0;j<ins;++j) if (in_edges[j]==out_edges[i] or mInv[in_edges[j]]==out_edges[i]) clash=true;
        }
        if (clash) continue;
        //TODO: now check support (check for in)
        uint64_t votes[ins][outs];
        for (auto &vi:votes) for (auto &vio:vi) vio=0;
        auto s1=oval.crosses_and_jumps(ssp.first);
        auto s2=oval.crosses_and_jumps(ssp.second);
        std::vector<uint64_t> shared_support;
        std::set_intersection(s1.begin(),s1.end(),s2.begin(),s2.end(),std::back_inserter(shared_support));
        for (auto spp:shared_support){
            for (auto i=0;i<ins;++i)
                for (auto o=0;o<outs;++o){
                    if (oval.mInformativePairs[spp].jumps_from_to(in_edges[i],out_edges[o])) ++votes[i][o];
                }
        }
        auto r_edge=mHBV.EdgeObjectIndexByIndexFrom(ssp.first,0);
        paintcr.push_back(r_edge);
        clash=false;
        int64_t outdest[ins];
        for (auto &o:outdest) o=-1;
        for (auto i=0;i<ins and not clash;++i)
            for (auto o=0;o<outs;++o){
                if (votes[i][o]){
                    if (outdest[i]!=-1) {
                        clash=true;
                        break;
                    }
                    else{
                        outdest[i]=out_edges[o];
                    }
                }
            }
        for (auto &od:outdest) if (-1==od) clash=true;
        for (auto oi=0;oi<outs;++oi) for (auto oi2=oi+1;oi2<outs;++oi2) if (outdest[oi]==outdest[oi2]) clash=true;
        if (clash) continue;

        //std::cout<<"repeat solved!  ("<<ins<<":"<<outs<<") ";
        for (auto i=0;i<ins;++i) {
            //TODO: is it right to assume the inverse will exist by definition? it should...
            if (in_edges[i]<mInv[outdest[i]]) canonical_separations.push_back({in_edges[i],r_edge,outdest[i]});
            //std::cout<<in_edges[i]<<"->"<<r_edge<<"->"<<outdest[i]<<"  ";
        }
        //std::cout<<std::endl;
        solved++;

    }
    std::cout<<Date()<<": "<<paintcr.size()<<" vertex pairs evaluated as repeats, "<<solved<<" solved"<<std::endl;
    GFADumpDetail("consecutive_support_to_evaluate",mHBV,mInv,paintcr);
    //TODO: separate vertices and re-route paths
    uint64_t sep=0;
    std::map<uint64_t,std::vector<uint64_t>> old_edges_to_new;
    for (auto p:canonical_separations){
        //std::cout<<"Separating "<<p[0]<<" "<<p[1]<<" "<<p[2]<<std::endl;
        if (old_edges_to_new.count(p.front()) > 0 or old_edges_to_new.count(p.back()) > 0) {
            continue;
        }

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
    RemoveUnneededVertices2(mHBV,mInv,mPaths);
    Cleanup(mHBV,mInv,mPaths);
    TestInvolution(mHBV,mInv);
    return std::set<uint64_t>();
}

void GraphImprover::migrate_readpaths(std::map<uint64_t,std::vector<uint64_t>> edgemap){
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
                    //if (mVerbose) std::cout<<"Warning, a path could not be updated, truncating it to its first element!!!!"<<std::endl;
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

std::map<uint64_t, std::vector<uint64_t>> GraphImprover::separate_path(std::vector<uint64_t> p) {

    //separates a path and its inverse, returns a list of new edges.

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
    mHBV.GiveEdgeNewToVx(p[0],mToRight[p[0]],current_vertex_fw);
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
        mToLeft.push_back(prev_vertex_fw);
        mToRight.push_back(current_vertex_fw);
        if (! old_edges_to_new.count(p[ei]))  old_edges_to_new[p[ei]]={};
        old_edges_to_new[p[ei]].push_back(nef);

        auto ner=mHBV.AddEdge(current_vertex_rev,prev_vertex_rev,mHBV.EdgeObject(mInv[p[ei]]));
        mToLeft.push_back(current_vertex_rev);
        mToRight.push_back(prev_vertex_rev);
        if (! old_edges_to_new.count(mInv[p[ei]]))  old_edges_to_new[mInv[p[ei]]]={};
        old_edges_to_new[mInv[p[ei]]].push_back(ner);

        mInv.push_back(ner);
        mInv.push_back(nef);
        mEdgeToPathIds.resize(mEdgeToPathIds.size()+2);
    }
    mHBV.GiveEdgeNewFromVx(p[p.size()-1],mToLeft[p[p.size()-1]],current_vertex_fw);
    mHBV.GiveEdgeNewToVx(mInv[p[p.size()-1]],mToRight[mInv[p[p.size()-1]]],current_vertex_rev);

    //TODO: cleanup new isolated elements and leading-nowhere paths.
    //for (auto ei=1;ei<p.size()-1;++ei) mHBV.DeleteEdges({p[ei]});
    return old_edges_to_new;

}

