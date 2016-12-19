//
// Created by Bernardo Clavijo (TGAC) on 14/12/2016.
//


#include "GraphImprover.h"
#include "Simplify.h"
#include "GapToyTools.h"
#include <GFADump.h>

void GraphImprover::improve_graph() {
    path_status(mPaths);
    OverlapValidator oval(mHBV,mInv,mPaths);
    oval.compute_overlap_support();
    //yoval.analyse_complex_overlaps();




    std::vector<int> paint;
    paint.reserve(mHBV.EdgeObjectCount());
    for (auto &tc:oval.find_unconnected_neighbours(1)) {
        int e=tc.e1;
        paint.push_back(e);
        e=tc.e2;
        paint.push_back(e);
    }
    GFADumpDetail("smallk_new_connections_detail",mHBV,mInv,paint);

    for (int e:oval.find_perfect_tips(1000,5)) paint.push_back(e);
    std::cout<<Date()<<": "<<paint.size()<<" perfect tips found"<<std::endl;

    mHBV.DeleteEdges(paint);
    Cleanup(mHBV,mInv,mPaths);
    graph_status(mHBV);
    path_status(mPaths);
    oval.compute_overlap_support();
    //oval.find_unconnected_neighbours(10);

}

std::set<uint64_t> GraphImprover::expand_cannonical_repeats(uint64_t min_support, uint64_t min_alternative_support) {
    //mark all repetitions that could be expanded (if solved!)
    OverlapValidator oval(mHBV,mInv,mPaths);
    oval.compute_overlap_support();
    auto ssvp=oval.shared_support_vertex_pairs();
    std::cout<<Date()<<": "<<ssvp.size()<<" vertex pairs with shared support found"<<std::endl;
    std::vector<int> paintcr;
    paintcr.reserve(mHBV.EdgeObjectCount());
    for (auto ssp:ssvp) paintcr.push_back(mHBV.EdgeObjectIndexByIndexFrom(ssp.first,0));
    GFADumpDetail("consecutive_support",mHBV,mInv,paintcr);

    paintcr.clear();
    for (auto ssp:ssvp){
        //How many "ins into v1"?
        auto ins=mHBV.ToSize(ssp.first);
        //How many "outs" into v2?
        auto outs=mHBV.FromSize(ssp.second);
        //std::cout<<"checking v1="<<ssp.first<<"(in: "<<ins<<") -> v2="<<ssp.second<<"(out:"<<outs<<")"<<std::endl;
        if (ins!=outs) continue;
        //All different?
        std::vector<uint64_t> in_edges(ins);
        std::vector<uint64_t> out_edges(outs);
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
        //TODO: now check support
        paintcr.push_back(mHBV.EdgeObjectIndexByIndexFrom(ssp.first,0));


    }
    std::cout<<Date()<<": "<<paintcr.size()<<" vertex pairs to evaluate as repeats"<<std::endl;
    GFADumpDetail("consecutive_support_to_evaluate",mHBV,mInv,paintcr);
    //TODO: separate vertices and re-route paths

    return std::set<uint64_t>();
}
