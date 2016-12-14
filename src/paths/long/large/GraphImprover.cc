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
    oval.analyse_complex_overlaps();
    std::vector<int> paint;
    paint.reserve(mHBV.EdgeObjectCount());
    for (int e:oval.find_perfect_tips(1000,5)) paint.push_back(e);
    std::cout<<Date()<<": "<<paint.size()<<" perfect tips found"<<std::endl;
    //GFADumpDetail("smallk_perfect_tips_detail",mHBV,mInv,paint);
    mHBV.DeleteEdges(paint);
    Cleanup(mHBV,mInv,mPaths);
    graph_status(mHBV);
    path_status(mPaths);
    oval.find_unconnected_neighbours(3);

}
