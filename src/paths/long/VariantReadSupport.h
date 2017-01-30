#ifndef VARIANT_READ_SUPPORT_H
#define VARIANT_READ_SUPPORT_H
#include "CoreTools.h"
#include "Basevector.h"
#include "Qualvector.h"
#include "paths/HyperBasevector.h"

// forward declaration
class read_place;
class ReadOriginTracker;
class VariantCallGroup;
class Variant;
class long_logging;

// Given the variants, calculate the probability by reconstructing a graph
// that bubble at EACH variant.
void FindVariantProb(const ReadOriginTracker* p_read_tracker,
        const vec<VariantCallGroup>& vcall_groups,
        const vecbasevector& Gplus, 
        const vec<int>& Gplus_ext, 
        std::map<Variant, vec<std::pair<double,double>>>& probs,
        const long_logging* logc);

void FindReadHomesBest(const vecbasevector& bases, const QualVecVec& quals,
        const HyperBasevector& bubble_graph,
        vec< vec< std::pair<int,int> > >* homes_index,
        vec<int> edges_to_show_supports,
        vec<std::tuple<int,int,read_place>>* edge_rid_place,
        int verbosity = 0,const bool bSafeFindPlaces=false);

#endif
