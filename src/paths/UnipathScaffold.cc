///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "CoreTools.h"
#include "graph/Digraph.h"
#include "paths/UnipathScaffold.h"
#include "graph/FindCells.h"



#include "graph/DigraphTemplate.h"
template digraphE<linklet>::digraphE();
template void digraphE<linklet>::DeleteEdgeFrom(int, int);
template void digraphE<linklet>::DeleteEdges(vec<int> const&);
template void digraphE<linklet>::DeleteEdgesAtVertex(int);
template linklet const& digraphE<linklet>::EdgeObjectByIndexFrom(int, int) const;
template linklet const& digraphE<linklet>::EdgeObjectByIndexTo(int, int) const;
template int digraphE<linklet>::EdgeObjectCount() const;
template int digraphE<linklet>::EdgeObjectIndexByIndexFrom(int, int) const;
template int digraphE<linklet>::EdgeObjectIndexByIndexTo(int, int) const;
template vec<linklet> const& digraphE<linklet>::Edges() const;
template vec<int> digraphE<linklet>::EdgesBetween(const int, const int) const;
template vec<int> digraphE<linklet>::EdgesBetween(const vec<int>&) const;
template void digraphE<linklet>::Initialize(const vec<vec<int> >&, const vec<vec<int> >&, const vec<linklet>&, const vec<vec<int> >&, const vec<vec<int> >&, const Bool);
template int digraphE<linklet>::InputFromOutputTo(int, int) const;
template int digraphE<linklet>::InputToOutputFrom(int, int) const;

template void digraphE<linklet>::PrettyDOT( std::ostream& out, 
     const vec<double>& lengths, const edge_label_info, Bool label_contigs, 
     Bool label_vertices, const vec<int>* componentsToPrint,
     const vec<String> *label_contigs_extra, const vec<int>* verticesToPrint,
     const vec<Bool>* dashed, const vec<Bool>* invisible,
     const vec<String>* edge_color, const vec<int>* pen_widths, const String,
     const double, const double, const double, const double ) const;

template void digraphE<linklet>::Reverse();

template
digraphE<linklet>::edge_label_info::edge_label_info(digraphE<linklet>::edge_label_info::ConstructorBehavior,
        unsigned char, unsigned char, vec<FeudalString<char, std::char_traits<char> > > const*);

template void digraphE<linklet>::ToLeft(vec<int, std::allocator<int> >&) const;
template void digraphE<linklet>::ToRight(vec<int, std::allocator<int> >&) const;
