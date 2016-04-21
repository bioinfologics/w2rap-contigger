///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "paths/Sepdev.h"
#include "graph/DigraphTemplate.h"

template void digraphE<sepdev>::ComponentEdges(vec<vec<int> >&) const;
template void digraphE<sepdev>::Initialize(const vec<vec<int> >&, const vec<vec<int> >&, const vec<sepdev>&, const vec<vec<int> >&, const vec<vec<int> >&, const Bool );
template void digraphE<sepdev>::readBinary(BinaryReader&);
template digraphE<sepdev> digraphE<sepdev>::Subgraph(vec<int> const&) const;
template void digraphE<sepdev>::ToLeft (vec<int>&) const;
template void digraphE<sepdev>::ToRight(vec<int>&) const;

template int digraphE<fsepdev>::AddEdge(const int, const int, const fsepdev&);
template void digraphE<fsepdev>::AddVertices(int);
template void digraphE<fsepdev>::ComponentEdges(vec<vec<int> >&) const;
template digraphE<fsepdev> digraphE<fsepdev>::Subgraph(vec<int> const&) const;
template void digraphE<fsepdev>::ToLeft (vec<int>&) const;
template void digraphE<fsepdev>::ToRight(vec<int>&) const;

template digraphE<Tsepdev<double> >::digraphE();
template void digraphE<Tsepdev<double> >::DeleteEdgeFrom(int, int);
template void digraphE<Tsepdev<double> >::DeleteEdgesAtVertex(int);
template const Tsepdev<double>& digraphE<Tsepdev<double> >::EdgeObject(int) const;
template Tsepdev<double> const& digraphE<Tsepdev<double> >::EdgeObjectByIndexFrom(int, int) const;
template Tsepdev<double>& digraphE<Tsepdev<double> >::EdgeObjectByIndexFromMutable(int, int);
template Tsepdev<double> const& digraphE<Tsepdev<double> >::EdgeObjectByIndexTo(int, int) const;
template int digraphE<Tsepdev<double> >::EdgeObjectCount() const;
template int digraphE<Tsepdev<double> >::EdgeObjectIndexByIndexFrom(int, int) const;
template int digraphE<Tsepdev<double> >::EdgeObjectIndexByIndexTo(int, int) const;
template vec<Tsepdev<double> > digraphE<Tsepdev<double> >::EdgeObjectsBetween(int, int) const;
template Bool digraphE<Tsepdev<double> >::EdgePaths(const int, const int, vec<vec<int> >&, const int, const int, const int) const;
template vec<Tsepdev<double> > const& digraphE<Tsepdev<double> >::Edges() const;
template vec<vec<int> > const& digraphE<Tsepdev<double> >::FromEdgeObj() const;
template vec<int> const& digraphE<Tsepdev<double> >::FromEdgeObj(int) const;
template void digraphE<Tsepdev<double> >::Initialize(const vec<vec<int> >&, const vec<vec<int> >&, const vec<Tsepdev<double> >&, const vec<vec<int> >&, const vec<vec<int> >&, const Bool);
template int digraphE<Tsepdev<double> >::InputFromOutputTo(int, int) const;
template int digraphE<Tsepdev<double> >::InputToOutputFrom(int, int) const;
template vec<vec<int> > const& digraphE<Tsepdev<double> >::ToEdgeObj() const;
template void digraphE<Tsepdev<double> >::readBinary(BinaryReader&);
template void digraphE<Tsepdev<double> >::writeBinary(BinaryWriter&) const;


template digraphE<Tsepdev<int> >::digraphE();
template int digraphE<Tsepdev<int> >::AddEdge(const int, const int, const Tsepdev<int>&);
template void digraphE<Tsepdev<int> >::Clear();
template void digraphE<Tsepdev<int> >::DeleteEdgeFrom(int, int);
template void digraphE<Tsepdev<int> >::DeleteEdges(const vec<int>&);
template void digraphE<Tsepdev<int> >::DeleteEdges(const vec<int>&, const vec<int>&);
template Tsepdev<int>const& digraphE<Tsepdev<int> >::EdgeObjectByIndexFrom(int, int) const;
template Tsepdev<int>& digraphE<Tsepdev<int> >::EdgeObjectByIndexFromMutable(int, int);
template Tsepdev<int> const& digraphE<Tsepdev<int> >::EdgeObjectByIndexTo(int, int) const;
template int digraphE<Tsepdev<int> >::EdgeObjectCount() const;
template int digraphE<Tsepdev<int> >::EdgeObjectIndexByIndexFrom(int, int) const;
template int digraphE<Tsepdev<int> >::EdgeObjectIndexByIndexTo(int, int) const;
template vec<Tsepdev<int> > digraphE<Tsepdev<int> >::EdgeObjectsBetween(int, int) const;
template vec<Tsepdev<int> > const& digraphE<Tsepdev<int> >::Edges() const;
template vec<int> digraphE<Tsepdev<int> >::EdgesBetween( const int, const int) const;
template vec<Tsepdev<int> >& digraphE<Tsepdev<int> >::EdgesMutable();
template vec<vec<int> > const& digraphE<Tsepdev<int> >::FromEdgeObj() const;
template vec<int> const& digraphE<Tsepdev<int> >::FromEdgeObj(int) const;
template int digraphE<Tsepdev<int> >::InputFromOutputTo(int, int) const;
template int digraphE<Tsepdev<int> >::InputToOutputFrom(int, int) const;
template vec<int> digraphE<Tsepdev<int> >::RemoveDeadEdgeObjects();
template vec<vec<int> > const& digraphE<Tsepdev<int> >::ToEdgeObj() const;
template vec<int> const& digraphE<Tsepdev<int> >::ToEdgeObj(int) const;
template void digraphE<Tsepdev<int> >::Used(vec<unsigned char>&) const;
template void digraphE<Tsepdev<int> >::writeBinary(BinaryWriter&) const;
template void digraphE<Tsepdev<int> >::Initialize(const int);
template void digraphE<Tsepdev<int> >::DeleteEdgesAtVertex(const int);
template void digraphE<Tsepdev<int> >::RemoveEdgelessVertices(const vec<int>&);
template void digraphE<Tsepdev<int> >::ReorderVertices(const vec<int>&);

template Bool digraphE<Tsepdev<double>>::EdgePaths( const vec<int>& left, 
     const vec<int>& right, const int v, const int w, vec< vec<int> >& paths, 
     const int max_copies, const int max_paths, const int max_iterations ) const;
