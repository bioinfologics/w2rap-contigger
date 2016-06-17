///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
//
//  Author: Neil Weisenfeld - Mar 10, 2014 - <crdhelp@broadinstitute.org>
//

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS



#ifndef BUBBLEFREEN50_H_
#define BUBBLEFREEN50_H_
#include "paths/HyperBasevector.h"

class BubbleFreeN50 {
  public:
    explicit BubbleFreeN50( HyperBasevector const& hbv, int min_len = 0 ) {
        vec<int> new_edges( hbv.EdgeObjectCount( ) );
        for ( int e = 0; e < hbv.EdgeObjectCount( ); e++ )
            new_edges[e] = hbv.EdgeLengthKmers(e);

        digraphE<int> size_graph( hbv, new_edges );

        // report original stats
        mPreNEdges = size_graph.Edges().size();
        mPreN50 = N50( size_graph.Edges() );

        // pop bubbles
        vec<int> to_delete;
        to_delete.reserve(size_graph.N());

        #pragma omp for
        for ( int vi = 0; vi < size_graph.N(); ++vi ) {
            if ( size_graph.FromSize(vi) == 2 &&
                    size_graph.From(vi)[0] == size_graph.From(vi)[1] ) {

                int iedge0 = size_graph.EdgeObjectIndexByIndexFrom(vi,1);
                int iedge1 = size_graph.EdgeObjectIndexByIndexFrom(vi,0);
                int sum = size_graph.EdgeObject(iedge0) +
                          size_graph.EdgeObject(iedge1);

                // change edge 0, delete edge 1
                size_graph.EdgeObjectMutable(iedge0) = sum / 2;
                to_delete.push_back( iedge1 );
            }
        }

        UniqueSort(to_delete);
        size_graph.DeleteEdges(to_delete);

        // std::cout << "cleanup..." << std::endl;
        for ( int i = 0; i < size_graph.N(); ++i ) {
            if ( size_graph.From(i).size() == 1 && size_graph.To(i).size() == 1 &&
                    size_graph.From(i)[0] != i ) {
                size_graph.JoinEdges(i, size_graph.EdgeObjectByIndexTo(i,0) +
                                     size_graph.EdgeObjectByIndexFrom(i,0) );
            }
        }

        vec<int> del2;
        for ( int e = 0; e < size_graph.EdgeObjectCount( ); e++ )
            if ( size_graph.EdgeObject(e) + hbv.K( ) - 1 < min_len )
                del2.push_back(e);
        size_graph.DeleteEdges(del2);

        // Clear out edges that have been removed from the graph
        // and vertices with no edges.
        size_graph.RemoveEdgelessVertices( );
        size_graph.RemoveDeadEdgeObjects( );

        // report N50
        mPostNEdges = size_graph.Edges().size();
        mPostN50 = N50( size_graph.Edges() ) + hbv.K( ) - 1;
    }

    int PreN50() {
        return mPreN50;
    }
    int PostN50() {
        return mPostN50;
    }
    size_t PreNEdges() {
        return mPreNEdges;
    }
    size_t PostNedges() {
        return mPostNEdges;
    }

  private:
    int mPreN50;
    int mPostN50;
    size_t mPreNEdges;
    size_t mPostNEdges;
};



#endif /* BUBBLEFREEN50_H_ */
