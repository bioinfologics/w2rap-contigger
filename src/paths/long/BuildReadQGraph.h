///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * BuildReadQGraph.h
 *
 *  Created on: Jan 22, 2014
 *      Author: tsharpe
 */

#ifndef PATHS_LONG_BUILDREADQGRAPH_H_
#define PATHS_LONG_BUILDREADQGRAPH_H_

#include "dvString.h"
#include "feudal/ObjectManager.h"
#include "feudal/PQVec.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"

inline void dumpHBV( std::ostream& out, const HyperBasevector& h ) {
    vec<int> to_left, to_right;
    h.ToLeft(to_left), h.ToRight(to_right);
    for ( int e = 0; e < h.EdgeObjectCount(); e++ ) {
        int v = to_left[e], w = to_right[e];
        h.EdgeObject(e).Print(out,
                              ToString(e)+" [vert_"+ToString(v)+"-->vert_"+ToString(w)+"]");
    }
}

void buildReadQGraph( vecbvec const& reads, ObjectManager<VecPQVec>& quals,
                      bool doFillGaps, bool doJoinOverlaps,
                      unsigned minQual, unsigned minFreq,
                      double minFreq2Fract, unsigned maxGapSize,
                      String const& refFasta,
                      bool useNewAligner, bool repathUnpathed,
                      HyperBasevector* pHBV, ReadPathVec* pPaths,
                      bool const VERBOSE = False );

void rePath( HyperBasevector const& hbv,
             vecbvec const& reads, VecPQVec const& quals,
             bool useNewAligner, bool verbose,
             HyperBasevector* pNewHBV, ReadPathVec* pPaths );

#endif /* PATHS_LONG_BUILDREADQGRAPH_H_ */
