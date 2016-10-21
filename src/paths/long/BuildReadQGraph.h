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

void buildReadQGraph( vecbvec const& reads, VecPQVec const& quals,
                        bool doFillGaps, bool doJoinOverlaps,
                        unsigned minQual, unsigned minFreq,
                        double minFreq2Fract, unsigned maxGapSize,
                        HyperBasevector* pHBV, ReadPathVec* pPaths, int _K, std::string workdir="",
                        std::string tmpdir="", unsigned char disk_batches=0);



#endif /* PATHS_LONG_BUILDREADQGRAPH_H_ */
