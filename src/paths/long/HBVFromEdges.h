/*
 * HBVFromEdges.h
 *
 *  Created on: Dec 13, 2013
 *      Author: tsharpe
 */

#ifndef PATHS_LONG_HBVFROMEDGES_H_
#define PATHS_LONG_HBVFROMEDGES_H_

#include "Basevector.h"
#include "Vec.h"
#include "paths/HyperBasevector.h"
#include "paths/HyperKmerPath.h"

// Takes a bunch of edges, canonicalizes their order, and builds an HBV from
// each edge and its RC (except palindromes).
// The two optional vec<int> args give you the translation from the original
// edge index to the HBV edge object index for the fwd and RC cases.
void buildHBVFromEdges( vecbvec const& edges, unsigned K, HyperBasevector* pHBV,
                             std::vector<int> &pFwdEdgeXlat, std::vector<int> &pRevEdgeXlat );

void buildHKPFromHBV( HyperBasevector const& hbv,
                        std::vector<int> const& fwdEdgeXlat,
                        std::vector<int> const& revEdgeXlat,
                        HyperKmerPath* pHKP );

#endif /* PATHS_LONG_HBVFROMEDGES_H_ */
