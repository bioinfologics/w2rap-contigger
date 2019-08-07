/*
 * \file LongReadsToPaths.h
 * \author tsharpe
 * \date Aug 21, 2012
 *
 * \brief
 */
#ifndef PATHS_LONG_LONGREADSTOPATHS_H_
#define PATHS_LONG_LONGREADSTOPATHS_H_

#include "Basevector.h"
#include "Vec.h"
#include "paths/HyperBasevector.h"
#include "paths/HyperKmerPath.h"

void LongReadsToPaths( vecbvec const& reads,
                            unsigned K, unsigned coverage,
                            HyperBasevector* pHBV,
                            HyperKmerPath* pHKP=nullptr,
                            vecKmerPath* pPaths=nullptr,
                            vecKmerPath* pPathsRC=nullptr,
                            vec<big_tagged_rpint>* pPathsDB=nullptr );

#endif /* PATHS_LONG_LONGREADSTOPATHS_H_ */
