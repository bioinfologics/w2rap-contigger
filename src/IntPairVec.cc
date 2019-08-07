/*
 * \file IntPairVec.cc
 * \author tsharpe
 * \date Jan 24, 2013
 *
 * \brief
 */

#include "IntPairVec.h"
#include "feudal/SmallVecDefs.h"
#include "feudal/OuterVecDefs.h"

template class SmallVec< IntPair, MempoolAllocator<IntPair> >;
template class OuterVec<IntPairVec>;
