/*
 * \file Intvector.cc
 * \author tsharpe
 * \date Aug 24, 2009
 *
 * \brief
 */
#include "Intvector.h"
#include "feudal/SmallVecDefs.h"
#include "feudal/OuterVecDefs.h"

template class SmallVec< int, MempoolAllocator<int> >;
template class OuterVec<IntVec>;

template class SmallVec< unsigned int, MempoolAllocator<unsigned int> >;
template class OuterVec<UIntVec>;

template class SmallVec< unsigned short, MempoolAllocator<unsigned short> >;
template class OuterVec<UShortVec>;

template class SmallVec< long, MempoolAllocator<long> >;
template class OuterVec<LongVec>;

template class SmallVec< unsigned long, MempoolAllocator<unsigned long> >;
template class OuterVec<ULongVec>;
