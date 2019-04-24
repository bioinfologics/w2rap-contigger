///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

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

template class SmallVec< int, std::allocator<int> >;
template class OuterVec<IntVec>;

template class SmallVec< unsigned int, std::allocator<unsigned int> >;
template class OuterVec<UIntVec>;

template class SmallVec< unsigned short, std::allocator<unsigned short> >;
template class OuterVec<UShortVec>;

template class SmallVec< long, std::allocator<long> >;
template class OuterVec<LongVec>;

template class SmallVec< unsigned long, std::allocator<unsigned long> >;
template class OuterVec<ULongVec>;
