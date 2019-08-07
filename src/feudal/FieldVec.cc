/*
 * \file FieldVec.cc
 * \author tsharpe
 * \date Mar 22, 2012
 *
 * \brief
 */

#include "feudal/FieldVecDefs.h"
template class FieldVec< 1, MempoolAllocator<unsigned char> >;
template class FieldVec< 2, MempoolAllocator<unsigned char> >;
template class FieldVec< 4, MempoolAllocator<unsigned char> >;
