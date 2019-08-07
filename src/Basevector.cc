/*
 * \file Basevector.cc
 * \author tsharpe
 * \date Sep 23, 2009
 *
 * \brief
 */
#include "Basevector.h"

void ReverseComplement( vecbasevector& vbv )
{
    vecbvec::iterator end(vbv.end());
    for ( vecbvec::iterator itr(vbv.begin()); itr != end; ++itr )
        itr->ReverseComplement();
}

#include "feudal/OuterVecDefs.h"
template class OuterVec<BaseVec>;
template class OuterVec<BaseVec,BaseVec::allocator_type>;
template class OuterVec< OuterVec<BaseVec,BaseVec::allocator_type>,
                         MempoolOwner<unsigned char> >;
