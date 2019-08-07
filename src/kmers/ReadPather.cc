/*
 * \file ReadPather.cc
 * \author tsharpe
 * \date Dec 13, 2011
 *
 * \brief
 */
#include "kmers/ReadPather.h"

unsigned const EdgeID::NULLVAL;
size_t const KDef::MAX_OFFSET;

#include "feudal/SmallVecDefs.h"
template class SmallVec< UnipathEvidence, MempoolAllocator<UnipathEvidence> >;

#include "feudal/OuterVecDefs.h"
template class OuterVec<UnipathEvidenceVec>;
