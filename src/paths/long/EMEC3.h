///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
//


#ifndef PATHS_LONG_EMEC3_H_
#define PATHS_LONG_EMEC3_H_

#include "feudal/BaseVec.h"
#include "simulation/ReadSimulatorSimpleCore.h"
#include "simulation/ReferenceIterator.h"
#include "Qualvector.h"
#include "dna/Bases.h"
#include "paths/long/ReadStack.h"

typedef vec<double > BaseProb;		// vector of base probabilities (always 4 elements)
typedef StackBaseVec _BaseVec;
typedef StackQualVec _QualVec;
typedef StackBaseVecVec _BaseVecVec;
typedef StackQualVecVec _QualVecVec;

/*
 * run_EMEC3 - main entrypoint
 *
 * call  - stack comprising the "founder" read (always position 0) and friends
 * callq - quality scores for reads in call
 * t, q  - output "true" read and updated quality score (zeros for edited positions)
 * pfriend - vector of friendship probabilities
 * trim_to - currently unused
 * debug_level - currently >0 for debugging
 *
 */
void run_EMEC3( const _BaseVecVec& call, const _QualVecVec& callq, _BaseVec& t, _QualVec& q, std::vector<double>& pfriend, int& trim_to, const unsigned debug_level = 0, const size_t serialno = 0 );

#endif /* PATHS_LONG_EMEC3_H_ */
