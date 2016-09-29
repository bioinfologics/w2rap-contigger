///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef ASSEMBLE_GAPS_H
#define ASSEMBLE_GAPS_H

#include "Basevector.h"
#include "CoreTools.h"
#include "feudal/PQVec.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "paths/long/large/GapToyTools.h"

void AssembleGaps2( HyperBasevector& hb, vec<int>& inv2, ReadPathVec& paths2, 
     VecULongVec& paths2_index, const vecbasevector& bases, VecPQVec const& quals,
     const String& work_dir, int K2_FLOOR,
     vecbvec& new_stuff, const Bool CYCLIC_SAVE,
     const int A2V, const int GAP_CAP, const int MAX_PROX_LEFT,
     const int MAX_PROX_RIGHT, const int MAX_BPATHS, const int pair_sample );

#endif
