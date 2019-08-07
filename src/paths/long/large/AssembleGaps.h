#ifndef ASSEMBLE_GAPS_H
#define ASSEMBLE_GAPS_H

#include "Basevector.h"
#include "CoreTools.h"
#include "feudal/PQVec.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "paths/long/large/GapToyTools.h"

void AssembleGaps2( HyperBasevector& hb, vec<int>& inv2, ReadPathVec& paths2, 
     const vecbasevector& bases, VecPQVec const& quals,
     const String& work_dir, std::vector<int>,
     vecbvec& new_stuff, const Bool CYCLIC_SAVE,
     const int A2V, const int MAX_PROX_LEFT,
     const int MAX_PROX_RIGHT, const int MAX_BPATHS, const int pair_sample );

#endif
