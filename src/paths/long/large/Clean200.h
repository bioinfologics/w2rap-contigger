#ifndef CLEAN_200_H
#define CLEAN_200_H

#include "CoreTools.h"
#include "feudal/PQVec.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"

void AnalyzeScores( const HyperBasevectorX& hb, const vec<int>& inv, const int e,
     const vec<vec<int>>& scores, vec<int>& to_delete );

void Clean200x( HyperBasevector& hb, vec<int>& inv, ReadPathVec& paths,
     const vecbasevector& bases, const VecPQVec& quals, const uint min_size);

void GetExtensions( const HyperBasevectorX& hb, const int v,
     const int max_exts, vec<vec<int>>& exts, int& depth );

#endif
