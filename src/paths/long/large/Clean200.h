///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef CLEAN_200_H
#define CLEAN_200_H

#include "CoreTools.h"
#include "feudal/PQVec.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"

void AnalyzeScores( const HyperBasevectorX& hb, const vec<int>& inv, const int e,
     const vec<vec<int>>& scores, vec<int>& to_delete, const int zpass, 
     const int verbosity, const int version );

void Clean200( HyperBasevector& hb, vec<int>& inv, ReadPathVec& paths,
     const vecbasevector& bases, const VecPQVec& quals, const int verbosity,
     const int version, const Bool REMOVE_TINY );

void Clean200x( HyperBasevector& hb, vec<int>& inv, ReadPathVec& paths,
     const vecbasevector& bases, const VecPQVec& quals, const int verbosity,
     const int version, const Bool REMOVE_TINY );

void GetExtensions( const HyperBasevectorX& hb, const int v,
     const int max_exts, vec<vec<int>>& exts, int& depth );

#endif
