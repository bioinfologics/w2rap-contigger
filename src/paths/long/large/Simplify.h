///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef GAPTOY_SIMPLIFY_H
#define GAPTOY_SIMPLIFY_H

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "CoreTools.h"
#include "feudal/PQVec.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"

void Simplify( const String& fin_dir, HyperBasevector& hb, vec<int>& inv, 
     ReadPathVec& paths, const vecbasevector& bases, const VecPQVec& quals,
     const int MAX_SUPP_DEL, const Bool TAMP_EARLY, const int MIN_RATIO2, 
     const int MAX_DEL2, const Bool PLACE_PARTNERS, 
     const Bool ANALYZE_BRANCHES_VERBOSE2, const String& TRACE_SEQ, 
     const Bool DEGLOOP, const Bool EXT_FINAL, const int EXT_FINAL_MODE,
     const Bool PULL_APART_VERBOSE, const vec<int>& PULL_APART_TRACE,
     const int DEGLOOP_MODE, const double DEGLOOP_MIN_DIST, 
     const Bool IMPROVE_PATHS, const Bool IMPROVE_PATHS_LARGE,
     const Bool FINAL_TINY, const Bool UNWIND3 );

#endif
