///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef MAKE_GAPS_H
#define MAKE_GAPS_H

#include "CoreTools.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"

void MakeGaps( HyperBasevector& hb, vec<int>& inv, ReadPathVec& paths,
     VecULongVec& edgeToPathIds, const int MIN_LINE, const int MIN_LINK_COUNT,
     const String& work_dir, const String& FIN, const Bool verbose,
     const Bool GAP_CLEANUP );

#endif
