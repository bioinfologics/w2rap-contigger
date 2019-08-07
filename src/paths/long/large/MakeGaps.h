#ifndef MAKE_GAPS_H
#define MAKE_GAPS_H

#include "CoreTools.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"

void MakeGaps(HyperBasevector& hb,
               vec<int>& inv,
               vec<vec<vec<vec<int>>>> &lines,
               vec<int> &npairs,
               ReadPathVec& paths,
               VecULongVec& edgeToPathIds,
               const int MIN_LINE,
               const int MIN_LINK_COUNT,
               const String& work_dir,
               const String& FIN,
               const Bool verbose,
               const Bool GAP_CLEANUP );

#endif
