///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef FINAL_FILES_H
#define FINAL_FILES_H

#include "CoreTools.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "paths/long/large/GapToyTools.h"
#include "paths/long/large/Lines.h"

// Build final assembly files, starting from the results of scaffolding.

void FinalFiles(
     const HyperBasevector& hb, const vec<int>& inv, const ReadPathVec& paths,
     const vec<String>& subsam_names, const vec<int64_t>& subsam_starts,
     const String& work_dir,
     const int MAX_CELL_PATHS, const int MAX_DEPTH,
     const vecbasevector& G);

#endif
