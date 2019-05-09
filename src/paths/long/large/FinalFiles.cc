///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "FastIfstream.h"
#include "FeudalMimic.h"
#include "math/Functions.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
//#include "paths/long/large/DiscoStats.h"
#include "paths/long/large/FinalFiles.h"
#include "paths/long/large/GapToyTools.h"
#include "paths/long/large/Lines.h"

// Build final assembly files, starting from the results of scaffolding.

void FinalFiles(const HyperBasevector &hb, const vec<int> &inv, const
ReadPathVec &paths, const String &work_dir, const String &prefix,
                const int MAX_CELL_PATHS, const int MAX_DEPTH, bool write_edges_file) {
     // Write some assembly files.

     TestInvolution(hb, inv);

     // Make scaffold lines.

     vec<vec<vec<vec<int>>>> linesx;

     FindLines(hb, inv, linesx, MAX_CELL_PATHS, MAX_DEPTH);
     SortLines(linesx, hb, inv);
     DumpLineFiles(linesx, hb, inv, paths, work_dir, prefix, write_edges_file);

}
