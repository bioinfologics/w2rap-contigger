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
     const String& work_dir,  const String& prefix,
     const int MAX_CELL_PATHS, const int MAX_DEPTH,
     bool write_edges_file=true);

#endif
