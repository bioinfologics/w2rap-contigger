#ifndef GAPTOY_SIMPLIFY_H
#define GAPTOY_SIMPLIFY_H

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "CoreTools.h"
#include "feudal/PQVec.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"

void graph_status(const HyperBasevector &hb);

void path_status(const ReadPathVec &paths);

void SimplifyDV(const String &fin_dir, HyperBasevector &hb, vec<int> &inv,
                ReadPathVec &paths, const vecbasevector &bases, const VecPQVec &quals);
void Simplify(const String &fin_dir, HyperBasevector &hb, vec<int> &inv,
              ReadPathVec &paths, const vecbasevector &bases, const VecPQVec &quals,
              const int MAX_SUPP_DEL, const Bool TAMP_EARLY, const int MIN_RATIO2,
              const int MAX_DEL2,
              const Bool ANALYZE_BRANCHES_VERBOSE2, const String &TRACE_SEQ,
              const Bool DEGLOOP, const Bool EXT_FINAL, const int EXT_FINAL_MODE,
              const Bool PULL_APART_VERBOSE, const vec<int> &PULL_APART_TRACE,
              const int DEGLOOP_MODE, const double DEGLOOP_MIN_DIST,
              const Bool IMPROVE_PATHS, const Bool IMPROVE_PATHS_LARGE,
              const Bool FINAL_TINY, const Bool UNWIND3);

void SimplifyEXP(const String &fin_dir, HyperBasevector &hb, vec<int> &inv,
            ReadPathVec &paths, const vecbasevector &bases, const VecPQVec &quals,
            const int MAX_SUPP_DEL, const Bool TAMP_EARLY, const int MIN_RATIO2,
            const int MAX_DEL2,
            const Bool ANALYZE_BRANCHES_VERBOSE2, const String &TRACE_SEQ,
            const Bool DEGLOOP, const Bool EXT_FINAL, const int EXT_FINAL_MODE,
            const Bool PULL_APART_VERBOSE, const vec<int> &PULL_APART_TRACE,
            const int DEGLOOP_MODE, const double DEGLOOP_MIN_DIST,
            const Bool IMPROVE_PATHS, const Bool IMPROVE_PATHS_LARGE,
            const Bool FINAL_TINY, const Bool UNWIND3, const bool RUN_PATHFINDER,
            const bool dump_pf_files, const bool VERBOSE_PATHFINDER);
#endif
