
///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef LONG_READ_PATCHER_OPTIMIZER_H
#define LONG_READ_PATCHER_OPTIMIZER_H


#include "paths/LongReadTools.h"



// This function optimizes bv_best_p against all of bvs.
// It assumes that the left and right paddings are the same for all bvs
// and that the paddings are >= than the Smith-Waterman radius nb_rad_SW.

void consensus_compute_padded(const BaseVecVec & bvs,
                              BaseVec * bv_best_p,
                              const int nb_pad_left,
                              const int nb_pad_right,
                              const int nb_rad_SW,
                              const unsigned verbosity,
                              const String label = "");


// This function optimizes bv_best_p against all of bvs.
// It will pad all the bvs and bv_best_p on the left and right
// with some hardcoded random sequence equal in size to the
// Smith-Waterman radius nb_rad_SW.

void consensus_compute(const BaseVecVec & bvs,
                       BaseVec * bv_best_p,
                       const int nb_rad_SW,
                       const unsigned verbosity,
                       const String label = "");





void patcher_optimal(const BaseVecVec & unibases,
                     const vec<GapPatcher> & patchers,
                     const size_t ip_best,
                     const int L,
                     const int nb_rad_SW,
                     const int sz_padding_min,
                     const unsigned i_gap,
                     GapPatcher0 * patcher0_opt_p,
                     const unsigned PATCH_VERBOSITY,
                     vec<double> * timers_p);

#endif
