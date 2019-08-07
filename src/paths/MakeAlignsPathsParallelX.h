/* MakeAlignsPathsParallelX: A module container for the function of the same
 * name, which is designed to provide a parallelized back-end to the pathing
 * code ReadsToPathsCoreX.  This module also includes the PathingProcessor class
 * as an aid to parallelization.
 *
 *
 * Josh Burton
 * March 2009
 *
 ******************************************************************************/



#ifndef MAKE_ALIGNS_PATHS_PARALLEL_X_H
#define MAKE_ALIGNS_PATHS_PARALLEL_X_H

#include "Basevector.h" // basevector
#include "Bitvector.h" // bitvector
#include "pairwise_aligners/BalancedMutmerGraph.h" // BMG





// Central function.  
// Parcelize kmers and run the PathingProcessor in parallel.

template<size_t I>
void MakeAlignsPathsParallelX(const size_t       K, 
                              const BaseVecVec & bases, 
                              VecBitVec        * kmer_chosen,
                              BMG<I>           * bmg, 
                              const String     & PARCEL_HEAD,
                              const size_t       NUM_PROCS,
                              const String       CHECKPOINT_HEAD);



#endif
