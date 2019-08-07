#ifndef FILL_PAIRS_H
#define FILL_PAIRS_H

#include "Basevector.h"
#include "PairsManager.h"

// FillPairs.  First truncate given reads when the multiplicity of a kmer
// drops below min_freq.  Path these truncated reads, and in cases where the 
// first and last kmer of a pair lie on a single unipath, fill in the sequence
// between them.
//
// This should be rewritten from scratch.

void FillPairs( const vecbvec& bases, const PairsManager& pairs,
                    const int min_freq, vecbvec& filled, bool newMethod);

#endif
