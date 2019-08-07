/// \file
/// This file defines the typedef "qualvector", which stores quality scores
/// between 0 and 255 as a vector of unsigned chars, and the typedef
/// "vecqualvector", which stores a vector of qualvectors.
/// \ingroup grp_quals

#ifndef QUALVECTOR
#define QUALVECTOR

#include "feudal/SerfVec.h"
#include "feudal/MasterVec.h"
#include "Charvector.h"
#include "dvString.h"
#include "Vec.h"
#include <algorithm>
#include <ostream>

/// Logical type for quality scores
typedef uint8_t qual_t;

/// Vector of quality scores, for example representing the quality of each base
/// in one read.
//TODO: BJ modify the vector definition
typedef std::vector<qual_t> QualVec;


/// Vector of vectors of quality scores, for example representing the quality
/// of each base in each read in a set of reads.
//TODO: BJ modify the vector definition
typedef std::vector<std::vector<qual_t>> QualVecVec;

#endif
