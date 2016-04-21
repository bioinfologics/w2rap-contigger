// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology

#ifndef RMR_H
#define RMR_H

#include "Alignment.h"
#include "math/Arith.h"
#include "Basevector.h"
#include "PackAlign.h"

// ReciprocalMatchRate.  For an alignment, form the list of perfectly aligning 
// segment lengths, compute its length-weighted mean, and take the reciprocal.  
// Technical note: we use the perfect segment length plus one, and adjacent 
// mismatches also cause a value of one to be entered in the list.

Float ReciprocalMatchRate( const alignment_plus& ap, const basevector& s1, 
     const basevector& s2 );

Float ReciprocalMatchRate( const align& a, const basevector& s1, 
     const basevector& s2 );

#endif
