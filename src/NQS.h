// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology

// NQS: compare two bases on two sequences, to see if they meet the NQS standard
// and whether they differ.

#ifndef NQS_H
#define NQS_H

#include "Alignment.h"
#include "Basevector.h"
#include "Bitvector.h"
#include "CoreTools.h"
#include "Qualvector.h"

const int NQS_radius = 5;

const int NQS_na = 0;   // not applicable (NQS not satisfied)
const int NQS_diff = 1; // bases different (and NQS satisfied)
const int NQS_same = 2; // bases same (and NQS satisfied)

int NQS( const basevector& rd1, const basevector& rd2, const qualvector& q1, 
     const qualvector& q2, int p1, int p2, int NQS_floor1, int NQS_floor2,
     int NQS_extra );

int NQS0( const basevector& rd1, const basevector& rd2, int p1, int p2, 
     int NQS_extra );

// NQS_vs_perf: same as NQS, but assumes second sequence is perfect.

int NQS_vs_perf( const basevector& rd1, const basevector& rd2, const qualvector& q1,
     int p1, int p2, int NQS_floor1, int NQS_floor2, int NQS_extra );

// CountNQS.  The optional argument restricted_to1 specifies the bases on rd1
// which should be looked at.  Ditto for restricted_to2.  For the third form, it is 
// assumed that quality score checking is already accounted for by restricted_to*.

void CountNQS( 
     /* inputs: */         const alignment& a, 
                           const basevector& rd1, const basevector& rd2, 
                           const qualvector& q1, const qualvector& q2, 
     /* heuristics: */     int NQS_floor1, int NQS_floor2, int NQS_extra,
     /* outputs: */        int& look, int& see );

void CountNQS(
     /* inputs: */         const align& a,
                           const basevector& rd1, const basevector& rd2,
                           const qualvector& q1, const qualvector& q2,
     /* heuristics: */     int NQS_floor1, int NQS_floor2, int NQS_extra,
     /* outputs: */        int& look, int& see,
     /* optional input: */ const bitvector& restricted_to1 = bitvector( ),
                           const bitvector& restricted_to2 = bitvector( ) );

void CountNQS( 
     /* inputs: */         const align& a, 
                           const basevector& rd1, const basevector& rd2, 
     /* heuristics: */     int NQS_extra,
     /* outputs: */        int& look, int& see,
     /* optional input: */ const bitvector& restricted_to1 = bitvector( ),
                           const bitvector& restricted_to2 = bitvector( ) );

// ConfirmBases: given an alignment between b1 and b2, mark bases on b1 as
// confirmed if they meet the following criteria:
// (a) length(b2) >= 400;
// (b) the alignment has score at most 100;
// (c) the base is centered in an 11-base gap-free block of perfect agreement in the
//     alignment, of quality >= 25 at all the bases, and quality >= 30 at the center
//     base, with the entire block situated at least 20 bases from the end of at
//     least one of the reads.

void ConfirmBases( bitvector& confirmed, const alignment_plus& ap, 
     const basevector& b1, const basevector& b2, const qualvector& q1,
     const qualvector& q2 );

void OverlapTransitivity( int id1, int id2, int id3, const alignment_plus& ap12,
     const alignment_plus& ap13, const alignment_plus& ap23,
     const String& source_dir, const vec<int>& idmap, const vecbasevector& bases,
     const vecbitvector& confirmed, const vec<int>& lengths, vec<int>& invalids,
     Bool& all_untrusted, Bool verbose );

#endif
