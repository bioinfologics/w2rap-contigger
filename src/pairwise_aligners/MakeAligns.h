// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 


#ifndef MAKEALIGNS_H
#define MAKEALIGNS_H

#include "system/Assert.h"
#include "Basevector.h"
#include "pairwise_aligners/MakeAlignsMethod.h"
#include "pairwise_aligners/MakeAlignsToCompare.h"
#include "dvString.h"
#include "system/Types.h"
#include "Vec.h"

/*
 * offset_A, offset_B: when the alignment is saved, the ids of the reads
 *   in the alignment are "corrected" by an offset value as follows:
 *   if ( id < N0 ) then id += offset_A; else id += offset_B.
 */



template<int I, int k, int BLOCKS_PER_NODE> void MakeAligns( int Passes, 
     int total_passes, const vecbasevector& EE, 
     to_compare which_to_compare, int max_clique, int max_badness, 
     String aligns_file, ostream& log, int max_errs = 31, int max_alignments = 255, 
     int min_mutmer = 0, int local_max_errs = 31, int stretch = 2, 
     int end_stretch = 2, Bool dump_multiplicities = False, int nstretch = 1, 
     int local_max_errs_done = 0, Bool avoid_promiscuous_kmers = False, 
     int cl = 20, Bool alt_method = False, int bandwidth = 0, 
     vec< vec< vec<int> > >* allowed_offsets = 0,
     int max_offset_discrep = 0, Bool process_frequent_kmers = False,
     String run_dir = "", int offset_A = 0, int offset_B = 0 );

template<int I, int k, int BLOCKS_PER_NODE> void MakeAligns( int Passes, 
     int total_passes, const vecbasevector& EE, to_compare which_to_compare, 
     makealigns_method *method, int max_clique, String aligns_file, ostream& log, 
     int max_alignments = 255, int min_mutmer = 0,
     Bool dump_multiplicities = False, Bool avoid_promiscuous_kmers = False,
     vec< vec< vec<int> > >* allowed_offsets = 0,
     int max_offset_discrep = 0, Bool process_frequent_kmers = False,
     String run_dir = "", int offset_A = 0, int offset_B = 0 );

#endif
