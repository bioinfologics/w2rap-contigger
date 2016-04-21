///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Utilities for for DisplayNhood and DisplayPlus.

#ifndef DISPLAY_TOOLS_H
#define DISPLAY_TOOLS_H

#include "CoreTools.h"
#include "kmers/KmerRecord.h"
#include "paths/HyperBasevector.h"

void ParseSeeds( const HyperBasevector& hb, const vec<int>& to_right,
     const String& SEEDS, const int RANDOM_SEED, const String& SEEDS_MINUS,
     vec<int>& seeds );

Bool ParseSeeds( const HyperBasevectorX& hb, const vec<int>& inv,
     const vec< triple<kmer<20>,int,int> >& kmers_plus,
     const vec<vec<vec<vec<int>>>>& lines,
     const vec<String>& genome_names, const vec< std::pair<int,ho_interval> >& ambint,
     Bool& ambflag,
     const vec< vec< std::pair<int,int> > >& hits, const String& SEEDS,
     const int RANDOM_SEED, const String& SEEDS_MINUS, vec<int>& seeds,
     const int max_seeds, std::ostream& tout );

void BuildNhood( const HyperBasevector& hb, const vec<int>& to_left,
     const vec<int>& to_right, const vec<int>& seeds, const int DEPTH, 
     vec<Bool>& invisible );

void BuildNhood( const HyperBasevectorX& hb, const vec<int>& seeds, 
     const int DEPTH, vec<Bool>& invisible, const int max_edges );

#endif
