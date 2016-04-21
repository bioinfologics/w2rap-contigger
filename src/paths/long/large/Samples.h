///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef BROAD_SAMPLES_H
#define BROAD_SAMPLES_H

// Test for hardcoded Broad-specific samples.

#include "CoreTools.h"

void Samples( String& species, const String& SAMPLE, String& X, String& EVALUATE,
     String& SELECT_FRAC, int& READS_TO_USE, const String& DATASET, String& BAM,
     vec<String>& subsam_names );

#endif
