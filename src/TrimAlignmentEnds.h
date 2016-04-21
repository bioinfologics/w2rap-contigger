// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 


#ifndef TRIM_ALIGNMENT_ENDS
#define TRIM_AIIGNMENT_ENDS

#include "Alignment.h"
#include "Basevector.h"
#include "system/Types.h"

void TrimAlignmentEnds( align& a, Bool RC, const basevector &read, 
     const basevector& contig, Bool trim_begin = True, Bool trim_end = True, 
     int par1 = 3, int par2 = 5, int par3 = 2 );

#endif
