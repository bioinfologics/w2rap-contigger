// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research

#ifndef PERFECTALIGNMENT
#define PERFECTALIGNMENT

#include "Basevector.h"
#include "PackAlign.h"
#include "Qualvector.h"
#include "Vec.h"

int PerfectAlignment( const basevector& rd1, 
  	              const qualvector& q1,
		      const basevector& rd2, 
		      const qualvector& q2, 
		      align& al, 
		      float& score,
                      Bool PERF_ALIGN_FAST_NO_MISMATCHES_IN = True,
                      Bool PERF_ALIGN_FAST_NO_MISMATCHES_OUT = True );

#endif // PERFECTALIGNMENT
