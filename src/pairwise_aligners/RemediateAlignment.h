// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 


#ifndef REMEDIATEALIGNMENT
#define REMEDIATEALIGNMENT

#include "Basevector.h"
#include "PackAlign.h"
#include "Qualvector.h"
#include "system/Types.h"
#include "Vec.h"

Bool RemediateAlignment( const basevector& rd1, const qualvector& q1,
     const basevector& rd2, const qualvector& q2, align& a, int& errors,
     float& score );

#endif
