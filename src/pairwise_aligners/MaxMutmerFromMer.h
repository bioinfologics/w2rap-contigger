#ifndef MAXMUTMERFROMMER
#define MAXMUTMERFROMMER

#include "Basevector.h"

void MaxMutmerFromMer( int& pos1, int& pos2, int& len, int& errors,
     const basevector& rd1, const basevector& rd2, Bool strict = False );

void MaxMutmerFromMerRev( int& pos1, int& pos2, int& len, int& errors,
     const basevector& rd1, const basevector& rd2, Bool strict = False );

#endif
