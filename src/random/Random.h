#ifndef RANDOM_RANDOM_H
#define RANDOM_RANDOM_H

#include "random/RNGen.h"
#include <stdlib.h>
#include <time.h>

inline long randomx() { return RNGen::random(); }
inline void srandomx( int x ) { RNGen::srandom(x); }
inline int randint( int u ) { return randomx()%u; }

// Although randomx returns a long, its actual range is that of a
// 31-bit int.  Therefore constructs like "randomx() % billion" are
// quite biased.  For this purpose, big_random is much better, as it
// produces a 62-bit random number:
inline long big_random() { return (randomx()<<31) | randomx();  }

#endif // RANDOM_RANDOM_H
