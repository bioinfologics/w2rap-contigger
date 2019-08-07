// This file defines class "Permutation".  It is ONE-BASED, for backward
// compatibility.

#ifndef PERMUTATION_H
#define PERMUTATION_H

#include "CoreTools.h"

class Permutation : public vec<int> {

     public:

     Permutation(int);   // construct the identity Permutation
     Permutation( ) { }

};

#endif
