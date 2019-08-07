#include "CoreTools.h"
#include "math/Permutation.h"

Permutation::Permutation(int n) : vec<int>(n)
{    for ( int i = 0; i < n; i++ )
          (*this)[i] = i+1;    }
