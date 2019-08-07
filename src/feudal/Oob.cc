/*
 * \file Oob.cc
 * \author tsharpe
 * \date Jul 2, 2009
 *
 * \brief Just a function to call FatalErr on out-of-bounds access.
 */
#include "feudal/Oob.h"
#include "system/System.h"
#include <cstdlib>

void OutOfBoundsReporter::oob( char const* className, size_t idx, size_t siz )
{
    FatalErr(className << " access out of bounds.  Tried to access " << idx << " when there were only " << siz << " elements.");
}
