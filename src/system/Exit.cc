/*
 * \file Exit.cc
 * \author tsharpe
 * \date Feb 6, 2009
 *
 * \brief A replacement for the ::exit(int) function.
 */
#include "system/Exit.h"
#include <stdlib.h>

namespace CRD
{

HOOKFUNC gExitHook;

void exit( int status )
{
    if ( gExitHook ) (*gExitHook)(status);

    if ( !status )
        ::exit(0);
    abort();
}

HOOKFUNC installExitHook( HOOKFUNC fHook )
{
    HOOKFUNC oldHook = gExitHook;
    gExitHook = fHook;
    return oldHook;
}

}
