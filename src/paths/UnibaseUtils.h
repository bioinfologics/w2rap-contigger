#ifndef __INCLUDE_paths_UnibaseUtils_h
#define __INCLUDE_paths_UnibaseUtils_h

#include "Basevector.h"
#include "CommonSemanticTypes.h"
#include "graph/Digraph.h"

/**
   Function: UnibaseInvolution

   For each <unibase>, identify its reverse complement.
 */
void UnibaseInvolution( const vecbasevector& unibases, vec< int >& toRc );
#endif
// __INCLUDE_paths_UnibaseUtils_h
