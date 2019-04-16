///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Various things go here, including all abnormal termination error messages for
// Discovar (not including internal errors).  This allows us to maintain some
// consistency between these messages.

#ifndef DISCOVAR_TOOLS_H
#define DISCOVAR_TOOLS_H

#include <feudal/CharString.h>

namespace DiscovarTools{

    // Exit 0 upon finding that the assembly is empty, or there are no reads.

    void ExitAssemblyEmpty( );
    void ExitPathsEmpty( );
    void ExitNoReads( );
    void ExitNoCorrectedReads( );
    void ExitShortReads( const String& additional_info="" );
}//namespace DiscovarTools

void SetThreads( uint& NUM_THREADS, bool malloc_per_thread_check = true );

#endif
