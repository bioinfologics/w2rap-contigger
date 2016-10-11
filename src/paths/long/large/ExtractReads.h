///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef EXTRACT_READS_H
#define EXTRACT_READS_H

#include "Basevector.h"
#include "CoreTools.h"
#include "feudal/PQVec.h"

void ExtractReads( String reads,
     const String& work_dir, vecbvec* pReads, VecPQVec* quals );

#endif
