///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
#include "paths/long/OverlapReads.h"
#include "VecUtilities.h"
#include "FeudalMimic.h"
#include <queue>

// ================================ static methods =============================

// If tail b1[len1-overlap:len1) is the same as head b2[0: overlap)


