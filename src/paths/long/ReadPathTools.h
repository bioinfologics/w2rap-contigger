///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// ReadPath Utilities

#ifndef READPATHTOOLS_H_
#define READPATHTOOLS_H_

#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"


// Validates a single read path against an HBV. Returns false if the path is invalid.

bool ValidateReadPath(const HyperBasevector& hbv, const vec<int>& to_left,
		      const vec<int>& to_right, const int offset,
		      const vec<int>& edge_list, String& message,
		      const int read_length = 0);

// Validates all read paths against an HBV. Returns false if any invalid paths are found.

bool ValidateAllReadPaths(const HyperBasevector& hbv, const ReadPathVec& readpaths );



#endif /* READPATHTOOLS_H_ */
