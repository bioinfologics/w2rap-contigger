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
