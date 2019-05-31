///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"


bool ValidateReadPath(const HyperBasevector& hbv, const vec<int>& to_left,
		      const vec<int>& to_right, const int offset,
		      const vec<int>& edge_list, String& message, 
		      const int read_length = 0) {

    int edge_count = hbv.EdgeObjectCount();
    uint K = hbv.K();

    // Check edge ids are in bounds
    for (auto edge_id : edge_list)
	if (edge_id >= edge_count) {
	    message = "ERROR - Invalid edge ID: " + ToString(edge_id);
	    return false;
	}
    
    // Validate path through graph
    for (uint i = 0 ; i < edge_list.size(); i++) {
	int edge = edge_list[i];
	int right_v = to_right[edge];
	if ((i < edge_list.size() - 1) && (to_left[edge_list[i+1]] != right_v)) {
	    message = "ERROR - no connection between edges : " + ToString(edge) 
		+ " and " + ToString(edge_list[i+1]);
	    return false;
	}
    }
    
    // Validate negative offset
    if (read_length != 0 && offset < 0 && -offset >= static_cast<int>(read_length)) {
	message = "ERROR - negative offset exceeds read length of: " 
	    + ToString(offset);
	return false;
    }
    
    // Validate positive offset
    int first_edge_length = hbv.EdgeLengthBases(edge_list[0]);
    if (offset < 0 && offset >= first_edge_length) {
	message = "ERROR - postive offset exceeds first edge length of: " 
	    + ToString(first_edge_length);
	return false;
    }

    // Check path through graph (only possible if we know the read length)
    if (read_length != 0) {
	int used_bases =  (offset < 0 ? -offset : 0);
	int edge_start_pos = (offset < 0 ? 0 : offset);
	
	for (size_t i = 0 ; i < edge_list.size(); i++) {

	    int remaining_bases = read_length - used_bases;
	    if ( remaining_bases == 0 ) {
		message = "ERROR - path extends beyond read onto edge: " 
		    + ToString(edge_list[i]);
		return false;
	    }

	    int edge_size = hbv.EdgeLengthBases(edge_list[i]);
	    int trimmed_edge_size = edge_size - ( i == edge_list.size() - 1 ? 0 : (K - 1)); 
	    int edge_end_pos = Min(edge_start_pos + remaining_bases, trimmed_edge_size); 
	    
	    if (edge_start_pos >= edge_end_pos) 
		edge_end_pos = edge_start_pos;

	    used_bases += (edge_end_pos - edge_start_pos);
	    
	    if (edge_start_pos == edge_end_pos)
		edge_start_pos = ( K - 1) - (edge_size - edge_start_pos);
	    else 
		edge_start_pos = 0;
	}
    }

    // Add test for ambigious overlap?

    message = "";
    return true;
}

bool ValidateAllReadPaths(const HyperBasevector& hbv, const ReadPathVec& readpaths ) {
  // Compute left and right indices
  vec<int> to_left, to_right;
  hbv.ToLeft(to_left);
  hbv.ToRight(to_right);

  String message;
  bool found_invalid = false;
  const int max_prints = 10;
  int count = 0;
#pragma omp parallel for shared(readpaths,count) private(message) 
  for (size_t i = 0; i < readpaths.size(); ++i) {
    const ReadPath& path = readpaths[i];
    if (!path.empty() 
        && ValidateReadPath(hbv, to_left, to_right, path.getOffset(), 
          vec<int>(path.begin(), path.end()), message ) == false) 
#pragma omp critical 
      {
        found_invalid = true;
        if ( ++count < max_prints ) 
          std::cout << "Path " << i << " = " << path.getOffset() << ":" << printSeq( path ) << "  " << message <<std::endl;
        if (count==max_prints)
          std::cout << "Too many broken paths, not printing anymore" <<std::endl;
      }
  }
  return !found_invalid;
}