///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef PLACE_READS0_H
#define PLACE_READS0_H

#include "Basevector.h"
#include "CoreTools.h"
#include "Qualvector.h"
#include "paths/long/EvalByReads.h"

void PlaceReads0( const HyperBasevector& hb, const vecbasevector& bases,
		  const vecqualvector& quals, vec< vec<read_place> >& PLACES, bool log = false );

void PlaceSingleRead0( const HyperBasevector& hb, const basevector& bases_single,
		       const qualvector& quals_single, vec<read_place>& PLACES_SINGLE, bool log = false );


#endif
