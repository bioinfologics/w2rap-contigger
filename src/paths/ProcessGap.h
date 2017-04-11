///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef PATHS__PROCESS_GAP_H
#define PATHS__PROCESS_GAP_H

#include "Basevector.h"
#include "PairsManager.h"
#include "efasta/EfastaTools.h"
#include "paths/LongReadTools.h"
//#include "paths/Ulink.h"
#include "paths/UnipathScaffold.h"
#include "paths/Uniseq.h"

void GetWalks( const int u1, const int u2, const int sep, const int dev,
     const vecbasevector& unibases, const int K, const vec<int>& to_rc, 
     const vec< vec< std::pair<int,int> > >& nextsx, const vec<int>& use,
     vec< vec< std::pair<int,int> > >& walks1, int& bad, const double max_dev_diff );


#endif
