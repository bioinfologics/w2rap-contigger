///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef UNSAT_H
#define UNSAT_H

#include "CoreTools.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"

void Unsat( const HyperBasevector& hb, const vec<int>& inv,
            const ReadPathVec& paths, vec< vec< std::pair<int,int> > >& xs,
            const String& work_dir, const int A2V );

#endif
