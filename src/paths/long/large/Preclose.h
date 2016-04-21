///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef PRECLOSE_H
#define PRECLOSE_H

#include "feudal/PQVec.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"

void Preclose( vecbvec const& reads, VecPQVec const& quals,
        String const& work_dir, HyperBasevector& hb, vec<int>& inv,
        ReadPathVec& paths, const Bool preclose_verbose );

#endif
