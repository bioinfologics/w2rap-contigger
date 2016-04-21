///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef IMPROVE_PATH_H
#define IMPROVE_PATH_H

#include "CoreTools.h"
#include "Qualvector.h"
#include "feudal/PQVec.h"
#include "kmers/KmerRecord.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"

class path_improver {

     public:

     enum path_status { old_better, same, new_better, indet };

     path_improver( )
     {    show_old_better = False;
          show_new_better = False;
          show_same = False;
          show_indet = False;
          print_align = False;    }

     Bool Logging( ) const 
     { return show_old_better || show_new_better || show_same || show_indet; }

     Bool show_old_better;
     Bool show_new_better;
     Bool show_same;
     Bool show_indet;
     Bool print_align;
};

void ImprovePaths( ReadPathVec& paths, const HyperBasevector& hb,
     const vec<int>& inv, const vecbasevector& bases, const VecPQVec& quals,

     // ids: if nonempty, assume that we're passing this subset of the reads;
     // used for logging only.  If you include a read, you must also include
     // its partner!

     const vec<int64_t>& ids, 

     const path_improver& pimp, const Bool IMPROVE_PATHS_LARGE, 
     const Bool BETSYBOB );

#endif
