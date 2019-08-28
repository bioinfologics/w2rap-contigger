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

     path_improver(bool show_old=True, bool show_same=True, bool show_new=True, bool show_indet=True) :
     show_old_better(show_old), show_same(show_same), show_new_better(show_new), show_indet(show_indet) {}

     Bool Logging( ) const 
     { return show_old_better || show_new_better || show_same || show_indet; }

     Bool show_old_better;
     Bool show_new_better;
     Bool show_same;
     Bool show_indet;
     Bool print_align = False;
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
