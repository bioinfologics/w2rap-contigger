#ifndef REMODEL_GAP_TOOLS_H
#define REMODEL_GAP_TOOLS_H

#include "Basevector.h"
#include "CoreTools.h"
#include "PairsManager.h"
#include "kmers/KmerRecord.h"



template<int K> void MakeKmerLookup0( const vecbasevector& tigs,
     vec< triple<kmer<K>,int,int> >& kmers_plus );

#endif
