///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef REMODEL_GAP_TOOLS_H
#define REMODEL_GAP_TOOLS_H

#include "Basevector.h"
#include "CoreTools.h"
#include "PairsManager.h"
#include "kmers/KmerRecord.h"



template<int K> void MakeKmerLookup0( const vecbasevector& tigs,
     vec< triple<kmer<K>,int,int> >& kmers_plus );

#endif
