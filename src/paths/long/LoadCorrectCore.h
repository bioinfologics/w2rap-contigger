///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef LOAD_CORRECT_CORE_H
#define LOAD_CORRECT_CORE_H

#include "Basevector.h"
#include "CoreTools.h"
#include "Qualvector.h"
#include "efasta/EfastaTools.h"
#include "paths/long/Logging.h"
#include "paths/long/LongProtoTools.h"
#include "PairsManager.h"

void CapQualityScores( QualVecVec& cquals, const vec<Bool>& done );


// zero all quality scores associated with corrections (for reads already in memory)
void ZeroCorrectedQuals( vecbasevector const& readsFile, vecbvec const& creads,
                            QualVecVec* pQuals );

void CorrectionSuite( vecbasevector& gbases, QualVecVec& gquals, PairsManager& gpairs,
     const long_heuristics& heur,
     //const long_logging& logc, const long_logging_control& log_control,
     vecbasevector& creads, VecEFasta& corrected, vec<int>& cid, 
     vec<pairing_info>& cpartner, const uint NUM_THREADS, const String& EXIT, 
     /*const double clock,*/ bool useOldLRPMethod/*, LongProtoTmpDirManager &tmp_mgr*/ );

void DefinePairingInfo( const PairsManager & gpairs, const vecbasevector& creads,
     const vec<Bool>& to_delete, vec<int>& cid, VecEFasta& corrected,
     vec<pairing_info>& cpartner/*, const long_logging& logc*/ );

#endif
