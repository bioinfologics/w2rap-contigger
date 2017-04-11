///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef LONG_PROTO_TOOLS_H
#define LONG_PROTO_TOOLS_H

// MakeDepend: library OMP

#include <omp.h>
#include <memory>
#include <unordered_map>

#include "Basevector.h"
#include "CoreTools.h"
#include "IntPairVec.h"
#include "Qualvector.h"
#include "efasta/EfastaTools.h"
#include "kmers/KmerRecord.h"
#include "paths/HyperEfasta.h"
#include "paths/long/Friends.h"
#include "paths/long/Heuristics.h"
#include "paths/long/Logging.h"
#include "paths/long/LongProtoTools.h"
#include "paths/long/PairInfo.h"
#include "paths/long/ultra/ConsensusScoreModel.h"
#include "system/System.h"
#include "paths/long/CreateGenome.h"

// MACROS.

#define REPORT_TIME( CLOCK, MSG ) \
     if (logc.PRINT_TIME_USED) std::cout << TimeSince(CLOCK) << " " << MSG << std::endl;

#define REPORT_TIMEX( CLOCK, MSG, OUT )                              \
     if (logc.PRINT_TIME_USED && omp_get_thread_num( ) == 0 ) \
          OUT << TimeSince( CLOCK, 1.3 ) << " " << MSG << std::endl;

#define DATE_MSG(MSG) std::cout << Date( ) << ": " << MSG << std::endl;

#define FAIL_MSG(MSG) FatalErr(MSG)

// Class ref_loc: record position of something that's aligned to a reference
// sequence id from start to stop, with rc2 = True if it's reversed relative
// to the reference.

class ref_loc {

     public:

     ref_loc( ) : id(-1) { }
     ref_loc( const int id, const int start, const int stop, const Bool rc2 )
          : id(id), start(start), stop(stop), rc2(rc2) { }

     int id;
     int start;
     int stop;
     Bool rc2;

     Bool Defined( ) const { return id >= 0; }

};

class pairing_status {

     public:

     int partner; // or -1 if unpaired
     int lib;     // or -1 if unpaired

};

class simple_map {

     public:

     vec< std::pair<String,int> > x;

     int operator[]( const String& s ) const
     {    for ( int j = 0; j < x.isize( ); j++ )
               if ( x[j].first == s ) return x[j].second;
          return 0;    }

};

class long_logging_control {

     public:

// it's too dangerous to have a constructor that leaves a gazillion pointers
// uninitialized
//     long_logging_control( ) { }

     long_logging_control( const ref_data& ref, vec<ref_loc>* readlocs, 
          const String& OUT_INT_HEAD, const String& VERB )
          : G(&ref.G), G3(&ref.G3), G3plus(&ref.G3plus), Gplus_ext(&ref.Gplus_ext), 
          is_circular(&ref.is_circular), ploidy(&ref.ploidy), LG(ref.LG), 
          Glocs(&ref.Glocs), G3locs(&ref.G3locs), G3pluslocs(&ref.G3pluslocs), 
          readlocs(readlocs), OUT_INT_HEAD(OUT_INT_HEAD)
          {    vec<String> m;
               ParseStringSet( VERB, m );
               for ( int j = 0; j < m.isize( ); j++ )
               {    ForceAssert( m[j].Contains( ":" ) );
                    verb.x.push( m[j].Before( ":" ), 
                         m[j].After( ":" ).Int( ) );    }    }

     const vecbasevector* G;
     const vecbasevector* G3;
     const vecbasevector* G3plus;
     const vec<int>* Gplus_ext;
     const vec<bool>* is_circular;
     const vec<double>* ploidy;
     int LG;
     const VecIntPairVec* Glocs;
     const VecIntPairVec* G3locs;
     const VecIntPairVec* G3pluslocs;
     vec<ref_loc>* readlocs;
     String OUT_INT_HEAD;
     simple_map verb;



};

void ReportPeakMem( const String msg = "" );

int SelectK2( const VecEFasta& corrected, const double K2frac,
     const long_logging& logc, const long_heuristics& heur );



#endif
