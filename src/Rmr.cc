///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Alignment.h"
#include "math/Arith.h"
#include "Basevector.h"
#include "CoreTools.h"
#include "math/Functions.h"
#include "PackAlign.h"
#include "Rmr.h"

Float ReciprocalMatchRate( const alignment_plus& ap, const basevector& s1, 
     const basevector& s2 )
{    vec<int> perf;
     perf.clear( );
     int p1, p2, n2 = s2.size( );
     Bool rc2 = ap.Rc2( );
     avector<int> gaps, lengths;
     ap.a.packalign::Unpack( p1, p2, gaps, lengths );
     for ( unsigned int j = 0; j < lengths.length; j++ )
     {    if ( gaps(j) > 0 ) p2 += gaps(j);
          if ( gaps(j) < 0 ) p1 -= gaps(j);
          int run = 0;
          for ( int x = 0; x < lengths(j); x++ )
          {    Bool mismatch;
               if ( !rc2 ) mismatch = ( s1[p1] != s2[p2] );
               else mismatch = ( s1[p1] != 3 - s2[n2-p2-1] );
               if (mismatch)
               {    perf.push_back( run + 1 );
                    run = 0;    }
               else ++run;
               ++p1; 
               ++p2;    }
          perf.push_back( run + 1 );    }
     return Float(1) / WeightedMean(perf);    }

Float ReciprocalMatchRate( const align& a, const basevector& s1, 
     const basevector& s2 )
{    vec<int> perf;
     perf.clear( );
     int n2 = s2.size( );
     int p1 = a.pos1( ), p2 = a.pos2( );
     for ( int j = 0; j < a.Nblocks( ); j++ )
     {    if ( a.Gaps(j) > 0 ) p2 += a.Gaps(j);
          if ( a.Gaps(j) < 0 ) p1 -= a.Gaps(j);
          int run = 0;
          for ( int x = 0; x < a.Lengths(j); x++ )
          {    Bool mismatch = ( s1[p1] != s2[p2] );
               if (mismatch)
               {    perf.push_back( run + 1 );
                    run = 0;    }
               else ++run;
               ++p1; 
               ++p2;    }
          perf.push_back( run + 1 );    }
     return Float(1) / WeightedMean(perf);    }
