///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "CoreTools.h"
#include "ParallelVecUtilities.h"
#include "paths/long/CreateGenome.h"
#include "paths/long/SupportedHyperBasevector.h"
#include "Map.h"


void SupportedHyperBasevector::FixWeights( const long_logging& logc )
{    double clock = WallClockTime( );

     for ( int i1 = 0; i1 < NPaths( ); i1++ )
     {    const vec<int>& p1 = Path(i1);
          if ( !InvDef( p1[0] ) ) continue;
          vec<int> p2;
          for ( int j = 0; j < p1.isize( ); j++ )
          {     if ( p1[j] < 0 ) p2.push_back( p1[j] );
               else p2.push_back( Inv( p1[j] ) );    }
          p2.ReverseMe( );
          int i2 = BinPosition( Paths( ), p2 );
          if ( i2 < 0 )
          {    std::cout << "\nAttempting to fix weights, found asymmetric path." << std::endl;
               std::cout << "path = " << printSeq( Path(i1) ) << std::endl;
               std::cout << "inv path = " << printSeq(p2) << std::endl << "Abort." << std::endl;
               TracebackThisProcess( );    }
          if ( WeightFw(i1) != WeightRc(i2) )
          {    fix64_6 w = Max( WeightFw(i1), WeightRc(i2) );
               if ( w > WeightRc(i2) ) 
               {    WeightRcMutable(i2) = w;
                    // Temporary if.
                    if ( WeightsFwOrigin( ).size( ) == WeightsFw( ).size( ) )
                         WeightRcOriginMutable(i2) = WeightFwOriginMutable(i1);    }
               else 
               {    WeightFwMutable(i1) = w;    
                    // Temporary if.
                    if ( WeightsFwOrigin( ).size( ) == WeightsFw( ).size( ) )
                         WeightFwOriginMutable(i1) = WeightRcOriginMutable(i2);    }
                         }
          if ( WeightFw(i2) != WeightRc(i1) )
          {    fix64_6 w = Max( WeightFw(i2), WeightRc(i1) );
               if ( w > WeightRc(i1) )
               {    WeightRcMutable(i1) = w;    
                    // Temporary if.
                    if ( WeightsFwOrigin( ).size( ) == WeightsFw( ).size( ) )
                         WeightRcOriginMutable(i1) = WeightFwOriginMutable(i2);    }
               else 
               {    WeightFwMutable(i2) = w;    
                    // Temporary if.
                    if ( WeightsFwOrigin( ).size( ) == WeightsFw( ).size( ) )
                         WeightFwOriginMutable(i2) = WeightRcOriginMutable(i1);    }
                         }    }

     for ( int i1 = 0; i1 < NPairs( ); i1++ )
     {    const vec<int> &p1 = PairLeft(i1), &q1 = PairRight(i1);
          if ( !InvDef( p1[0] ) || !InvDef( q1[0] ) ) continue;
          vec<int> p2, q2;
          for ( int j = 0; j < p1.isize( ); j++ )
          {     if ( p1[j] < 0 ) p2.push_back( p1[j] );
               else p2.push_back( Inv( p1[j] ) );    }
          for ( int j = 0; j < q1.isize( ); j++ )
          {     if ( q1[j] < 0 ) q2.push_back( q1[j] );
               else q2.push_back( Inv( q1[j] ) );    }
          p2.ReverseMe( ), q2.ReverseMe( );
          int i2 = BinPosition( Pairs( ), std::make_pair( q2, p2 ) );
          if ( i2 < 0 )
          {    std::cout << "\nAttempting to fix weights, found asymmetric pair." << std::endl;
               std::cout << "pair = " << printSeq( PairLeft(i1) ) << " ... "
                    << printSeq( PairRight(i1) ) << std::endl;
               std::cout << "inv path = " << printSeq(q2) << " ... " 
                    << printSeq(p2) << std::endl << "Abort." << std::endl;
               TracebackThisProcess( );    }
          if ( PairData(i1) != PairData(i2) )
          {    PairDataMutable(i2) = PairData(i1);    }    }

     REPORT_TIME( clock, "used fixing weights" );    }


