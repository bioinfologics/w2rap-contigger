// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
//


#include "layout/MainArrays.h"
#include "pairwise_aligners/SmithWatScore.h"
#include <fstream>
#include <iostream>
#include <iomanip>

//  double TwoDscores[maxlength+1+(maxbandwidth+1)*2][1+(maxbandwidth+1)*2];

ostream& operator<<( ostream &out, SmithWatScore &SWs ) {

    int first_a = max( 0, SWs.beg_a_ - SWs.bandwidth_ ) -1 ;
    int last_a = min( (int) SWs.length_a_, SWs.end_a_ + SWs.bandwidth_ + 1 ) ;

    int first_b = max( 0, SWs.beg_b_ - SWs.bandwidth_ ) -1  ;
    int last_b = min( (int) SWs.length_b_, SWs.end_b_ + SWs.bandwidth_ + 1 ) ;

    last_a = min( last_a + 2, 24 ) ;
    last_b = min( last_b + 2, 24 ) ;

    out << "first_a,last_a = " << first_a << " " << last_a << std::endl ;
    out << "first_b,last_b = " << first_b << " " << last_b << std::endl ;

    for( int i = first_a; i < last_a+2; i++ ) {
        //    out << " SmithWateScore.cc  SWs.firstB( i ) = " <<  SWs.firstB( i ) << " i = " << i << std::endl ;
        for ( int j = first_b; j < SWs.firstB( i ); j++ ) out << "   " ;
        for ( int j = 1; j < SWs.bandwidth_ - i - 1; j++ ) out << "  " ;
        for ( int j = SWs.firstB( i ) -1; j <= min( last_b, SWs.lastB( i ) )+1; j++ ) {
            //      out << " SmithWatScore.cc: i = " << i << " j = " << j << std::endl ;
            if ( SWs.Legal( i, j ) == false ) out << "  SWs.Legal( i, j ) == false   i = "
                                                      << i << " j = " << j << std::endl ;
            out << setprecision(2) << setw(4) << min(999,max(-999,(int)SWs( i, j ))) ;
            //      out << setprecision(2) << setw(4) << ((SWs( i, j ))) ;
            if ( ( ( i == SWs.beg_a_ ) && ( j == SWs.beg_b_ ) ) ||
                    ( ( i == SWs.end_a_ ) && ( j == SWs.end_b_ ) ) ) out <<"#" ;
            else {
                if ( ( j%5 == 0 ) ) {
                    if ( i%5 == 0 ) {
                        if ( ( i%25 == 0 ) && ( j%25 == 0 ) ) out << "@" ;
                        else out << "," ;
                    } else  out << "." ;
                } else {
                    out << " " ;
                }
            }
        }
        out << std::endl ;
    }
    out << " a # is printed at beg_a, beg_b and at end_a, end_b " << std::endl ;
    out << " a , is printed after each element for which both indices are multiples of 5 " << std::endl ;
    out << " a @ is printed after each element for which both indices are multiples of 25 " << std::endl ;
    out << " a . is printed after each element for which j is a multiple of 5 but i is not " << std::endl ;
    out << " -999 means unitialized; 999 means set to MAX_COST " << std::endl ;
    out << " SWs.beg_a_ = " << SWs.beg_a_  << " SWs.end_a_ = " << SWs.end_a_
        << " SWs.beg_b_ = " << SWs.beg_b_  << " SWs.end_b_ = " << SWs.end_b_  << std::endl ;
    return     out << std::endl ;
}

