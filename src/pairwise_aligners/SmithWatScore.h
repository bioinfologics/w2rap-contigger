// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
//


//   Header: /seq/wga/ArachneCVS/Arachne/SmithWatScore.h,v 1.13 2000/07/13 13:40:22 kstanley Exp
//
// class SmithWatScore.h
//
//   It would nice to make this into a template, so that we could handle ints
//   doubles and even enums with the same template.
//
#include "system/Assert.h"
#include "math/Functions.h"
#include "math/Arith.h"

const int maxlength = 10000 ;
const int maxbandwidth = 100 ;
//const int lda_ = 1+(maxbandwidth+1)*2 ;
const int lda_ = maxlength+1+(maxbandwidth+1)*2;



//
//  Invariants:
//    beg_a and end_a are inclusive
//    first_a and last_a are inclusive
//    first_a-1 and last_a-1 are legal and set to MAX_COST
//
class SmithWatScore {

    //
    //  I hope to get rid of maxlength and maxbandwidth in the near future,
    //  automatically allocating an array of length x 1+bandwidth*2
    //
  private:
    int length_ ;
    int bandwidth_ ;
    int beg_a_ ;
    int beg_b_ ;
    int end_a_ ;
    int end_b_ ;
    int length_a_ ;
    int length_b_ ;
    int array_refs_;

#define ALLOW_ANY_LENGTH_OR_BANDWIDTH
#ifdef  ALLOW_ANY_LENGTH_OR_BANDWIDTH
    vector < vector<Float> > TwoDscores ;
#else
    Float TwoDscores[maxlength+1+(maxbandwidth+1)*2][2*(maxbandwidth+1)+1];
#endif

  public:
    SmithWatScore( int length, int bandwidth, int beg_a, int end_a, int length_a,
                   int beg_b, int end_b, int length_b ) {
        length_ = max( end_a - beg_a , end_b - beg_b ) + 1 ;

#ifndef ALLOW_ANY_LENGTH_OR_BANDWIDTH
        if ( bandwidth > maxbandwidth ) std::cout << " bandwidth = " << bandwidth << std::endl ;
        Assert( length_ <= maxlength ) ;
        Assert( bandwidth <= maxbandwidth ) ;
#endif
        bandwidth_ = bandwidth ;
        beg_a_ = beg_a ;
        end_a_ = end_a ;
        length_a_ = length_a ;
        beg_b_ = beg_b ;
        end_b_ = end_b ;
        length_b_ = length_b ;
        int ArrayHeight = length + 1 + ( bandwidth + 1 ) * 2 ;
        int ArrayWidth = 1 + ( bandwidth + 1 ) * 2 ;
        array_refs_ = 0 ;

#ifdef ALLOW_ANY_LENGTH_OR_BANDWIDTH
        TwoDscores.resize( ArrayHeight );
        for ( int i = 0 ; i < ArrayHeight ; i ++ )   {
            TwoDscores[i].resize(ArrayWidth);
        }
#endif

#if 0
        //  The following is very inefficient, but it gets the job done for debug
        for( int i = 0 ; i < maxlength+1+maxbandwidth*2 ; i++ ) {
            for( int j = 0; j < 1+maxbandwidth*2; j++ ) TwoDscores[i][j] = -87654.321 ;
        }
#endif

    }
    ~SmithWatScore( ) {}

    // for( j = firstA(i); i <= lastA(i); i++ )
    //   ... Set SmithWatScore( i, j )
    //   Note that SmithWatScore( i, j ) depends on SmithWatScore( i-1, j-1 )
    //   Hence we must allow room for those when computing firstA().

    //  SmithWat( firstA( b ), b ) to  SmithWat( lastA( b ), b ) is the legal range
    //  for this banded SmithWaterman

    int firstA( int b ) {
        return std::max( {b + beg_a_ - beg_b_ - bandwidth_ ,
                          beg_a_ - bandwidth_, 0
                         } ) ;
    }
    int lastA( int b ) {
        return std::min( {b + beg_a_ - beg_b_ + bandwidth_ ,
                          length_ + beg_a_ + bandwidth_ , length_a_ - 1
                         } ) ;
    }
    int firstB( int a ) {
        return std::max( {a + beg_b_ - beg_a_ - bandwidth_ + 1,
                          beg_b_ - bandwidth_, 0
                         } ) ;
    }
    int lastB( int a ) {
        return std::min( {a + beg_b_ - beg_a_ + bandwidth_ ,
                          length_ + beg_b_ + bandwidth_ , length_b_ - 1
                         } ) ;
    }

    inline Float& operator() ( int x, int y ) {
        int truex = x - beg_a_ + bandwidth_ ;
        int truey = y - beg_b_ + bandwidth_ ;
#ifdef CHECK
        if( ( ( truex >= -1 ) && ( truey >= -1 ) && ( Abs( truex - truey ) <= bandwidth_ + 1 )
                && ( truex < length_ + 2 * bandwidth_ ) ) == false ) {
            cout << " SmithWatScore.h truex = " << truex << " truey = " << truey << " bandwidth = " << bandwidth_ << " x = " << x << " y = " << y << std::endl ;
            Assert( false ) ;
        }
        Assert( truex >= -1 );
        Assert( truey >= -1 );
        Assert( Abs( truex - truey ) <= bandwidth_ + 1 ) ;

        Assert( truex < length_ + 2 * bandwidth_ ) ;
#endif
        //    array_refs_++ ;
        return TwoDscores[ truex+1][ ( truey-truex+bandwidth_ ) ];
    } ;

    int ArrayRefs() {
        return array_refs_ ;
    } ;
    int Length() {
        return length_ ;
    } ;
    int Bandwidth() {
        return bandwidth_ ;
    } ;

    //  This is just for debugging
    bool Legal( int x, int y ) {
        int truex = x - beg_a_ + bandwidth_ ;
        int truey = y - beg_b_ + bandwidth_ ;

        //    Assert( false ) ;
        return( ( truex >= -1 ) && ( truey >= -1 ) && ( Abs( truex - truey ) <= bandwidth_ + 1)
                && ( truex < length_ + 2 * bandwidth_ ) ) ;
    } ;

    friend ostream& operator<<( ostream &out, SmithWatScore &SWs );


} ;
