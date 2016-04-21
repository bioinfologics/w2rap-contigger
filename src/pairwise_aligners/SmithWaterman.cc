// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research


//  #define NDEBUG      // This seems to increase execution time by 10%
//
//  See SmithWaterman.doc for Ken's Feb 20, 2001 thoughts on how to speed this up for
//  use in ReadsToAligns
//
//  This is the common calling sequence:
//
//  bool SmithWaterman( bandwidth, seq1, seq2, qual1, qual2, align &al )
//    Performs a SmithWaterman returning true if a better alignment is found.
//    al is replaced with the better alignment, if any.
//
//  The following calling sequence is used in TestSW and MakeTemplates:
//
//  bool SmithWaterman( bandwidth, seq1, seq2, qual1, qual2, beg_a, end_a, beg_b, end_b, read_template &rt )
//    Performs a SmithWaterman returning a read_template and modifying beg_a,
//    end_a, beg_b and end_b.  (This callinq sequence always returns false).
//
//  The following calling sequences are no longer used:
//
//  bool SmithWaterman( bandwidth, seq1, seq2, qual1, qual2, beg_a, end_a, beg_b, end_b,
//    read_template &rt, align )
//    Performs a SmithWaterman returning a read_template and an alignment,
//    while modifying beg_a, end_a, beg_b and end_b.  Returns true if a better alignment is
//    found.
//
//  bool SmithWaterman( bandwidth, seq1, seq2, qual1, qual2, beg_a, end_a, beg_b, end_b )
//    Performs a SmithWaterman modifying beg_a,
//    end_a, beg_b and end_b.  This always returns false.
//
//  RELEASE NOTES (Feb 20, 2001 )
//    Known bugs:
//      0) The cost of allowing arbitrary length and bandwidth is
//         a significant (say 30%) increase in execution time.
//         At present, I have left in the arbitrary length and
//         bandwidth option but would not be surprised if
//         it has been taken out by the time you read this note.
//         You can switch this yourself by commenting (or uncommenting)
//         ALLOW_ANY_LENGTH_OR_BANDWIDTH in SmithWatScore.h
//      1) The scoring system needs work, specifically:
//         A)  The current scoring system occasionally causes an alignment that
//         begins or ends with a deletion rather than a mutation.
//         B)  The scoring system differs from the system used by ScoreAlignment
//         hence optimal for one does not mean optimal for the other.
//         C)  The scoring system can lead to an insertion followed immediately by a
//         deletion.  The right fix is to fix the scoring system.  However, a simple
//         peephole optimizer could fix this after the fact by looking for inserionts
//         followed immediately by a deletion and replace them both by a mutation.
//         D)  Pow is expensive
//         E)  On the other hand, eliminating pow from the scoring system does not
//         not the quality significantly.  Hence, the scoring sytem does not
//         seem to be critical.
//      2)  It is grossly inefficient as presently used.  As called by ReadsToAligns,
//      it treats stretches of 400 matching bases the same as stretches than contain
//      errors.  We have never seen an improvement that involved bases more than
//      a few bases away from known errors.  Hence, one way to improve performance would
//      be to only run SmithWaterman only on stretches that include errors (perhaps allowing
//      say 10 good bases past the last error.)
//      3)  Testing is woefully inadequate and there is no automated testing.
//      4)  The limitation of alignment to short vectors appears to be a problem for large alignments.
//          It could be that alignments with more than 255 gaps are not interesting, but we are running
//          into alignments with more than 255 gaps, even in h_infl simulation.

// Please please do not include MainArrays.h.  It is unmodular and
// causes problems.  I claim that I have removed all of the unmodular aspects of
// MainArrays.h - what is causing problems?  Ken
#include "layout/MainArrays.h"
#include "PrintAlignment.h"

#include "math/Arith.h"


#include "pairwise_aligners/SmithWatScore.h"
#include "pairwise_aligners/SmithWaterman.h"
#include <cmath>
#include "layout/common.h"
#include "system/System.h"
#include "ScoreAlignment.h"
#include "Qualvector.h"

const Float MAX_COST = 1234567E0 ;

// This is correct _only_ for 8-byte IEEE doubles: // NONSENSE NOW!
#if !defined(DBL_MAX_10_EXP)
#define DBL_MAX_10_EXP		10
#endif

//
//  These should be an enum, but that reuires turning class SmithWatScore into a template.
//  Match, Mutation, Deletion and Insertion are the only valid entries in SWptr().
//
const int Match = 14 ;
const int Mutation = 15 ;
const int Deletion = 16 ;   // Deletion means an element in seq1 is not in seq2
const int Insertion = 17 ;

const Float MatchValue = -1.0 ;
//
//  The following functions compute SmithWaterman costs based on the
//  underlying Phred scores.
//
//  AdjustedScore first takes makes sure that no base has a score more
//  than twice either of its neighbors scores, then it exponentiates the
//  score.
//
//  MutateCost takes the minumum of the scores of the two bases involved.
//  We assume that the scores passed to MutateCost() and InDelCost have
//  already been adjusted (i.e. passed to AdjustedScore).
//
//  InDelCost takes the minimum of the score of the base being deleted or
//  of either of the two bases adjacent between which a base is being
//  inserted.
//
//
//

template <class T> inline T AdjustedScore(const T& x1, const T& x2, const T& x3)
{
  // XXX: This extra Min is to prevent an overflow.  There are better ways to
  // make sure this doesn't happen:
  return  Pow10( Min(T(DBL_MAX_10_EXP), Min( x2, Float(2) * Min(x1, x3))) / Float(25) );
  // return  Min( x2, 2 * Min( x1, x3 ) ) ; // This seems to work just about as well in practice.
}

template <class T> inline T MutateCost(const T& x, const T& y)
{ return Min( x, y ) ; }

template <class T> inline T InDelCost(const T& x, const T& y1, const T& y2)
{ return Float(2) * Min( x, Min ( y1, y2 ) ) ; }

  //
  //  The rest of this code should not have to change.  Changes made to AdjustedScore(),
  //  MutateCost() and InDelCost() will impact the SmithWaterman code.
  //
  //  AdjustedA and AdjustedB are shifted up by 1 because I don't know how to allow
  //  Adjusted[-1]
  //
  //  i.e. AdjustedA[i+1] is based on qual1[i]
  //  i.e. AdjustedB[i+1] is based on qual2[i]
  //

bool InternalSmithWaterman( const int bandwidth,
			    const basevector &seq1,  // Basevectors store bases as 4-bit entities
			    const basevector &seq2,
			    const qualvector &qual1,  // Quality Scores for seq1
			    const qualvector &qual2,
			    int &beg_a, //  The start of the match (within seq1) MODIFIED
			    int &end_a,
			    int &beg_b,
			    int &end_b,
			    bool create_template,
			    read_template &read_temp,
			    bool create_alignment,
			    align &new_alignment) {
  Assert( qual1.size() == seq1.size()) ;
  Assert( qual2.size() == seq2.size()) ;

  int Alength = seq1.size() ;
  int Blength = seq2.size() ;
  Assert( Alength >= 2 ) ;
  Assert( Blength >= 2 ) ;
  Assert( beg_a >= 0 ) ;
  Assert( beg_a <= end_a ) ;
  Assert( end_a <= Alength-1 ) ;
  Assert( beg_b >= 0 ) ;
  Assert( beg_b <= end_b ) ;
  Assert( end_b <= Blength-1 ) ;
  Assert( beg_a == 0 || beg_b == 0 ) ;
  Assert( end_a == Alength-1 || end_b == Blength-1 ) ;

  int length = min( end_a - beg_a, end_b - beg_b ) + 2 * bandwidth ;
  SmithWatScore SWscore( length, bandwidth, beg_a, end_a, seq1.size(), beg_b, end_b,
			 seq2.size() ) ;
  SmithWatScore SWptr( length, bandwidth, beg_a, end_a, seq1.size(), beg_b, end_b,
		       seq2.size() ) ;
  vec< Float > AdjustedA(  Alength+2 );
  vec< Float > AdjustedB(  Blength+2 );

  int first_a = max( 0, beg_a - bandwidth ) ;
  int first_b = max( 0, beg_b - bandwidth ) ;
  int last_a = min( (int) seq1.size() - 1, end_a + bandwidth ) ;
  last_a =  SWscore.lastA( end_b ) ;

  AdjustedA[0] = AdjustedScore( MAX_COST, MAX_COST, MAX_COST ) ;
  AdjustedA[1] = AdjustedScore( MAX_COST, Float(qual1[0]), Float(qual1[1]) ) ;
  AdjustedA[Alength] = AdjustedScore( Float(qual1[Alength-2]),
       Float(qual1[Alength-1]), MAX_COST ) ;
  AdjustedA[Alength+1] = AdjustedScore( MAX_COST, MAX_COST, MAX_COST ) ;
  for ( int i = 1; i < Alength - 1 ; i++ )
    AdjustedA[i+1] = AdjustedScore( Float(qual1[ i-1 ]), Float(qual1[ i ]),
         Float(qual1[ i+1 ]) ) ;

  AdjustedB[0] = AdjustedScore( MAX_COST, MAX_COST, MAX_COST ) ;
  AdjustedB[1] = AdjustedScore( MAX_COST, Float(qual2[ 0 ]), Float(qual2[ 1 ]) ) ;
  AdjustedB[Blength] = AdjustedScore( Float(qual2[ Blength-2 ]),
       Float(qual2[ Blength-1 ]), MAX_COST ) ;
  AdjustedB[Blength+1] = AdjustedScore( MAX_COST, MAX_COST, MAX_COST ) ;
  for ( int i = 1; i < Blength - 1 ; i++ )
    AdjustedB[i+1] = AdjustedScore( Float(qual2[ i-1 ]), Float(qual2[ i ]),
         Float(qual2[ i+1 ]) ) ;

  //
  //  This is the SmithWaterman code itself.
  //

  for ( int j = SWscore.firstB( first_a-1 ) -1 ; j <= SWscore.lastB( first_a-1 )+1; j++ ) {
    SWscore( first_a - 1, j ) = 0 ;
  }

  Float MUTcost  =-13 ;
  int number_of_spots = 0;
  for ( int i = first_a; i <= last_a ; i++ ) {
    if ( SWscore.firstB(i) == first_b )
      SWscore( i, SWscore.firstB( i ) - 1 ) = 0 ;
    else
      SWscore( i, SWscore.firstB( i ) - 1 ) = MAX_COST ;
    for ( int j = SWscore.firstB( i ); j <= SWscore.lastB( i ); j++ ) {
      number_of_spots ++ ;
      Float DELcost = SWscore( i-1, j ) + InDelCost( AdjustedA[ 1+ i ] , AdjustedA[ 1+ i-1 ] , AdjustedB[ 1+ j ] ) ;
      Float INScost = SWscore( i, j-1 ) + InDelCost( AdjustedB[ 1+ j ] , AdjustedB[ 1+ j-1 ] , AdjustedA[ 1+ i ] ) ;
      if ( j == 0 ) DELcost = MAX_COST ;  // ludge
      if ( i == 0 ) INScost = MAX_COST ; // Kludge
      if ( seq1[ i ] == seq2[ j ] ) {
	MUTcost = SWscore( i-1, j-1 ) + MatchValue;

	Assert( INScost > Float(-8000) ) ; // Kludge - must eliminate after debugging
	Assert( DELcost > Float(-8000) ) ;
	Assert( MUTcost > Float(-8000) ) ;
	//	Assert( MUTcost <=  min( DELcost, INScost ) ) ; // if true we needn't compute DELcost, INscost on this the common path.
	SWptr( i , j )  = Float(Match);
      }
      else {
	MUTcost = SWscore( i-1, j-1 ) + MutateCost( AdjustedA[ 1+ i ], AdjustedB[ 1+ j ] );
	SWptr( i , j )  = Float(Mutation);
      }
      if ( MUTcost <= min( DELcost, INScost ) ) {
	SWscore( i , j ) = MUTcost ;
      } else {
	if ( DELcost < INScost ) {
	  SWscore( i , j ) = DELcost ;
	  SWptr( i , j )  = Float(Deletion);
	} else {
	  SWscore( i , j ) = INScost ;
	  SWptr( i , j )  = Float(Insertion);
	}
      }
    }
    SWscore( i, SWscore.lastB( i )+1 ) = MAX_COST ;
  }

  //  SWscore and SWptr print out nicely when they involve no more than say 30-50 bases
  //  cout << std::endl << std::endl << " SWscore at the end " << std::endl << std::endl << SWscore << std::endl << std::endl ;
  //  cout << std::endl << std::endl << " SWptr at the end " << std::endl << std::endl << SWptr << std::endl << std::endl ;


  //
  //  Compute Read_Template
  //
  //    Find the lowest Swith Waterman score that includes all of either end_a or end_b
  //    i.e. SWptr( * , end_b ) and SWptr( end_a, * )
  //
  Float BestScore = MAX_COST ;
  int BestI = -1 ;
  int BestJ = -1 ;
  for ( int i = SWscore.firstA( end_b ); i <= SWscore.lastA( end_b ); i++ ) {
    if ( SWscore( i, end_b ) < BestScore ) {
      BestScore = SWscore( i, end_b ) ;
      BestI = i ;
      BestJ = end_b ;
    }
  }
  for ( int j = SWscore.firstB( end_a ); j <= SWscore.lastB( end_a ); j++ ) {
    if ( SWscore( end_a, j ) < BestScore ) {
      BestScore = SWscore( end_a, j ) ;
      BestI = end_a ;
      BestJ = j ;
    }
  }

  Assert( BestI <= Alength - 1 ) ;
  Assert( BestJ <= Blength - 1 ) ;


  //
  //    Compute Read_Template ( i.e. RTlength, RTquals, RTerrors ) ;
  //    and/or new_alignment
  //
  int i =  0 ;
  int j =  0 ;
  int RTlength = 0 ;
  qualvector RTquals( 2 * ( length + bandwidth ) ) ;  // conservatively large
  qualvector RTerrors( 2 * ( length + bandwidth ) ) ;  // conservatively large

  avector<int> gaps( 255 )  ;
  avector<int> lengths( 255 )  ;

  int currentgap = 0 ;
  bool gapover = false ;
  int numgaps = 0 ;
  int currentlength = 0 ;

  while ( ( BestI - i ) >= 0 && ( BestJ - j ) >= 0 && numgaps < 254 ) {
    switch( (int) SWptr( BestI - i , BestJ - j ) ) {
    case( (int) Match ) :
      //cout << "MATCH\n" ;
      RTquals[ RTlength ] = qual2[ BestJ - j ] ;
      j = j + 1;
      i = i + 1;
      RTerrors[ RTlength ] = BASE_CORRECT ;
      gapover = true;
      currentlength = currentlength + 1 ;
      break;
    case( (int) Mutation ) :
      //cout << "Mutation\n" ;
      RTquals[ RTlength ] = qual2[ BestJ - j ] ;
      j = j + 1;
      i = i + 1;
      RTerrors[ RTlength ] = BASE_MUT ;
      gapover = true;
      currentlength = currentlength + 1 ;
      break;
    case( (int) Deletion ) :
      //cout << "Deletion\n" ;
      RTquals[ RTlength ] = IRRELEVANT_QUAL ;
      i = i + 1;
      RTerrors[ RTlength ] = BASE_DEL ;
      //      if ( currentlength > 0 ) {  // This is how it was
      if ( currentlength > 0 || currentgap > 0 ) {
	gaps( numgaps ) = currentgap ;
	lengths( numgaps ) = currentlength ;
	currentgap = 0 ;
	currentlength = 0 ;
	gapover = false ;
        numgaps = numgaps + 1 ;
      }
      currentgap = currentgap - 1;
      break;
    case( (int) Insertion ) :
      //cout << "Insertion\n" ;
      RTquals[ RTlength ] = qual2[ BestJ - j ] ;
      j = j + 1;
      RTerrors[ RTlength ] = BASE_INS ;
      if ( currentlength > 0 || currentgap < 0 ) {
	gaps( numgaps ) = currentgap ;
	lengths( numgaps ) = currentlength ;
	currentgap = 0 ;
	currentlength = 0 ;
	gapover = false ;
        numgaps = numgaps + 1 ;
      }
      Assert( currentgap >= 0 ) ;
      currentgap = currentgap + 1;
      break;
    default:
      Assert( false ) ;
      break;
    }
    RTlength = RTlength + 1 ;
  }
  qualvector RTqualsTrue( RTlength ) ;
  vec<unsigned char> RTerrorsTrue( RTlength ) ;
  for ( int ii = 0 ; ii < RTlength; ii ++ ) {
    RTqualsTrue[ RTlength - ii - 1 ] = RTquals[ii] ;
    RTerrorsTrue[ RTlength - ii - 1 ] = RTerrors[ii] ;
  }
  if ( create_template ) {
    read_template RTemplate( RTlength, RTqualsTrue, RTerrorsTrue ) ;
    read_temp = RTemplate ;
  }
  // Assert( RTerrors[ RTlength - 1 ] == BASE_CORRECT ) ;
  beg_a = BestI - i + 1 ;
  end_a = BestI ;
  beg_b = BestJ - j + 1;  // Kludge  - isn't always correct
  end_b = BestJ ;
  Assert( beg_a >= 0 ) ;
  Assert( end_a >= 0 ) ;

  //
  //  Compute Template
  //
  gaps( numgaps ) = currentgap ;
  lengths( numgaps ) = currentlength ;
  numgaps = numgaps + 1 ;

  avector<int> newgaps( numgaps )  ;
  avector<int> newlengths( numgaps )  ;

  for ( int i = 1; i < numgaps ; i++ ) {
    newgaps( i ) = gaps( numgaps - i ) ;
    newlengths( i ) = lengths( numgaps - i - 1 ) ;
  }
  newgaps(0 ) = 0 ;
  newlengths(0 ) = lengths(numgaps-1 ) ;

  //  I don't understand why occasionally it makes more sense to have
  //  a gap at the end instead of a mutation.  Some artifact of the scoring system.
  //  In any case, it happens and the following four lines deals with it.
  if ( gaps(0) != 0 ) {
    //  Assert( false ) ;  // Does this happen on my home machine?  No.  But it does on molybdenum.  Ugh!
    //  Note:  Alignment_plus appears not to like newlengths( numgaps ) = 0
    newlengths.resize( numgaps + 1 );
    newgaps.resize( numgaps + 1 );

    newlengths( numgaps ) = 0 ;
    newgaps( numgaps ) = gaps(0) ;
    numgaps = numgaps + 1;
  }

  bool ret = false;
  if ( create_alignment ) {
    align align_ab( beg_a, beg_b, newgaps, newlengths ) ;

    //
    //  Score the two alignments
    //
    Float NewScore = ScoreAlignment( align_ab, seq1, qual1, seq2, qual2 ) ;
    Float OldScore = ScoreAlignment( new_alignment, seq1, qual1, seq2, qual2 ) ;


    if ( ( OldScore - NewScore) > (OldScore + NewScore) / Float(10.0) ) {
    }

    if ( NewScore < OldScore && numgaps < 254 ) {
      new_alignment = align_ab ;
      ret = true ;
    }
  }

  return ret ;
}



bool SmithWaterman( int bandwidth,
		      basevector &seq1,
	       basevector &seq2,
	       qualvector &qual1,
	       qualvector &qual2,
	       int &beg_a,
	       int &end_a,
	       int &beg_b,
	       int &end_b  ) {
  read_template read_temp ;
  align new_alignment ;
  return InternalSmithWaterman( bandwidth, seq1, seq2, qual1, qual2, beg_a, end_a, beg_b, end_b,
				false, read_temp, false, new_alignment ) ;
}

bool SmithWaterman( const int bandwidth,
		      const basevector &seq1,
		      const basevector &seq2,
		      const qualvector &qual1,
		      const qualvector &qual2,
		      int &beg_a,
		      int &end_a,
		      int &beg_b,
		      int &end_b,
		      read_template &read_temp ) {
  align new_alignment ;
  return InternalSmithWaterman( bandwidth, seq1, seq2, qual1, qual2, beg_a, end_a, beg_b, end_b,
				true, read_temp, false, new_alignment ) ;
}

bool SmithWaterman( const int bandwidth,
		      const basevector &seq1,
		      const basevector &seq2,
		      const qualvector &qual1,
		      const qualvector &qual2,
		      int &beg_a,
		      int &end_a,
		      int &beg_b,
		      int &end_b,
		      read_template &read_temp,
		      align &new_alignment ) {
  return InternalSmithWaterman( bandwidth, seq1, seq2, qual1, qual2, beg_a, end_a, beg_b, end_b,
				true, read_temp, true, new_alignment ) ;
}

const int MAXINSERTIONS = 3 ;
bool SmithWaterman( const int bandwidth,
		      const basevector &seq1,
		      const basevector &seq2,
		      const qualvector &qual1,
		      const qualvector &qual2,
		      align &align ) {
  int beg_a = align.pos1() ;
  int beg_b = align.pos2() ;
  int end_a = min( align.Pos1() + MAXINSERTIONS, (int) seq1.size() - 1) ;
  int end_b = min( align.Pos2() + MAXINSERTIONS, (int) seq2.size() - 1) ;

  int internal_bandwidth = max( bandwidth, abs( ( end_a - beg_a ) - ( end_b - beg_b ) ) + 1 ) ;

  internal_bandwidth = min( internal_bandwidth, maxbandwidth ) ;


  ForceAssert( seq1.size() == qual1.size() ) ;
  ForceAssert( seq2.size() == qual2.size() ) ;
  int length = min( end_a - beg_a + 1, end_b - beg_b + 1) ;

  if ( length > 0 ) {
    read_template read_temp ;
    return InternalSmithWaterman( internal_bandwidth, seq1, seq2, qual1, qual2, beg_a, end_a, beg_b, end_b,
				  false, read_temp, true, align ) ;
  } else {
    Assert( false ) ;
    return false ;
  }
}
