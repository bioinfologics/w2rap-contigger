///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef PRINTALIGNMENT
#define PRINTALIGNMENT

#include "Alignment.h"
#include "Basevector.h"
#include "CoreTools.h"
#include "Qualvector.h"

void PrintBlanks( std::ostream& out, int n );

template<class BASEVEC>
void PrintBases( std::ostream& out, const BASEVEC& rd, int from, int to );

template<class BASEVEC1, class BASEVEC2>
void PrintVisualAlignment( Bool abbreviate, std::ostream& out, const BASEVEC1& rd1, 
     const BASEVEC2& rd2, const align& a, 
     const qualvector& scores1 = qualvector(0), 
     const qualvector& scores2 = qualvector(0), 
     int begin = 0, Bool one_frame = false, int min_score_to_abbrev = 0,
     Bool abbeviate_poor = False, float min_fract_poor = 2.0,
     Bool abbreviate_good = False, float max_fract_good = 0.05,
     Bool print_heads_and_tails = True, const Bool CtoT_special = False,
     const int pw = 80 );

// This version deletes spurious blank lines.

template<class BASEVEC1, class BASEVEC2>
void PrintVisualAlignmentClean( Bool abbreviate, std::ostream& out, const BASEVEC1& rd1, 
     const BASEVEC2& rd2, const align& a, 
     const qualvector& scores1 = qualvector(0), 
     const qualvector& scores2 = qualvector(0), 
     int begin = 0, Bool one_frame = false, int min_score_to_abbrev = 0,
     Bool abbeviate_poor = False, float min_fract_poor = 2.0,
     Bool abbreviate_good = False, float max_fract_good = 0.05,
     Bool print_heads_and_tails = True, const Bool CtoT_special = False,
     const int pw = 80 )
{
     std::ostringstream outx;
     PrintVisualAlignment( abbreviate, outx, rd1, rd2, a, scores1, scores2,
          begin, one_frame, min_score_to_abbrev, abbeviate_poor, min_fract_poor,
          abbreviate_good, max_fract_good, print_heads_and_tails, CtoT_special, pw );
     std::istringstream in( outx.str( ) );
     String line;
     vec<String> lines;
     while(1)
     {    getline( in, line );
          if ( in.fail( ) ) break;
          lines.push_back(line);    }
     for ( int i = 0; i < lines.isize( ); i++ )
     {    out << lines[i] << "\n";
          int j;
          for ( j = i + 1; j < lines.isize( ); j++ )
               if ( WhiteSpaceFree( lines[j] ).size( ) > 0 ) break;
          for ( int k = 0; k < Min( 1, j - i - 1 ); k++ )
               out << "\n";
          i = j - 1;    }    }

template<class BASEVEC1, class BASEVEC2>
void PrintVisualAlignment( Bool rd2_is_rc, Bool abbreviate, std::ostream& out, 
     const BASEVEC1& rd1, BASEVEC2 rd2, const align& a, 
     const qualvector& scores1 = qualvector(0), qualvector scores2 = qualvector(0),
     int begin = 0, Bool one_frame = false, int min_score_to_abbrev = 0,
     Bool abbeviate_poor = False, float min_fract_poor = 2.0,
     Bool abbreviate_good = False, float max_fract_good = 0.05, const int pw = 80 );

template<class BASEVEC1, class BASEVEC2>
inline void PrintVisualAlignment( Bool abbreviate, std::ostream& out, 
     const BASEVEC1& rd1, const BASEVEC2& rd2, const alignment& a, 
     const qualvector& scores1 = qualvector(0), 
     const qualvector& scores2 = qualvector(0),
     int begin = 0, Bool one_frame = false, int min_score_to_abbrev = 0,
     Bool abbreviate_poor = False, float min_fract_poor = 2.0, const int pw = 80 )
{    PrintVisualAlignment( abbreviate, out, rd1, rd2, align(packalign(a)),
          scores1, scores2, begin, one_frame, min_score_to_abbrev,
          abbreviate_poor, min_fract_poor, pw );    }

template<class BASEVEC1, class BASEVEC2>
inline void PrintVisualAlignment( Bool rd2_is_rc, Bool abbreviate, std::ostream& out, 
     const BASEVEC1& rd1, BASEVEC2 rd2, const alignment& a, 
     const qualvector& scores1 = qualvector(0),
     qualvector scores2 = qualvector(0),
     int begin = 0, Bool one_frame = false, int min_score_to_abbrev = 0,
     Bool abbreviate_poor = False, float min_fract_poor = 2.0, const int pw = 80 )
{    PrintVisualAlignment( rd2_is_rc, abbreviate, out, rd1, rd2, 
          align(packalign(a)), scores1, scores2, begin, one_frame,
          min_score_to_abbrev, abbreviate_poor, min_fract_poor, pw );    }

void PrintAlignment( std::ostream& out, const basevector& rd1, 
     const basevector& rd2, const alignment& a );
void PrintAlignment( Bool rd2_is_rc, std::ostream& out, const basevector& rd1, 
     basevector rd2, const alignment& a );
void PrintReadWithScores( basevector& B, qualvector& Q, std::ostream& out,
     const int pw = 80 );
void PrintErrorsInAlignment(std::ostream& out, const basevector& rd1, 
     const basevector& rd2, const alignment& a, const qualvector& scores1,
     const qualvector& scores2 );

#endif
