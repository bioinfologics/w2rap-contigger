// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 


// I kludged this terribly by copying a subset of the old alignment class definition 
// here, now renamed alignment_old.  Blame me.  -Dave

// Copied stuff begins here.

#include "system/Assert.h"
#include "Basevector.h"
#include "math/Functions.h"
#include "pairwise_aligners/Mutmer.h"
#include "ShortVector.h"
#include "system/Types.h"
#include "Vec.h"

class gaplen {

     public:
  
     gaplen( ) { }

     gaplen(int gap, int length) : gap_(gap), length_(length) { }

     int Gap( ) const { return gap_; }
     int Length( ) const { return length_; }

     void Unpack(int& gap, int& length)
     {    gap = gap_;
          length = length_;    }

     private:

     int gap_, length_;

};

class alignment_old {

     public:

     shortvector<gaplen> block;

     alignment_old( ) { }

     void SetPosPosErrors( int pos1, int pos2, int errors )
     {    pos1_ = pos1;

          pos2_ = pos2;
          errors_ = errors;    }

     int pos1( ) const { return pos1_; }

     int pos2( ) const { return pos2_; }

     int Pos1( ) const
     {    shortvector<int> gaps(0), lengths(0);
          int pos1, pos2, errors;
          Unpack( pos1, pos2, errors, gaps, lengths );
          for ( int j = 0; j < gaps.length; j++ )
          {    if ( gaps(j) < 0 ) pos1 -= gaps(j);
               pos1 += lengths(j);    }
          return pos1;    }

     int Pos2( ) const
     {    shortvector<int> gaps(0), lengths(0);
          int pos1, pos2, errors;
          Unpack( pos1, pos2, errors, gaps, lengths );
          for ( int j = 0; j < gaps.length; j++ )
          {    if ( gaps(j) > 0 ) pos2 += gaps(j);
               pos2 += lengths(j);    }
          return pos2;    }

     alignment_old( int pos1, int pos2, int errors, const shortvector<int>& gaps,
          const shortvector<int>& lengths ) : block( gaps.length )
     {    SetPosPosErrors( pos1, pos2, errors );
          for ( int i = 0; i < block.length; i++ )
               block(i) = gaplen( gaps(i), lengths(i) );    }

     void Unpack( int& pos1, int& pos2, int& errors, shortvector<int>& gaps,
          shortvector<int>& lengths ) const
     {    pos1 = pos1_;
          pos2 = pos2_;
          errors = errors_;
          gaps.Setsize( block.length );
          lengths.Setsize( block.length );
          for ( int i = 0; i < block.length; i++ )
               block(i).Unpack( gaps(i), lengths(i) );    }

     int Errors( ) const { return errors_; }

     private:

     int pos1_, pos2_, errors_;

};

int Bandwidth( alignment_old& a )
{    const int add_to_bandwidth = 8; // heuristic
     int low = 0, high = 0, gap_total = 0;
     for ( unsigned int l = 0; l < a.block.length; l++ )
     {    gap_total += a.block(l).Gap( );
          low = Min( low, gap_total );
          high = Max( high, gap_total );    }
     return Max( Abs(low), Abs(high) ) + add_to_bandwidth;    }

// Copied stuff ends here.

/*
 * TODO
 *  Fix errors_ with ActualErrors(...);
 */ 


/*
 * Cambridge, March 12, 2001 (last revised: March 29, 2001.)
 *
 * PerfectAlignment.cc
 *
 * It tries to improve a given initial alignment.
 *
 * The idea is we start with an inital rough alignment, which may or may
 * not be improved by any algorithm. Applying Smith-Waterman would produce
 * the desired effect, but it would be too slow. Therefore this is what I
 * do: I look at the initial alignment, and concentrate on its "bad regions".
 * 
 * A bad region is a subsequence of the alignment with mismatches and/or
 * gaps, and it depends on the parameter h_min_matches. Each badregion is
 * characterized by a starting and ending point (on both the reads), so that:
 * there are exactly h_min_matches correct matches at both ends of the local
 * alignment; and there are no subsequences with 2*h_min_matches or more
 * consecutive matches inside the local alignment (the factor 2 is used to
 * avoid the overlapping of the badregions.)
 *
 * Each bad region is then turned into a local alignment, and eventually
 * improved with Smith-Waterman.
 */
#include "Alignment.h"
#include "Basevector.h"
#include "math/Functions.h"
#include "Overlap.h"
#include "PackAlign.h"
#include "pairwise_aligners/PerfectAlignment.h"
#include "ScoreAlignment.h"
#include "ShortVector.h"
#include "system/System.h"
#include "system/Types.h"
#include "Vec.h"
#include "pairwise_aligners/SmithWaterman.h"
#include "PrintAlignment.h"

/*
 * Defines.
 *
 *  PERF_ALIGN_FAST_PRINT:
 *   prints many messages (for debugging).
 *
 *  PERF_ALIGN_FAST_NO_MISMATCHES_IN:  (now a parameter)
 *   do not even try to change a local alignment, if it contains one ore
 *   more mismatches.
 *
 *  PERF_ALIGN_FAST_NO_MISMATCHES_OUT:  (now a parameter)
 *   disregard a local alignment change, if it is returned with one or more
 *   mismatches. This may be used in combination with  PERF_ALIGN_FAST_NO_
 *   MISMATCHES_IN, for instance: by using both, you only address 2+ gaps
 *   (no mismatches), and if you do not want to change the local alignment
 *   if a new mismatch has been brought in.
 */
// #define PERF_ALIGN_FAST_PRINT


/*
 * Heuristic.
 *
 *  h_min_matches:
 *   number of consecutive matches.
 *
 *  h_allowed_mismatches:
 *   maximum number of mismatches in a local alignment (do not try to
 *   improve alignment if there are less than or equal to h_allowed_mismatches
 *   mismatches.)
 */
// Value of h_min_matches changed from 8 to 12, 3/27/03, DJ.
const int h_min_matches = 12;
const int h_allowed_mismatches = 3;



/*
 * PerfectAlignment
 *
 * Given initial reads, quality scores, alignment, and alignment score,
 * it looks for badregions in the alignment, it builds local alignments
 * around them, and it calls Smith-Waterman on (some of) the local
 * alignments.
 *
 * Legenda:
 *  rd1: first read;
 *  q1: quality scores for rd1;
 *  rd2: second read;
 *  q2: quality scores for rd2;
 *  al: initial alignment;
 *  score: alignment score.
 *
 * Return:
 *  It ovewrites the alignment and the alignment score with the new data
 *  (if improved!) and it returns an integer which absolute value is
 *  the number of badregions which have been accessed by Smith-Waterman
 *  (a few of the badregions are obviously not improvable, and are
 *  therefore not tocuhed.) Remark: the returned value should never
 *  be negative, but if the score after the local Smith-Waterman is worst
 *  than the original score, then I do not touch alignment and alignment
 *  score, and return negative.
 */
int PerfectAlignment( const basevector& rd1, 
			    const qualvector& q1,
			    const basevector& rd2, 
			    const qualvector& q2, 
			    align& al, 
			    float& score,
                            Bool PERF_ALIGN_FAST_NO_MISMATCHES_IN,
                            Bool PERF_ALIGN_FAST_NO_MISMATCHES_OUT )
{

     // Generate an old-style alignment, al_old.

     shortvector<int> gaps_old( al.Nblocks( ) ), lengths_old( al.Nblocks( ) );
     for ( int i = 0; i < al.Nblocks( ); i++ )
     {    gaps_old(i) = al.Gaps(i);
          lengths_old(i) = al.Lengths(i);    }
     alignment_old al_old( al.pos1( ), al.pos2( ), 0, gaps_old, lengths_old );
     

#ifdef PERF_ALIGN_FAST_PRINT
  PrintVisualAlignment(True, cout, rd1, rd2, al, q1, q2);      
#endif

  shortvector<gaplen>& block = al_old.block;
  int nBlocks = block.length;
  vec<int> gap(nBlocks);
  vec<int> length(nBlocks);
  for (int ii=0; ii<nBlocks; ++ii)
    block(ii).Unpack( gap[ii], length[ii] );

  // badregion_start (_end): vector of start (respectively end) points
  //  of bad regions. Remark: each badregion starts and end with h_min_matches
  //  consecutive correct matches, and there are no other 2*h_min_matches or
  //  longer consecutive matches between start and end. The idea is that each
  //  badregion is a localized area where there are mismatches, with a buffer
  //  of correct matches at the sides.
  // n_errors: number of mismatches (gaps are not counted.)
  vec<int> badregion_start_rd1;
  vec<int> badregion_end_rd1;
  vec<int> badregion_start_rd2;
  vec<int> badregion_end_rd2;
  vec<int> n_errors;

  // counter: number of consecutive exact matches, at this point (included.) If
  //  counter = -1, then we stepped onto a mismatch or gap.
  int counter = 0;
  int loc_errors = 0;
  bool inside_badregion = false;
  int pos1 = al_old.pos1();
  int pos2 = al_old.pos2();
  for (int ii=0; ii<nBlocks; ++ii)
    {
      if ( gap[ii] > 0 )
	pos2 += gap[ii];
      else
	pos1 += -gap[ii];
      
      int initial_pos1 = pos1;
      int initial_pos2 = pos2;

      if ( ii > 0 ) 
	counter = -1;
      for (int jj=0; jj<length[ii]; ++jj)
	{
	  // Update counter.
	  if ( as_base(rd1[pos1]) == as_base(rd2[pos2]) )
	    {
	      if ( counter < 0 )
		counter = 1; 
	      else
		++counter;
	    }
	  else
	    {
	      counter = -1;
	      ++loc_errors;
	    }

	  // Four cases: inside_badregion true or false, counter<0 or >0. Do
	  //  something only in two cases.
	  int  loc_pos1, loc_pos2;
	  
	  if ( false == inside_badregion && counter < 0 )
	    {
	      // if nBlock=0, check if we are too close to the beginning
	      //  of alignment.
	      if ( 0 == ii && ( pos1 - h_min_matches < initial_pos1 ) )
		{
		  loc_pos1 = initial_pos1;
		  loc_pos2 = initial_pos2;
		}
	      else
		{
		  loc_pos1 = pos1 - h_min_matches;
		  loc_pos2 = pos2 - h_min_matches;
		}
	      
	      badregion_start_rd1.push_back( loc_pos1 );
	      badregion_start_rd2.push_back( loc_pos2 );
	      inside_badregion = true;
	    }
	  
	  if ( true == inside_badregion && counter > 0 )
	    {
	      if ( 2*h_min_matches == counter )
		{
		  badregion_end_rd1.push_back( pos1 - h_min_matches );
		  badregion_end_rd2.push_back( pos2 - h_min_matches );
		  inside_badregion = false;
		  n_errors.push_back( loc_errors );
		  loc_errors = 0;
		  counter -= h_min_matches;
		}
	    }

	  ++pos1;
	  ++pos2;
	} // next step in alignment.
      
      // At the end of a gaplength. Beware: now pos1 and pos2 are on the first step
      //  of the next gaplen!
      if ( ii < nBlocks - 1 )
	{
	  if ( false == inside_badregion )
	    {
	      int start1, start2;

	      if ( h_min_matches > length[ii] )
		{
		  start1 = pos1 - length[ii];
		  start2 = pos2 - length[ii];
		}
	      else
		{
		  start1 = pos1 - h_min_matches;
		  start2 = pos2 - h_min_matches;
		}		

	      badregion_start_rd1.push_back( start1 );
	      badregion_start_rd2.push_back( start2 );
	      inside_badregion = true;	      
	    }
	}
      else
	{
	  // Last gaplength.
	  if ( true == inside_badregion )
	    {
	      badregion_end_rd1.push_back( pos1 - 1 );
	      badregion_end_rd2.push_back( pos2 - 1 );
	      n_errors.push_back( loc_errors );
	    }
	}
    } // next gaplength.
  
  // The five vectors below must have the same size.
  int n_badregions =
    (int) ( badregion_start_rd1.size() + badregion_end_rd1.size() +
	    badregion_start_rd2.size() + badregion_end_rd1.size() +
	    n_errors.size() ) / 5; 
  
  ForceAssert ( (int) badregion_start_rd1.size() == n_badregions ||
		(int) badregion_end_rd1.size() == n_badregions ||
		(int) badregion_start_rd2.size() == n_badregions ||
		(int) badregion_end_rd2.size() == n_badregions ||
		(int) n_errors.size() == n_badregions );
  
  // Store score for original alignment,
  Float original_score = ScoreAlignment(al, rd1, q1, rd2, q2);

  // loc_*: for each badregion there are local alignment, local read1, local read2,
  //  and corresponding local quality scores. 
  vec< alignment_old > loc_al(n_badregions);
  vec< basevector > loc_rd1(n_badregions);
  vec< basevector > loc_rd2(n_badregions);
  vecqualvector loc_q1(n_badregions);
  vecqualvector loc_q2(n_badregions);

  // Must resize shortvector block for each local alignment.
  for (int kk=0; kk<n_badregions; ++kk)
    loc_al[kk].block.resize(0);

  // Build local alignments from the original alignment.
  //  kk: the local alignment;
  //  loc_pos1: position on read1 (it moves from the initial point of
  //   one gaplen to the next.)
  int kk = 0;
  int loc_pos1 = al_old.pos1();
  int loc_gap = 0;
  int loc_length = 0;
  for (int ii=0; ii<nBlocks; ++ii)
    {
      // There might be no badregions.
      if ( kk >= n_badregions )
	break;

      int tot_length_rd1;
      if ( gap[ii] < 0 )
	tot_length_rd1 = -gap[ii] + length[ii];
      else
	tot_length_rd1 = length[ii];

      if ( ii > 0 )
	{
	  loc_gap = gap[ii];
	  
	  if ( badregion_end_rd1[kk] < loc_pos1 + tot_length_rd1 )
	    {
	      if ( gap[ii] < 0 )
		loc_length = badregion_end_rd1[kk] - loc_pos1 + gap[ii] + 1;
	      else
		loc_length = badregion_end_rd1[kk] - loc_pos1 + 1;
	      loc_al[kk].block.Append(gaplen(loc_gap, loc_length));

	      if ( ++kk >= n_badregions )
		break;
	    }
	  else
	    {
	      loc_length = length[ii];
	      loc_al[kk].block.Append(gaplen(loc_gap, loc_length));
	      
	      loc_pos1 += tot_length_rd1;
	      continue;
	    }
	}

      while ( badregion_start_rd1[kk] < loc_pos1 + tot_length_rd1 )
	{
	  loc_gap = 0;
	  
	  if ( badregion_end_rd1[kk] < loc_pos1 + tot_length_rd1 )
	    {
	      loc_length = badregion_end_rd1[kk] - badregion_start_rd1[kk] + 1;
	      loc_al[kk].block.Append(gaplen(loc_gap, loc_length));

	      if ( ++kk >= n_badregions )
		break;
	    }
	  else
	    {
	      loc_length = loc_pos1 + tot_length_rd1 - badregion_start_rd1[kk];
	      loc_al[kk].block.Append(gaplen(loc_gap, loc_length));

	      loc_pos1 += tot_length_rd1;
	      break;
	    }
	}
    }

  // Check badregions and eventually apply Smith-Waterman.
  int n_skipped_badregions = 0;
  for (int kk=0; kk<n_badregions; ++kk)
    {
      // Unpack vector of gaplens.
      int n_loc_al_blocks = loc_al[kk].block.length;
      avector<int> loc_al_gap(n_loc_al_blocks);
      avector<int> loc_al_length(n_loc_al_blocks);
      for (int ii=0; ii<n_loc_al_blocks; ++ii)
	loc_al[kk].block(ii).Unpack(loc_al_gap(ii), loc_al_length(ii) );
 
      // tot_badness is the sum of mismatches and gaps.
      int tot_badness = n_errors[kk];
      for (int ii=0; ii<n_loc_al_blocks; ++ii)
	tot_badness += abs(loc_al_gap(ii));

      // local alignment initial score.
      basevector& seq1 = loc_rd1[kk];
      basevector& seq2 = loc_rd2[kk]; 
      qualvector& qual1 = loc_q1[kk];
      qualvector& qual2 = loc_q2[kk];
      alignment_old& align_old = loc_al[kk];
      
      // List of cases we want to skip.
      {	  
	// Don't do tot_badness = 1.
	if ( tot_badness < 2 )
	  {
	    ++n_skipped_badregions;
	    continue;
	  }

	// Skip three mismatches or less, with no gaps.
	if ( 1 == n_loc_al_blocks )
	  if ( n_errors[kk] < h_allowed_mismatches + 1 )
	    {
	      ++n_skipped_badregions;
	      continue;
	    }
	
	// Skip three mismatches or less, with one gaps.
	if ( 2 == n_loc_al_blocks )
	  if ( loc_al_gap(1) < 2 && n_errors[kk] < h_allowed_mismatches + 1 )
	    {
	      ++n_skipped_badregions;
	      continue;
	    }
	// Skip one or more mismatches. Remark: this is a local "overkill",
	//  since it skips all local alignments with one or more mismatches
	//  (it just leaves the 2+ gaps.)
	if ( PERF_ALIGN_FAST_NO_MISMATCHES_IN && n_errors[kk] > 0 )
	  {
	    ++n_skipped_badregions;
	    continue;
	  }

      } // End of list of cases to skip (go with Smith-Waterman.)
      
      // Resize reads and quality scores.
      int loc_size1 = badregion_end_rd1[kk] - badregion_start_rd1[kk] + 1;
      int loc_size2 = badregion_end_rd2[kk] - badregion_start_rd2[kk] + 1;	
      loc_rd1[kk].resize( loc_size1 );
      loc_rd2[kk].resize( loc_size2 );
      loc_q1[kk].resize( loc_size1 );
      loc_q2[kk].resize( loc_size2 );

      for (int ii=0; ii<loc_size1; ++ii)
 	{
	  loc_rd1[kk].Set(ii, rd1[ badregion_start_rd1[kk] + ii]);
	  loc_q1[kk][ii] = q1[badregion_start_rd1[kk] + ii];
 	}
      
      for (int ii=0; ii<loc_size2; ++ii)
 	{
	  loc_rd2[kk].Set(ii, rd2[ badregion_start_rd2[kk] + ii]);
	  loc_q2[kk][ii] = q2[badregion_start_rd2[kk] + ii];
 	}

      // Set pos1, pos2, and errors.
      loc_al[kk].SetPosPosErrors(0, 0, n_errors[kk]);

      // Set bandwidth and call Smith Waterman.
      int bandwidth = Bandwidth(align_old);

      align ali( 0, 0, loc_al_gap, loc_al_length );

#ifdef PERF_ALIGN_FAST_PRINT
      cout << "Before local S-W:";
      PrintVisualAlignment(True, cout, seq1, seq2, ali, qual1, qual2);      
#endif
      int old_Pos1 = ali.Pos1();
      int old_Pos2 = ali.Pos2();
      Float loc_old_score = ScoreAlignment(ali, seq1, qual1, seq2, qual2);
      SmithWaterman(bandwidth, seq1, seq2, qual1, qual2, ali);    
      Float loc_new_score = ScoreAlignment(ali, seq1, qual1, seq2, qual2);
#ifdef PERF_ALIGN_FAST_PRINT      
      cout << "After local S-W:";
      PrintVisualAlignment(True, cout, seq1, seq2, ali, qual1, qual2);
      cout << "old score: " << loc_old_score << "\n";
      cout << "new score: " << loc_new_score << "\n\n";
#endif      

      // Unpack the new style alignment (align).

      int pos1q = ali.pos1( ), pos2q = ali.pos2( );
      loc_al_gap = ali.Gaps( );
      loc_al_length = ali.Lengths( );
      n_loc_al_blocks = ali.Nblocks( );

      // If alignment is not improved.
      if ( loc_new_score >= loc_old_score )
	continue;

      // Check if loc_al[kk] pos' and Pos' did not change.
      if ( pos1q != 0 ||
	   pos2q != 0 ||
	   ali.Pos1() != old_Pos1 ||
	   ali.Pos2() != old_Pos2 )
	{
#ifdef PERF_ALIGN_FAST_PRINT
	  cout << std::endl
	       << "Local alignment pos and/or Pos changed, throw away change."
	       << std::endl;
	  cout << "\tpos1\tpos2\tPos1\tPos2"
	       << std::endl
	       << "Old:\t"
	       << "0" << "\t" 
	       << "0" << "\t" 
	       << old_Pos1 << "\t" 
	       << old_Pos2
	       << std::endl
	       << "New:\t" 
	       << pos1q << "\t" 
	       << pos2q << "\t" 
	       << ali.Pos1() << "\t" 
	       << ali.Pos2() 
	       << std::endl << std::endl;
#endif
	  // For now just don't do anything.
	  continue;
	}

      //  Disregard changes which bring in mismatches (we want only to
      //   move gaps, for now.)
      if (PERF_ALIGN_FAST_NO_MISMATCHES_OUT)
      {
      int t_errors = 0;
      int t_pos1 = 0;
      int t_pos2 = 0;

      for (int tt=0; tt<n_loc_al_blocks; ++tt)
	{
	  if ( loc_al_gap(tt) < 0 )
	    t_pos1 += -loc_al_gap(tt);
	  else
	    t_pos2 += loc_al_gap(tt);

	  for (int mm=0; mm<loc_al_length(tt); ++mm)
	    {
	      if ( as_base(seq1[t_pos1]) != as_base(seq2[t_pos2]) )
		++t_errors;
	      
	      ++t_pos1;
	      ++t_pos2;
	    }
	}

      loc_al[kk].SetPosPosErrors(0, 0, t_errors);      
      if ( t_errors > 0 )
	{
	  // For now just don't do anything.
	  continue;
	}
      }

      // Replace the old alignment with the new one. BEWARE: this is not
      //  efficient (we should replace all loc_al's at once, and not
      //  one by one. Anyways, there are three alignments in play:
      //   al: initial alignment;
      //   loc_al[kk]: local alignment kk;
      //   t_al: a temporary alignment (the result.)
      alignment_old t_al;
      t_al.block.resize(0);

      // t_al_pos1: can only jump from the beginning (first read) of one
      //  gaplen of t_al to the next.
      int t_al_pos1 = al_old.pos1();
      int temp_pos1 = al_old.pos1();

      // jj: gaplen jj on the alignment.
      int jj_loc_al = 0;	      
      int jj_al = 0;
      int temp_jj_al = 0;

      // Set t_al = al until badregion_start_rd1[kk]..
#ifdef PERF_ALIGN_FAST_PRINT
      cout << "Replacing alignment - Step 1 (before local alignment.)\n";
#endif
      bool do_continue = true;
      while ( do_continue && ( jj_al < nBlocks ) )
	{
	  int tot_length = length[jj_al];
	  if ( gap[jj_al] < 0 )
	    tot_length += -gap[jj_al];
	  
	  if ( t_al_pos1 + tot_length <= badregion_start_rd1[kk] )
	    {
#ifdef PERF_ALIGN_FAST_PRINT
	      cout << "Adding gaplength " 
		   << gap[jj_al] << "/" <<  length[jj_al] << "\n";
#endif
	      t_al.block.Append( gaplen(gap[jj_al], length[jj_al]) );
	      t_al_pos1 += tot_length;
	      ++jj_al;
	    }
	  else
	    {
	      do_continue = false;
	      
	      // temp_pos1 (temp_jj_al) is the pos on rd1 at the beginning of
	      //  (is the) gaplen where badregion_end_rd1[kk] is.
	      temp_jj_al = jj_al;
	      temp_pos1 = t_al_pos1;
	      bool temp_do_continue = true;
	      while ( temp_do_continue && ( temp_jj_al < nBlocks - 1 ) )
		{
		  int temp_length = length[temp_jj_al];
		  if  ( gap[temp_jj_al] < 0 )
		    temp_length += -gap[temp_jj_al];

		  if ( temp_pos1+temp_length <= badregion_end_rd1[kk] )
		    {
		      temp_pos1 += temp_length;
		      ++temp_jj_al;
		      
		    }
		  else
		    temp_do_continue = false;
		}
	    }
	}
      
      // We are now on the (changed) local alignment: plug it in t_al.
#ifdef PERF_ALIGN_FAST_PRINT
	  cout << "Replacing alignment - Step 2 (local alignment.)\n";
#endif       

      if ( gap[jj_al] < 0 )
	t_al_pos1 += -gap[jj_al];
      int initial_length = badregion_start_rd1[kk] - t_al_pos1;

      int loc_length = length[temp_jj_al];
      if ( gap[temp_jj_al] < 0 )
	loc_length += -gap[temp_jj_al];      
      int final_length = temp_pos1 + loc_length - badregion_end_rd1[kk] - 1;
      
      jj_loc_al = 0;
      while ( jj_loc_al < n_loc_al_blocks )
	{
	  int add_gap;
	  int add_length;

	  add_length = loc_al_length(jj_loc_al);
	  if ( 0 == jj_loc_al )
	    {
	      add_gap = gap[jj_al];
	      add_length += initial_length;
	    }
	  else
	    add_gap = loc_al_gap(jj_loc_al);
	  
	  if ( n_loc_al_blocks - 1 == jj_loc_al )
	    add_length += final_length;

#ifdef PERF_ALIGN_FAST_PRINT
	  cout << "Adding gaplength " 
	       << add_gap << "/" << add_length << "\n";
#endif
	  t_al.block.Append( gaplen( add_gap, add_length ) );
	  ++jj_loc_al;
	}
     
      // Final step in the buildibg of t_al,  t_al=al from here on.
#ifdef PERF_ALIGN_FAST_PRINT
	  cout << "Replacing alignment - Step 3 (after local alignment.)\n";
#endif
      for (int jj=temp_jj_al + 1; jj<nBlocks; ++jj)
	{
#ifdef PERF_ALIGN_FAST_PRINT
	  cout << "Adding gaplength " 
	       << gap[jj] << "/" << length[jj] << "\n";
#endif
	  t_al.block.Append(gaplen(gap[jj], length[jj]) );
	}
  
      // BEWARE: al_errors_ has not been updated!
      //  Replace al with t_al
      al_old.block.resize(0);
      nBlocks = t_al.block.length;
      gap.resize(nBlocks);
      length.resize(nBlocks);

      vec<int> temp_gap(nBlocks);
      vec<int> temp_length(nBlocks);
      for (int ii=0; ii<nBlocks; ++ii)
	{
	  t_al.block(ii).Unpack( temp_gap[ii], temp_length[ii] );
	  al_old.block.Append( t_al.block(ii) );
	  gap[ii] = temp_gap[ii];
	  length[ii] = temp_length[ii];
	}

      for (int ii=0; ii<nBlocks; ++ii)
	al_old.block(ii).Unpack( gap[ii], length[ii] );

      // Now replace the new style alignment al by the old style alignment al_old.

      avector<int> newgaps(gap.size( )), newlengths(length.size( ));
      for ( unsigned int i = 0; i < gap.size( ); i++ )
      {    newgaps(i) = gap[i];
           newlengths(i) = length[i];    }
      al.Set( al_old.pos1( ), al_old.pos2( ), newgaps, newlengths );
      
#ifdef PERF_ALIGN_FAST_PRINT
      cout << "Alignment replaced.";
      PrintVisualAlignment(True, cout, rd1, rd2, al, q1, q2);      
#endif
      
    } // next badregion.
  
  // Calculate new score (if need be.)
  int n_improved_badregions = n_badregions - n_skipped_badregions;
  Float new_score;
  if (  n_improved_badregions > 0 )
    new_score = ScoreAlignment(al, rd1, q1, rd2, q2);
  else
    new_score = original_score;

  // Return > 0 iff the alignment has been improved.
  int ret_value;
  if ( new_score < original_score )
    ret_value = n_improved_badregions;
  else
    ret_value = -n_improved_badregions;

  return ret_value;
}
