// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology

#include "Alignment.h"
#include "math/Arith.h"
#include "Basevector.h"
#include "Bitvector.h"
#include "CoreTools.h"
#include "math/Functions.h"
#include "IndexedAlignmentPlusVector.h"
#include "NQS.h"
#include "PackAlign.h"
#include "Qualvector.h"
#include "Rmr.h"

int NQS( const basevector& rd1, const basevector& rd2, const qualvector& q1, 
     const qualvector& q2, int p1, int p2, int NQS_floor1, int NQS_floor2,
     int NQS_extra )
{    if ( q1[p1] < NQS_floor1 ) return NQS_na;
     if ( q2[p2] < NQS_floor1 ) return NQS_na;
     for ( int y = 1; y <= NQS_radius; y++ )
     {    if ( q1[p1-y] < NQS_floor2 || q1[p1+y] < NQS_floor2 
               || q2[p2-y] < NQS_floor2 || q2[p2+y] < NQS_floor2 )
          {    return NQS_na;    }    }
     int ndiff = 0;
     for ( int y = -NQS_radius; y <= NQS_radius; y++ )
          if ( y != 0 && rd1[p1+y] != rd2[p2+y] ) ++ndiff;
     if ( ndiff > NQS_extra ) return NQS_na;
     if ( rd1[p1] == rd2[p2] ) return NQS_same;
     return NQS_diff;    }

int NQS0( const basevector& rd1, const basevector& rd2, int p1, int p2, 
     int NQS_extra )
{    int ndiff = 0;
     for ( int y = -NQS_radius; y <= NQS_radius; y++ )
          if ( y != 0 && rd1[p1+y] != rd2[p2+y] ) ++ndiff;
     if ( ndiff > NQS_extra ) return NQS_na;
     if ( rd1[p1] == rd2[p2] ) return NQS_same;
     return NQS_diff;    }

int NQS_vs_perf( const basevector& rd1, const basevector& rd2, const qualvector& q1, 
     int p1, int p2, int NQS_floor1, int NQS_floor2, int NQS_extra )
{    if ( q1[p1] < NQS_floor1 ) return NQS_na;
     for ( int y = 1; y <= NQS_radius; y++ )
          if ( q1[p1-y] < NQS_floor2 || q1[p1+y] < NQS_floor2 ) return NQS_na;
     int ndiff = 0;
     for ( int y = -NQS_radius; y <= NQS_radius; y++ )
          if ( y != 0 && rd1[p1+y] != rd2[p2+y] ) ++ndiff;
     if ( ndiff > NQS_extra ) return NQS_na;
     if ( rd1[p1] == rd2[p2] ) return NQS_same;
     return NQS_diff;    }

void CountNQS( 
     /* inputs: */     const alignment& a, 
                       const basevector& rd1, const basevector& rd2, 
                       const qualvector& q1, const qualvector& q2, 
     /* heuristics: */ int NQS_floor1, int NQS_floor2, int NQS_extra,
     /* outputs: */    int& look, int& see )
{
     look = 0, see = 0;
     static avector<int> gaps, lengths;
     int pos1, pos2, errors;
     a.Unpack( pos1, pos2, errors, gaps, lengths );
     unsigned int j, p1 = pos1, p2 = pos2;
     for ( j = 0; j < gaps.length; j++ )
     {    if ( gaps(j) > 0 ) p2 += gaps(j);
          if ( gaps(j) < 0 ) p1 -= gaps(j);
          for ( int x = 0; x < lengths(j); x++ )
          {    if ( x >= NQS_radius && x < lengths(j) - NQS_radius )
               {    int n = NQS( rd1, rd2, q1, q2, p1, p2, NQS_floor1, NQS_floor2, 
                         NQS_extra );
                    if ( n == NQS_diff ) ++see;
                    if ( n == NQS_same || n == NQS_diff ) ++look;    }
               ++p1;
               ++p2;    }    }
     return;    }

void CountNQS(
     /* inputs: */         const align& a,
                           const basevector& rd1, const basevector& rd2,
                           const qualvector& q1, const qualvector& q2,
     /* heuristics: */     int NQS_floor1, int NQS_floor2, int NQS_extra,
     /* outputs: */        int& look, int& see,
     /* optional input: */ const bitvector& restricted_to1,
                           const bitvector& restricted_to2 )
{
     look = 0, see = 0;
     int p1 = a.pos1( ), p2 = a.pos2( );
     for ( int j = 0; j < a.Nblocks( ); j++ )
     {    if ( a.Gaps(j) > 0 ) p2 += a.Gaps(j);
          if ( a.Gaps(j) < 0 ) p1 -= a.Gaps(j);
          for ( int x = 0; x < a.Lengths(j); x++ )
          {    if ( restricted_to1.capacity( )
                    && ( !restricted_to1[p1] || !restricted_to2[p2] ) )
               {    ++p1; ++p2;
                    continue;    }
               if ( x >= NQS_radius && x < a.Lengths(j) - NQS_radius )
               {    int n = NQS( rd1, rd2, q1, q2, p1, p2, NQS_floor1, NQS_floor2,
                         NQS_extra );
                    if ( n == NQS_diff ) ++see;
                    if ( n == NQS_same || n == NQS_diff ) ++look;    }
               ++p1; ++p2;    }    }
     return;    }

void CountNQS( 
     /* inputs: */         const align& a, 
                           const basevector& rd1, const basevector& rd2, 
     /* heuristics: */     int NQS_extra,
     /* outputs: */        int& look, int& see,
     /* optional input: */ const bitvector& restricted_to1,
                           const bitvector& restricted_to2 )
{
     look = 0, see = 0;
     int p1 = a.pos1( ), p2 = a.pos2( );
     for ( int j = 0; j < a.Nblocks( ); j++ )
     {    if ( a.Gaps(j) > 0 ) p2 += a.Gaps(j);
          if ( a.Gaps(j) < 0 ) p1 -= a.Gaps(j);
          for ( int x = 0; x < a.Lengths(j); x++ )
          {    if ( restricted_to1.capacity( )
                    && ( !restricted_to1[p1] || !restricted_to2[p2] ) )
               {    ++p1; ++p2;
                    continue;    }
               if ( x >= NQS_radius && x < a.Lengths(j) - NQS_radius )
               {    int n = NQS0( rd1, rd2, p1, p2, NQS_extra );
                    if ( n == NQS_diff ) ++see;
                    if ( n == NQS_same || n == NQS_diff ) ++look;    }
               ++p1; ++p2;    }    }
     return;    }

// ConfirmBases: given an alignment between b1 and b2, mark bases on b1 as
// confirmed if they meet the following criteria:
// (a) length(b2) >= 400;
// (b) the alignment has score at most 100;
// (c) the base is centered in an 11-base gap-free block of perfect agreement in the
//     alignment, of quality >= 25 at all the bases, and quality >= 30 at the center 
//     base, with the entire block situated at least 20 bases from the end of at 
//     least one of the reads;
// (d) rmr(alignment) <= 1.0%.
// The first and last 10 bases are never confirmed.

void ConfirmBases( bitvector& confirmed, const alignment_plus& ap, 
     const basevector& b1, const basevector& b2, const qualvector& q1,
     const qualvector& q2 )
{    int n1 = b1.size( ), n2 = b2.size( );
     if ( Float(ap.score) > Float(100) || n2 < 400 ) return;
     if ( ReciprocalMatchRate( ap, b1, b2 ) > Float(0.01) ) return;
     static vec<Bool> match, match_plus;
     match.resize_and_set( n1, False );
     match_plus.resize_and_set( n1, False );
     static basevector rd2rc;
     if ( ap.Rc2( ) ) rd2rc.ReverseComplement(b2);
     const basevector &rd1 = b1, &rd2 = ( !ap.Rc2( ) ? b2 : rd2rc );
     static qualvector q2rc;
     if ( ap.Rc2( ) ) q2rc.SetToReverseOf(q2);    
     const qualvector &Q1 = q1, &Q2 = ( !ap.Rc2( ) ? q2 : q2rc );
     const align& a = ap.a;
     int j, p1 = a.pos1( ), p2 = a.pos2( );
     for ( j = 0; j < a.Nblocks( ); j++ )
     {    if ( a.Gaps(j) > 0 ) p2 += a.Gaps(j);
          if ( a.Gaps(j) < 0 ) p1 -= a.Gaps(j);
          for ( int x = 0; x < a.Lengths(j); x++ )
          {    if ( Q1[p1] >= 25 && Q2[p2] >= 25 )
               {    int dist1 = Min( p1, n1 - p1 ), dist2 = Min( p2, n2 - p2 );
                    if ( ( dist1 >= 20 || dist2 >= 20 ) && rd1[p1] == rd2[p2] ) 
                    {    match[p1] = True;    
                         if ( Q1[p1] >= 30 && Q2[p2] >= 30 )
                              match_plus[p1] = True;    }    }
               ++p1; ++p2;    }    }
     static vec<int> dist_to_diff;
     dist_to_diff.resize(n1);
     int dc = 0;
     for ( int z = 0; z < n1; z++ )
     {    if ( !match[z] ) dc = 0;
          dist_to_diff[z] = dc;
          if ( match[z] ) ++dc;    }
     dc = 0;
     for ( int z = n1 - 1; z >= 0; z-- )
     {    if ( !match[z] ) dc = 0;
          dist_to_diff[z] = Min( dc, dist_to_diff[z] );
          if ( match[z] ) ++dc;    }
     for ( int z = 10; z < n1 - 10; z++ )
          if ( dist_to_diff[z] >= 5 && match_plus[z] ) confirmed.Set( z, True );    }

// OverlapTransitivity: given three reads and pairwise alignments between them, test
// for failure of transitivity, as described below.
//
// We look for failure of transitivity of SNPS, as follows.  Suppose reads a and b
// have a SNP, at base x on a and at base y on b.  Suppose that there is a third
// read c, and we have alignments (a,c) and (b,c) both with score <= 100, and that
// the three alignments have consistent orientation.  Further suppose that the 
// first alignment maps x to a confirmed base on c, and that it maps the 11-base 
// neighborhood of x perfectly to c.  Suppose that the analogous statement holds 
// for the second alignment.  Then the SNP is declared invalid.
//
// In a second transitivity test, suppose reads a, b, and c have pairwise alignments
// between each other, and that the orientations are consistent.
//
//                    --------------------- a
//                                -------------------------------b
//                         -------------------------- c
//
// We assume that a and b do not both extend off an end of c.  Then we compute the
// predicted overlap between a and b (via c), and compare this to the direct overlap
// between a and b.  If the difference exceeds 50, and 
// rmr(A,B) >= 2 Max{ rmr(A,C), rmr(B,C) }, then we discard all SNPs for the
// alignment of A and B.  Note that this is potentially misleading, because what we
// really mean is that we think the alignment is wrong: probably a and b are from
// nearby parts of the genome, but do not overlap as given.
//
// Upon entry, ap23 should be reversed if ap12.Rc2( ).
//
// Upon entry reads may be in memory or not.  If idmap[x] >= 0, then
// bases[ idmap[x] ] is read x.  Otherwise, we load the read from 
// source_dir/reads.fastb.
//
// Action: if untrusted SNPs are found, their positions on id1 are appended to
// "invalids".  However, if the entire alignment is untrusted (as determined by the
// first transitivity test), we instead report via all_untrusted.

void OverlapTransitivity( int id1, int id2, int id3, const alignment_plus& ap12, 
     const alignment_plus& ap13, const alignment_plus& ap23, 
     const String& source_dir, const vec<int>& idmap, const vecbasevector& bases,
     const vecbitvector& confirmed, const vec<int>& lengths, vec<int>& invalids, 
     Bool& all_untrusted, Bool verbose )
{    
     // Set up.

     all_untrusted = False;
     int len1 = lengths[id1], len2 = lengths[id2];
     static basevector rd1, rd2, rd3;
     static vecbasevector rd1v, rd2v, rd3v;
     rd1v.clear( ), rd2v.clear( ), rd3v.clear( );
     Bool reads_loaded = False;
     String fastb = source_dir + "/reads.fastb";

     // Test for inconsistent orientations.

     if ( ap12.Rc2( ) ^ ap13.Rc2( ) ^ ap23.Rc2( ) ) return;

     // Do second transitivity test.

     if ( ap13.a.pos1( ) == 0 || ap23.a.pos1( ) == 0 )
     {    if ( ap13.a.Pos1( ) == len1 || ap23.a.Pos1( ) == len2 )
          {    int indirect_overlap 
                    = IntervalOverlap( ap13.a.pos2( ), ap13.a.Pos2( ), 
                         ap23.a.pos2( ), ap23.a.Pos2( ) );
               int direct_overlap = ap12.a.Pos1( ) - ap12.a.pos1( );
               int delta = Abs( indirect_overlap - direct_overlap );
               if ( delta >= 50 )
               {    if ( idmap[id1] < 0 ) rd1v.ReadOne( fastb, id1 );
                    if ( idmap[id2] < 0 ) rd2v.ReadOne( fastb, id2 );
                    if ( idmap[id3] < 0 ) rd3v.ReadOne( fastb, id3 );
                    if ( idmap[id1] < 0 ) rd1 = rd1v[0];
                    else rd1 = bases[ idmap[id1] ];
                    if ( idmap[id2] < 0 ) rd2 = rd2v[0];
                    else rd2 = bases[ idmap[id2] ];
                    if ( idmap[id3] < 0 ) rd3 = rd3v[0];
                    else rd3 = bases[ idmap[id3] ];
                    if ( ap12.Rc2( ) ) rd2.ReverseComplement( );
                    if ( ap13.Rc2( ) ) rd3.ReverseComplement( );
                    reads_loaded = True;
                    static align a, a1, a2;
                    a.UnpackFrom(ap12.a); 
                    a1.UnpackFrom(ap13.a); 
                    a2.UnpackFrom(ap23.a);
                    Float rmrp12 = Float(100) * ReciprocalMatchRate( a, rd1, rd2 );
                    Float rmrp13 = Float(100) * ReciprocalMatchRate( a1, rd1, rd3 );
                    Float rmrp23 = Float(100) * ReciprocalMatchRate( a2, rd2, rd3 );
                    if ( rmrp12 >= Float(2) * rmrp13 && rmrp12 >= Float(2) * rmrp23 )
                    {    if (verbose)
                         {    cout << "\noverlap discrepancy:\n";
                              PRINT4( delta, id1, id2, id3 );    
                              PRINT3( rmrp12, rmrp13, rmrp23 );    }
                         all_untrusted = True;
                         return;    }    }    }    }

     // Test scores for first transitivity test.

     if ( ap13.score > 100 || ap23.score > 100 ) return;

     // Define c1, c2, c3, rd1, rd2, rd3.

     static bitvector c2rc, c3rc;
     if ( ap12.Rc2( ) )
     {    c2rc = confirmed[id2];
          c2rc.ReverseMe( );    }
     if ( ap13.Rc2( ) )
     {    c3rc = confirmed[id3];
          c3rc.ReverseMe( );    }
     const bitvector& c1 = confirmed[id1];
     const bitvector& c2 = ( !ap12.Rc2( ) ? confirmed[id2] : c2rc );
     const bitvector& c3 = ( !ap13.Rc2( ) ? confirmed[id3] : c3rc );
     if ( !reads_loaded )
     {    if ( idmap[id1] < 0 ) rd1v.ReadOne( fastb, id1 );
          if ( idmap[id2] < 0 ) rd2v.ReadOne( fastb, id2 );
          if ( idmap[id3] < 0 ) rd3v.ReadOne( fastb, id3 );
          if ( idmap[id1] < 0 ) rd1 = rd1v[0];
          else rd1 = bases[ idmap[id1] ];
          if ( idmap[id2] < 0 ) rd2 = rd2v[0];
          else rd2 = bases[ idmap[id2] ];
          if ( idmap[id3] < 0 ) rd3 = rd3v[0];
          else rd3 = bases[ idmap[id3] ];
          if ( ap12.Rc2( ) ) rd2.ReverseComplement( );
          if ( ap13.Rc2( ) ) rd3.ReverseComplement( );    }

     // Perform first transitivity test.
     
     static align a;
     a.UnpackFrom(ap12.a);
     int p1 = a.pos1( ), p2 = a.pos2( );
     for ( int j = 0; j < a.Nblocks( ); j++ )
     {    if ( a.Gaps(j) > 0 ) p2 += a.Gaps(j);
          if ( a.Gaps(j) < 0 ) p1 -= a.Gaps(j);
          for ( int x = 0; x < a.Lengths(j); x++ )
          {    if ( c1[p1] && c2[p2] 
                    && x >= NQS_radius && x < a.Lengths(j) - NQS_radius )
               {    if ( NQS0( rd1, rd2, p1, p2, 1 ) == NQS_diff )
                    {    static align a1, a2;
                         a1.UnpackFrom(ap23.a);
                         Bool invalid1 = False, invalid2 = False;
                         int p1_1 = a1.pos1( ), p2_1 = a1.pos2( );
                         for ( int j = 0; j < a1.Nblocks( ); j++ )
                         {    if ( a1.Gaps(j) > 0 ) p2_1 += a1.Gaps(j);
                              if ( a1.Gaps(j) < 0 ) p1_1 -= a1.Gaps(j);
                              if ( p1_1 + a1.Lengths(j) <= p1 )
                              {    p1_1 += a1.Lengths(j);
                                   p2_1 += a1.Lengths(j);
                                   continue;    }
                              if ( p1_1 > p1 ) break;
                              int x = p1 - p1_1;
                              p1_1 += x;
                              p2_1 += x;
                              if ( x >= NQS_radius && x < a1.Lengths(j) - NQS_radius
                                   && c3[p2_1] )
                              {    Bool mismatch = False;
                                   for ( int k = -NQS_radius; k <= NQS_radius; k++ )
                                   {    if ( rd1[p1_1+k] != rd3[p2_1+k] )
                                        {    mismatch = True;
                                             break;    }    }
                                   if ( !mismatch ) invalid1 = True;    }
                              break;    }
                         if ( !invalid1 ) 
                         {    ++p1; ++p2;
                              continue;    }
                         a2.UnpackFrom(ap23.a);
                         int p1_2 = a2.pos1( ), p2_2 = a2.pos2( );
                         for ( int j = 0; j < a2.Nblocks( ); j++ )
                         {    if ( a2.Gaps(j) > 0 ) p2_2 += a2.Gaps(j);
                              if ( a2.Gaps(j) < 0 ) p1_2 -= a2.Gaps(j);
                              if ( p1_2 + a2.Lengths(j) <= p2 )
                              {    p1_2 += a2.Lengths(j);
                                   p2_2 += a2.Lengths(j);
                                   continue;    }
                              if ( p1_2 > p2 ) break;
                              int x = p2 - p1_2;
                              p1_2 += x;
                              p2_2 += x;
                              if ( x >= NQS_radius && x < a2.Lengths(j) - NQS_radius
                                   && c3[p2_2] )
                              {    Bool mismatch = False;
                                   for ( int k = -NQS_radius; k <= NQS_radius; k++ )
                                   {    if ( rd2[p1_2+k] != rd3[p2_2+k] )
                                        {    mismatch = True;
                                             break;    }    }
                                   if ( !mismatch ) invalid2 = True;    }
                              break;    }
                         if (invalid2)
                         {    invalids.push_back(p1);
                              if (verbose)
                              {    cout << "\nsee SNP at " << id1 << "." << p1 
                                        << " vs " << id2 << "." << p2 
                                        << ", invalidated by " << id3 
                                        << "\n";    }    }    }    }
               ++p1; ++p2;    }    }    }
