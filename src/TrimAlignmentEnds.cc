// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 


#include "Alignment.h"
#include "Basevector.h"
#include "math/Functions.h"
#include "ShortVector.h"
#include "TrimAlignmentEnds.h"
#include "system/Types.h"
#include "Vec.h"

// Trim back the ends of the alignment, if there is junk there,
// according to the following scheme, using the default par1 = 3, par2 = 5, par3 = 2.
//
// Walk through the alignment, starting at an end.
// while(1)
// {   Do we see par1 correct bases in a row?
//     {    yes:
//          {    Of the next par2 bases (following the par1), 
//               are at least par3 correct?
//               {    yes:
//                    {    walk backwards over any correct base;
//                         break;    }
//                    no: trim off par1 bases;    }    }
//          no: trim off par1 bases;    }    }

// The values of trim_begin and trim_end determine which end(s) of the alignment
// are trimmed.

// If RC = True, then the alignment is between rc(the first sequence) and the
// second sequence.

void TrimAlignmentEnds( align& a, Bool RC, const basevector &read, 
     const basevector& contig, Bool trim_begin, Bool trim_end, int par1, int par2,
     int par3 )
{
  // First find the error positions, relative to the beginning of
  // the alignment.

  vec<Bool> error_positions;
  error_positions.resize( a.Pos1( ) - a.pos1( ) );
  for ( unsigned int j = 0; j < error_positions.size( ); j++ )
    error_positions[j] = False;

  int pos1x = a.pos1( ), pos2x = a.pos2( );
  const avector<int> &gaps = a.Gaps( ), &lengths = a.Lengths( );
  int nblocks = a.Nblocks( );
  int p1 = pos1x, p2 = pos2x;
  for ( int j = 0; j < nblocks; j++ )
  {   
    if ( gaps(j) > 0 )
    { 
      if ( p1 - pos1x < (int) error_positions.size( ) )
	error_positions[ p1 - pos1x ] = True;
      // if ( p1 + 1 - pos1x < (int) error_positions.size( ) )
      //      error_positions[ p1 + 1 - pos1x ] = True;
      p2 += gaps(j);
    }
    if ( gaps(j) < 0 )
    {
      for ( int k = p1; k < p1 - gaps(j); k++ )
	if ( k - pos1x < (int) error_positions.size( ) )
	  error_positions[ k - pos1x ] = True;
      p1 -= gaps(j);    }
    for ( int x = 0; x < lengths(j); x++ )
    { 
      if ( ( !RC && read[p1] != contig[p2] )
	   || ( RC && 3 - read[ read.size( ) - p1 - 1 ] != contig[p2] ) )
      { 
	if ( p1 - pos1x < (int) error_positions.size( ) )
	  error_positions[ p1 - pos1x ] = True;    
      }
      ++p1;
      ++p2;   
    }   
  }

  // Compute trim amount for beginning of alignment.

  int good_in_a_row = 0, j0;
  for ( j0 = 0; j0 < (int) error_positions.size( ); j0++ )
  { 
    if ( !error_positions[j0] ) 
      ++good_in_a_row;
    else
      good_in_a_row = 0;
    if ( good_in_a_row == par1 )
    { 
      good_in_a_row = 0;
      int errors_in_five = 0;
      for ( int k = j0 + 1; 
	    k < Min( (int) error_positions.size( ), j0 + par2 + 1 ); k++ )
	if ( error_positions[k] )
	  ++errors_in_five;
      if ( errors_in_five <= par2 - par3 )
      {   
	for ( ; j0 >= 0; j0-- )
	  if ( error_positions[j0] ) 
	    break;
	++j0;
	break;
      }
    }
  }

  // Compute trim amount for end of alignment.

  good_in_a_row = 0; 
  int j1;
  for ( j1 = (int) error_positions.size( ) - 1; j1 >= 0; j1-- )
  { 
    if ( !error_positions[j1] ) 
      ++good_in_a_row;
    else
      good_in_a_row = 0;
    if ( good_in_a_row == par1 )
    {    
      good_in_a_row = 0;
      int errors_in_five = 0;
      for ( int k = j1 - 1; k >= Max( 0, j1 - par2 ); k-- )
	if ( error_positions[k] )
	  ++errors_in_five;
      if ( errors_in_five <= par2 - par3 ) 
      {
	for ( ; j1 < (int) error_positions.size( ); j1++ )
	  if ( error_positions[j1] ) 
	    break;
	--j1;
	break;
      }
    }    
  }

  // Trim the alignment.  

  if ( a.Pos1( ) - a.pos1( ) - j0 
       - ( (int) error_positions.size( ) - j1 - 1 ) < 20 ) 
    return;

  if ( j0 > 0 && trim_begin ) 
    TrimAlignmentFront( a, j0 );
  if ( j1 < (int) error_positions.size( ) - 1 && trim_end )
  {    
    a.ReverseThis( read.size( ), contig.size( ) );
    TrimAlignmentFront( a, (int) error_positions.size( ) - j1 - 1 );
    a.ReverseThis( read.size( ), contig.size( ) );    
  } 
}
