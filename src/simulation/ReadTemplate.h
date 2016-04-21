// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 


#include <iostream>
#include "Qualvector.h"
#include "Vec.h"


#ifndef READ_TEMPLATE_H
#define READ_TEMPLATE_H

const unsigned char IRRELEVANT_QUAL = 97 ; 

const unsigned char BASE_CORRECT    = 0;
const unsigned char BASE_MUT        = 1;
const unsigned char BASE_INS        = 2;
const unsigned char BASE_DEL        = 3;

//
// Class read_template
// This class contains all the information about a real read
//   needed in order to replicate its pattern of errors in a
//   simulated read.
//
// More specifically, the length, quality scores, and error pattern
//   (with 0 == correct, 1 == base substitution, 2 == base insertion,
//   and 3 == base deletion) are stored in the read_template.
// Moreover, the length of trimmed sequence to the left and right
//   is stored.
// 
// Later, in CreateReadsFromTemplates, an array of templates
//   constructed from Whitehead reads, is used to create random
//   simulated reads that exhibit the error patterns of real reads.
//

const int max_read_template_length  = 10000; // No template for a read longer than 10 Kb (there is no such read).

class read_template {

 public:

  //
  // Consttuctors/ Trivial destructor
  //

  read_template() :
    length_ ( 0 ) {}
  
  read_template( int length,
		 int left_trim,
		 int right_trim,
		 qualvector quals,
		 vec< unsigned char > errors ) :
    length_ ( length     ),
    l_trim_ ( left_trim  ),
    r_trim_ ( right_trim ),
    quals_  ( quals      ),
    errors_ ( errors     ) {}
  
  read_template( int length,
		 qualvector quals,
		 vec< unsigned char > errors ) :
    length_ ( length     ),
    l_trim_ ( 0          ),
    r_trim_ ( 0          ),
    quals_  ( quals      ),
    errors_ ( errors     ) {}
  
  ~read_template() {}

  //--------------------------------------------

  //
  // Const accessors of members
  //
  unsigned char operator[]( int i ) const { return errors_[ i ];                } // The ith error
  unsigned char Qual      ( int i ) const { return quals_[ i ];                 } // The ith quality
  int           UntrimmedLength()   const { return length_;                     } // Length of the untrimmed part of the read
  int           TotalLength()       const { return length_ + l_trim_ + r_trim_; } // Total read length
  int           LeftTrim()          const { return l_trim_;                     } // Length of the left trimmed part
  int           RightTrim()         const { return r_trim_;                     } // Length of the right trimmed part

  //-----------------------------------------------------------------------------

  
  //
  // Modifiers
  //

  void SetTrim( int left, int right ) { // Setting the left and right trimmed lengths (in effect,trimming the read!)
    l_trim_ = left;
    r_trim_ = right;
  }

  void SetQual( int i, int new_qual ) { quals_[ i ] = new_qual; } // Resetting the quality score at position i


  //
  // This function generates sequence for a read, given a genome, and given this template
  //   of read error pattern. It goes to the position in the 'genome' specified by 'shift'
  //   and generates sequence, putting the exact same error pattern in the generated sequence
  //   as in the template.
  //
  void GenerateSequence( char *read,         // The generated read
			 const char *genome, // The genome
			 int shift,          // Shift from the beginning of genome
			 bool rc,            // Orientation with respect to genome
			 int &mut,        // number of mutations
			 int &ins,        // number of insertions
			 int &del,        // number of deletions
			 int &consec_mut, // number of consecutive mutations
			 int &consec_ins, // number of consecutive insertions
			 int &consec_del, // number of consecutive deletions
			 int &qual40_mut, // # mutations above qual 40
			 int &qual40_ins  // # insertions above qual 40
			 );
  //-----------------------------------------------------------------


  friend istream& operator>>( istream &in,  read_template &r );
  friend ostream& operator<<( ostream &out, read_template &r );

 private:
  int length_; // read length
  int l_trim_; // amount of the original read that was trimmed from the left
  int r_trim_; // amount of the original read that was trimmed from the right

  qualvector quals_;  // vector of read Phred qualities
  vec< unsigned char > errors_; // vector of errors ( see above for meaning of 0,1,2,3)

  // IMPORTANT: When errors_[ i ] == BASE_DEL, that means that the ( i+1 )th position
  //            of this read template corresponds to a base in the source sequence that
  //            does not appear in the read. In that case, the quality score could be
  //            ignored.

};


istream& operator>>( istream &in,  read_template &r );
ostream& operator<<( ostream &out, read_template &r );

#endif
