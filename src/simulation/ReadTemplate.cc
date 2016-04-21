// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 


#include <assert.h>
#include "random/Random.h"
#include "simulation/ReadTemplate.h"
#include "charTranslations.h"

void read_template::GenerateSequence( char *r, 
				      const char *genome,
				      int shift, // leftmost nucl of read in the genome.
				      bool rc,
				      int &mut,
				      int &ins,
				      int &del,
				      int &consec_mut,
				      int &consec_ins,
				      int &consec_del,
				      int &qual40_mut,
				      int &qual40_ins ) {
  
  int curr_read_pos   = 0;
  int curr_genome_pos = shift;

  unsigned char prev = BASE_CORRECT;

  int total_length = TotalLength(); // length_ + l_trim_ + r_trim_

  int i = ( rc ? total_length - 1 : 0 );

  while ( i >= 0 && i < total_length ) {

    if ( i >= length_ + l_trim_ || // i is in the region of r_trim_
	 i <  l_trim_ ) {          // i is in the region of l_trim_

      r[ curr_read_pos++ ] = 'X';
      ++curr_genome_pos;      
    }
    else {
      int j = i - l_trim_;

      Assert( j >= 0 && j < length_ );
     
      switch ( errors_[ j ] ) {
      case BASE_CORRECT:
	r[ curr_read_pos++ ] = genome[ curr_genome_pos++ ];
	prev = BASE_CORRECT;
	break;
	
      case BASE_MUT:
	if ( quals_[ j ] >= 40 )
	  ++qual40_mut;

	r[ curr_read_pos++ ] = num2char[ (char2num[ (int)genome[ curr_genome_pos++ ] ] + 1 + randint( 3 ) ) % 4 ];  // Breaks cxx - no randint 
	if ( prev == BASE_MUT )
	  ++consec_mut;
	prev = BASE_MUT;
	++mut;
	break;
	
      case BASE_INS:      
	if ( quals_[ j ] >= 40 )
	  ++qual40_ins;
	
	r[ curr_read_pos++ ] = num2char[ randint( 4 ) ];  // Breaks cxx - no randint
	if ( prev == BASE_INS )
	  ++consec_ins;
	prev = BASE_INS;
	++ins;
	break;
	
      case BASE_DEL:
	curr_genome_pos++;
	if ( prev == BASE_DEL )
	  ++consec_del;
	prev = BASE_DEL;      
	++del;
	break;
	
      default:
	Assert( false );
      }
    }
    if ( rc ) --i;
    else      ++i;

  }
  r[ curr_read_pos ] = '\0';
  assert( curr_read_pos + del == total_length );
  if ( (int) strlen(r ) != curr_read_pos )
    cout << curr_read_pos << std::endl << total_length << std::endl
	 << "strlen is: " << strlen( r ) << std::endl << r << std::endl << *this << std::endl;
  assert( (int)strlen( r ) + del == total_length );
  
}

istream& operator>>( istream &in, read_template &r ) {
  r.quals_.clear();
  r.errors_.clear();
  
  in >> r.length_;
  Assert( r.length_ > 0 && r.length_ < max_read_template_length );

  in >> r.l_trim_;
  in >> r.r_trim_;
  // Assert( r.l_trim_ + r.r_trim_ <= r.length_ );
  
  char buf[100];
  unsigned char temp_char;
  for ( int i = 0; i< r.length_; ++i ) {
    in >> buf;
    r.quals_.push_back( atoi( buf ) );
  }
  for ( int i = 0; i< r.length_; ++i ) {
    in >> buf;
    temp_char = (unsigned char)atoi( buf );
    Assert( buf[ 1 ] == '\0' );
    Assert( temp_char <= 3 );
    r.errors_.push_back( temp_char );
  }
  return in;
}

ostream& operator<<( ostream &out, read_template &r ) {

  out << r.length_ << std::endl;
  out << r.l_trim_ << " ";
  out << r.r_trim_ << std::endl;

  for ( int i = 0; i< r.length_; ++i )
    out << (int) r.quals_[ i ] << ' ';
  out << std::endl;

  for ( int i = 0; i< r.length_; ++i )
    out << (int) r.errors_[ i ] << ' ';
  out << std::endl;

  return out;
}
