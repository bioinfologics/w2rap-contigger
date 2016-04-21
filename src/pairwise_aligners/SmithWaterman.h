// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 


#ifndef SMITHWATERMAN
#define SMITHWATERMAN

#include "Alignment.h"
#include "Basevector.h"
#include "PackAlign.h"
#include "Qualvector.h"
#include "simulation/ReadTemplate.h"

//  
// class SmithWaterman.h
//

bool SmithWaterman( const int bandwidth, 
		      const basevector &seq1, 
		      const basevector &seq2, 
		      const qualvector &qual1,
		      const qualvector &qual2,
		      int &beg_a, 
		      int &end_a, 
		      int &beg_b, 
		      int &end_b  ) ;

bool SmithWaterman( const int bandwidth,
		      const basevector &seq1, 
		      const basevector &seq2, 
		      const qualvector &qual1,
		      const qualvector &qual2,
		      int &beg_a, 
		      int &end_a, 
		      int &beg_b, 
		      int &end_b,
		      read_template &read_temp ) ;

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
		      align &al ) ;

bool SmithWaterman( const int bandwidth,
		      const basevector &seq1, 
		      const basevector &seq2, 
		      const qualvector &qual1,
		      const qualvector &qual2,
		      align &al ) ;

inline bool SmithWaterman( const int bandwidth,
		      const basevector &seq1, 
		      const basevector &seq2, 
		      const qualvector &qual1,
		      const qualvector &qual2,
		      int &beg_a, 
		      int &end_a, 
		      int &beg_b, 
		      int &end_b,
		      read_template &read_temp,
		      alignment &al )
{    packalign p = al;
     align a( p );
     bool answer = SmithWaterman( bandwidth, seq1, seq2, qual1, qual2, beg_a, end_a,
          beg_b, end_b, read_temp, a );
     al = alignment( packalign(a), 0 );
     return answer;    }

inline bool SmithWaterman( const int bandwidth,
		      const basevector &seq1, 
		      const basevector &seq2, 
		      const qualvector &qual1,
		      const qualvector &qual2,
		      alignment &al )
{    packalign p = al;
     align a( p );
     bool answer = SmithWaterman( bandwidth, seq1, seq2, qual1, qual2, a );
     al = alignment( packalign(a), 0 );
     return answer;    }

#endif
