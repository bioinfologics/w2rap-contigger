// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
//


#include "Alignment.h"
#include "math/Functions.h"
#include "pairwise_aligners/GenAlignments.h"
#include "Kclock.h"
#include "kmers/KmerRecord.h"
#include "pairwise_aligners/Mutmer.h"
#include "pairwise_aligners/MutmerGraph.h"
#include "ShortVector.h"
#include "dvString.h"
#include "system/System.h"
#include "system/Types.h"
#include "Vec.h"

// #include "gzstream/gzstream.h"

template<int I, int k, int BLOCKS_PER_NODE> void
GenAlignments(
    const vecbasevector& EE,
    mutmer_graph<I, BLOCKS_PER_NODE>& M,
    int N0,
    String& aligns_file,
    ostream* log,
    int min_mutmer,
    int max_alignments,
    makealigns_method *method,
    int offset_A,
    int offset_B )
{
  longlong N = EE.size( );
  ostream* s;
  bool is_gzip_file = false;

  if (aligns_file.Contains(".gz", -1))
    is_gzip_file = true;

  procbuf* pb;
  if ( is_gzip_file ) {
    String command = "gzip -1 > " + aligns_file; 
    pb = new procbuf( command.c_str(), ios::out );
    s = new ostream( pb );
  }
  else {
    pb = 0;
    s = new ofstream(aligns_file.c_str(), ios::out|ios::binary);
  }

  ostream& aligns_out = *s;

  vec< mutmer_read_idX<I> > mid(Max(M.Counts()));
  vec<align> aligns(max_alignments);
  vec<int> errors(max_alignments);

  int aligns_length, j;
  basevector rcrd2;
  int id2_last = -1;
  for ( int l = 0; l < N; l++ ) {
    int count = M.All(l, mid);
    sort(mid.begin( ), mid.begin( ) + count);
    for ( int i = 0; i < count; i++ ) {
      for ( j = i+1; j < count; j++ ) {
	if (mid[i].ReadIdRc( ) != mid[j].ReadIdRc()) break;
      }

      int max_mutmer_len = 0;
      vec<mutmer> mm(j-i);
      for ( int r = 0; r < j-i; r++ ) {
	int pos1, pos2, len, e;
	mid[i+r].Unpack( pos1, pos2, len, e );
	mm[r].SetFrom( pos1, pos2, len, e );
	max_mutmer_len = Max( max_mutmer_len, len );
      }

      Bool RC = mid[i].Rc( );
      int id1 = l, id2 = mid[i].ReadId( );
      Bool succeed = False;
      int cur_min_mutmer = min_mutmer;

//       PRINT2( id1, id2 );

      while( !succeed ) {
	if ( !RC ) {
	  succeed = method->MutmersToAlign( mm, k, EE[id1], EE[id2], aligns, 
	      errors, aligns_length, cur_min_mutmer, log );
	}

	else {
          if ( id2 != id2_last )
	  {    rcrd2.Setsize(EE[id2].size());
	       rcrd2.ReverseComplement(EE[id2]);
               id2_last = id2;    }
	  //succeed = MutmersToAlign( mm, k, EE[id1], EErc[id2], aligns, 
	  succeed = method->MutmersToAlign( mm, k, EE[id1], rcrd2, aligns, 
	      errors, aligns_length, cur_min_mutmer, log );
	}

	if ( cur_min_mutmer > max_mutmer_len )
	  break;
       
	cur_min_mutmer += 15;
      }

      if ( aligns_length > 0 ) {
	BinWrite( aligns_out, aligns_length );

	// Fix id1 and id2.
	int save_id1 = id1;
	if ( save_id1 < N0 )
	  save_id1 += offset_A;
	else
	  save_id1 += offset_B;

	int save_id2 = id2;
	if ( save_id2 < N0 )
	  save_id2 += offset_A;
	else
	  save_id2 += offset_B;

	// Write alignment.
	for ( int q = 0; q < aligns_length; q++ ) {
	  aligns[q].Write( aligns_out, save_id1, save_id2, RC, errors[q] );
	}
      }
      i = j - 1;
    }
  }
  
  delete s;
  delete pb;
}

#define INSTANTIATE_GENALIGNMENTS(_i,_k,_bpn) \
template void GenAlignments<_i,_k,_bpn>(				\
    const vecbasevector& EE, mutmer_graph<_i, _bpn>& M, int N0,		\
    String& aligns_file, ostream* log, int min_mutmer, int max_alignments, \
    makealigns_method *method, int offset_A, int offset_B )

#define INSTANTIATE_GENALIGNMENTS_FOR_K(_K, dummy) \
  INSTANTIATE_GENALIGNMENTS(1,_K,12); \
  INSTANTIATE_GENALIGNMENTS(2,_K,12); \
  INSTANTIATE_GENALIGNMENTS(1,_K,50); \
  INSTANTIATE_GENALIGNMENTS(2,_K,50)

FOR_ALL_K( INSTANTIATE_GENALIGNMENTS_FOR_K,unused);

