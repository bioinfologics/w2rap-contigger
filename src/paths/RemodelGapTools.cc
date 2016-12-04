///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "Basevector.h"
#include "ParallelVecUtilities.h"
#include "kmers/KmerRecord.h"
#include "paths/RemodelGapTools.h"

template<int K> void MakeKmerLookup0( const vecbasevector& unibases,
     vec< triple<kmer<K>,int,int> >& kmers_plus )
{    vec<int64_t> starts;
     starts.push_back(0);
     for ( size_t i = 0; i < unibases.size( ); i++ )
     {    const basevector& u = unibases[i];
          starts.push_back( starts.back( ) + Max( 0, u.isize( ) - K + 1 ) );    }
     kmers_plus.resize( starts.back( ) );
     #pragma omp parallel for
     for ( size_t i = 0; i < unibases.size( ); i++ )
     {    const basevector& u = unibases[i];
          kmer<K> x;
          for ( int j = 0; j <= u.isize( ) - K; j++ )
          {    int64_t r = starts[i] + j;
               x.SetToSubOf( u, j );
               kmers_plus[r].first = x;
               kmers_plus[r].second = i;
               kmers_plus[r].third = j;    }    }
     ParallelSort(kmers_plus);    }

template void MakeKmerLookup0( const vecbasevector& tigs,
     vec< triple< kmer<12>,int,int> >& kmers_plus );

template void MakeKmerLookup0( const vecbasevector& tigs,
     vec< triple< kmer<16>,int,int> >& kmers_plus );

template void MakeKmerLookup0( const vecbasevector& tigs,
     vec< triple< kmer<20>,int,int> >& kmers_plus );

template void MakeKmerLookup0( const vecbasevector& tigs,
     vec< triple< kmer<24>,int,int> >& kmers_plus );

template void MakeKmerLookup0( const vecbasevector& tigs,
     vec< triple< kmer<28>,int,int> >& kmers_plus );

template void MakeKmerLookup0( const vecbasevector& tigs,
     vec< triple< kmer<32>,int,int> >& kmers_plus );

template void MakeKmerLookup0( const vecbasevector& tigs,
     vec< triple< kmer<40>,int,int> >& kmers_plus );

template void MakeKmerLookup0( const vecbasevector& tigs,
     vec< triple< kmer<48>,int,int> >& kmers_plus );

template void MakeKmerLookup0( const vecbasevector& tigs,
     vec< triple< kmer<60>,int,int> >& kmers_plus );

template void MakeKmerLookup0( const vecbasevector& tigs,
     vec< triple< kmer<80>,int,int> >& kmers_plus );

template void MakeKmerLookup0( const vecbasevector& tigs,
     vec< triple< kmer<88>,int,int> >& kmers_plus );

template void MakeKmerLookup0( const vecbasevector& tigs,
     vec< triple< kmer<96>,int,int> >& kmers_plus );

template void MakeKmerLookup0( const vecbasevector& tigs,
     vec< triple< kmer<100>,int,int> >& kmers_plus );

template void MakeKmerLookup0( const vecbasevector& tigs,
     vec< triple< kmer<128>,int,int> >& kmers_plus );

template void MakeKmerLookup0( const vecbasevector& tigs,
     vec< triple< kmer<200>,int,int> >& kmers_plus );

template void MakeKmerLookup0( const vecbasevector& tigs,
     vec< triple< kmer<320>,int,int> >& kmers_plus );

template void MakeKmerLookup0( const vecbasevector& tigs,
     vec< triple< kmer<400>,int,int> >& kmers_plus );


