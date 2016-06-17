///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include <sys/types.h>

#include <algorithm>

#include "Basevector.h"
#include "CoreTools.h"
#include "math/Functions.h"
#include "pairwise_aligners/GenAlignments.h"
#include "kmers/GetNextKmerPair.h"
#include "Kclock.h"
#include "kmers/KmerRecord.h"
#include "pairwise_aligners/MakeAligns.h"
#include "pairwise_aligners/MutmerGraph.h"
#include "PackAlign.h"
#include "pairwise_aligners/ProcessFrequentKmers.h"
#include "ShortVector.h"
#include "kmers/SortKmers.h"

// The following documentation is out of date.
//
// MakeAligns is invoked by:
//
//     MakeAligns<I>( Passes, total_passes, EE, which_to_compare, max_clique,
//                    max_badness, aligns_file, log, max_errs, max_alignments )
//
// where:
//
// I = 1 or 2,
//
// Passes is either 1 or 10 or 100,
//
// EE is a vector of basevectors,
//
// See MakeAligns.h for which_to_compare.
//
// Use I = 1 if EE[i] < 1024 for all i.  Otherwise, use I = 2.
//
// Don't increase max_alignments past 255 if you plan to store them in a shortvector.


template<int I, int k, int BLOCKS_PER_NODE> void MakeAligns( int Passes,
        int total_passes, const vecbasevector& EE,
        to_compare which_to_compare, int max_clique, int max_badness,
        String aligns_file, ostream& log, int max_errs, int max_alignments,
        int min_mutmer, int local_max_errs, int stretch, int end_stretch,
        Bool dump_multiplicities, int nstretch, int local_max_errs_done,
        Bool avoid_promiscuous_kmers, int cl, Bool alt_method, int bandwidth,
        vec< vec< vec<int> > >* allowed_offsets, int max_offset_discrep,
        Bool process_frequent_kmers, String run_dir, int offset_A, int offset_B )

{
    makealigns_orig_method *orig_method_ptr = new makealigns_orig_method;
    makealigns_alt_method  *static_alt_method_ptr  = new makealigns_alt_method;
    makealigns_method *method_ptr = 0;

    if ( alt_method ) {
        static_alt_method_ptr->SetMaxErrs( max_errs );
        static_alt_method_ptr->SetBandwidth( bandwidth );

        method_ptr = static_alt_method_ptr;
    } else {
        orig_method_ptr->SetMaxBadness( max_badness );
        orig_method_ptr->SetMaxErrs( max_errs );
        orig_method_ptr->SetLocalMaxErrs( local_max_errs );
        orig_method_ptr->SetLocalMaxErrsDone( local_max_errs_done );
        orig_method_ptr->SetStretch( stretch );
        orig_method_ptr->SetEndStretch( end_stretch );
        orig_method_ptr->SetNStretch( nstretch );
        orig_method_ptr->SetCl( cl );
        orig_method_ptr->SetMaxAlignCtorCalls( 1000000 );

        method_ptr = orig_method_ptr;
    }

    MakeAligns<I,k,BLOCKS_PER_NODE>( Passes, total_passes, EE, which_to_compare,
                                     method_ptr, max_clique,
                                     aligns_file, log, max_alignments, min_mutmer,
                                     dump_multiplicities, avoid_promiscuous_kmers,
                                     allowed_offsets, max_offset_discrep,
                                     process_frequent_kmers, run_dir, offset_A, offset_B );
}


template<int I, int k, int BLOCKS_PER_NODE> void MakeAligns( int Passes,
        int total_passes, const vecbasevector& EE, to_compare which_to_compare,
        makealigns_method *method, int max_clique, String aligns_file, ostream& log,
        int max_alignments, int min_mutmer, Bool dump_multiplicities, Bool avoid_promiscuous_kmers,
        vec< vec< vec<int> > >* allowed_offsets, int max_offset_discrep,
        Bool process_frequent_kmers, String run_dir, int offset_A, int offset_B ) {
    int N = EE.size( );
    int N0 = which_to_compare.N0;
    ForceAssert( N0 > 0 || which_to_compare.compare_type == ALL_VS_ALL );

    kclock pass1;

    if ( I == 1 ) {
        for ( int i = 0; i < N; i++ )
            if ( EE[i].size( ) >= 1024 )
                FatalErr( "I'm sorry, but MakeAligns<1> has been passed a read "
                          << "of length " << EE[i].size( ) << ".  For reads of "
                          << "length >= 1024, you need to use MakeAligns<2>." );
    }

    mutmer_graph<I, BLOCKS_PER_NODE> M(N);

    vec<int> rid(N);
    for ( int i = 0; i < N; i++ )
        rid[i] = i;
    longlong S_init = 0;
    for ( size_t l = 0; l < EE.size( ); l++ )
        if ( EE[l].size( ) >= k ) S_init += EE[l].size( ) - k;
    S_init += S_init/4;

    Assert( Passes == 1 || Passes == 10 || Passes == 100 );
    if ( Passes == 10 ) S_init /= 7;
    if ( Passes == 100 ) S_init /= 33;

    if ( S_init > 3000000000u )
        FatalErr( "MakeAligns: the value of S is " << S_init
                  << ".  This is dangerously large, because it represents\nthe "
                  << "size of a vector which is to be stored in an unsigned int,\n"
                  << "and this vector could grow. "
                  << "You may need to rewrite some code to get it to work." );

    unsigned int S = S_init;

    unsigned int max_S = 0, total_S = 0;

    {
        vec< kmer_record<k,I> > R(S);
        vec<int> multiplicities;
        log << "pass " << flush;
        for ( int pass = 0; pass < total_passes; pass++ ) {
            dummy<1> d1;
            dummy<10> d10;
            dummy<100> d100;
//	       cerr << "pass #" << pass << " of " << total_passes << std::endl;
            if ( Passes == 1 ) SortKmers( d1, EE, rid, pass, R, S );
            else if ( Passes == 10 ) SortKmers( d10, EE, rid, pass, R, S );
            else SortKmers( d100, EE, rid, pass, R, S );
            Assert( R.size( ) >= S );
            if (process_frequent_kmers)
                ProcessFrequentKmers( EE, R, S, max_clique, run_dir );
            max_S = max(S, max_S);
            total_S += S;
            // gcc 4.2 warns pos1 may be used uninit'd, but in reality, this
            // will never happen.  So, we init to 0 to suppress the warn.
            int read_id1(-1), read_id2, pos1=0, pos2;
            while(1) {
                if ( which_to_compare.compare_type != FIRST_VS_SECOND )
                    GetNextKmerPair<k,I>( R, S, read_id1, read_id2, pos1, pos2,
                                          max_clique, multiplicities, dump_multiplicities,
                                          avoid_promiscuous_kmers );
                else GetNextKmerPair<k,I>( N0, R, S, read_id1, read_id2, pos1,
                                               pos2, max_clique, multiplicities, dump_multiplicities,
                                               avoid_promiscuous_kmers );
                if ( read_id1 < 0 ) break;
                if (dump_multiplicities) continue;

//    		    cerr << "got next kmer pair"
// 			 << " id1: " << read_id1 << " id2: " << read_id2
// 			 << std::endl;


                if ( N0 > 0 ) {
                    if ( which_to_compare.compare_type == FIRST_VS_ALL
                            && read_id1 >= N0 && read_id2 >= N0 ) continue;

                    if ( allowed_offsets ) {
                        int p1 = pos1;
                        int p2 = pos2;

                        int rc = (p1 < 0) ^ (p2 < 0);
                        if ( p1 < 0 ) p1 = -p1;
                        if ( p2 < 0 ) p2 = -p2;
                        --p1;
                        --p2;

                        if ( read_id1 < N0 ) {
                            if (rc) p1 = EE[read_id1].size( ) - k - p1;
                        } else {
                            swap( p1, p2 );
                            if (rc) p1 = EE[read_id2].size( ) - k - p1;
                        }

                        int offset = p1 - p2;
                        vec<int>& v = (read_id1 < N0)
                                      ? (*allowed_offsets)[read_id1][read_id2 - N0]
                                      : (*allowed_offsets)[read_id2][read_id1 - N0];
                        unsigned int i;
                        for ( i = 0; i < v.size( ); i++ )
                            if ( Abs(offset - v[i]) <= max_offset_discrep )
                                break;
                        if ( i == v.size( ) ) continue;
                    }
                }

                M.MergeKmer( pos1, pos2, k, read_id1, read_id2,
                             EE[read_id1], EE[read_id2] );
            }
            if ( pass == 97 ) log << "1";
            else if ( pass % 10 == 8 ) log << (pass/10 + 1) % 10;
            else if ( pass % 10 == 9 ) log << "0 ";
            else log << ".";
            flush(log);
        }
        log << "\n" << flush;
        if (dump_multiplicities) {
            sort( multiplicities.begin( ), multiplicities.end( ) );
            Ofstream( mult, aligns_file + ".mult" );
            for ( unsigned int i = 0; i < multiplicities.size( ); i++ )
                mult << multiplicities[i] << "\n";
            return;
        }
    }

    log << pass1.UserTime( ) << " seconds used in phase 1\n";
    log << "largest vector to sort has " << max_S << " entries" << std::endl;

    if ( total_S == 0 ) {
        cout << "MakeAligns: I've been passed reads of length";
        for ( size_t l = 0; l < EE.size( ); l++ )
            cout << " " << EE[l].size( );
        cout << "\nI can't work with this data.\n";
        ForceAssert( 0 == 1 );
    }

    log << "This is " << 100.0 * double(max_S) / double(total_S) <<
        "% of what would be needed for a single pass.\n\n";

    kclock pass2;

    // Generate file containing alignments:
    GenAlignments<I,k,BLOCKS_PER_NODE>( EE, M, N0, aligns_file, &log,
                                        min_mutmer, max_alignments, method,
                                        offset_A, offset_B);

    log << pass2.UserTime( ) << " seconds used in phase 2\n";
}

#define INSTANTIATE_MAKEALIGNS(I, K, BLOCKS_PER_NODE)  \
template void MakeAligns<I, K, BLOCKS_PER_NODE>( int Passes, \
     int total_passes, const vecbasevector& EE, to_compare which_to_compare, \
     makealigns_method *method, int max_clique, String aligns_file, ostream& log, \
     int max_alignments, int min_mutmer, Bool dump_multiplicities, Bool avoid_promiscuous_kmers, \
     vec< vec< vec<int> > > *allowed_offsets, int max_offset_discrep, \
     Bool process_frequent_kmers, String run_dir, int offset_A, int offset_B ); \
\
template void MakeAligns<I, K, BLOCKS_PER_NODE>( int Passes, \
     int total_passes, const vecbasevector& EE, \
     to_compare which_to_compare, int max_clique, int max_badness, \
     String aligns_file, ostream& log, int max_errs, int max_alignments, \
     int min_mutmer, int local_max_errs, int stretch, int end_stretch, \
     Bool dump_multiplicities, int nstretch, int local_max_errs_done, \
     Bool avoid_promiscuous_kmers, int cl, Bool alt_method, int bandwidth,\
     vec< vec< vec<int> > > *allowed_offsets, int max_offset_discrep, \
     Bool process_frequent_kmers, String run_dir, int offset_A, int offset_B ) \

#define INSTANTIATE_MAKEALIGNS_FOR_K(K,dummy) \
   INSTANTIATE_MAKEALIGNS(1,K,12); \
   INSTANTIATE_MAKEALIGNS(2,K,12); \
   INSTANTIATE_MAKEALIGNS(1,K,50); \
   INSTANTIATE_MAKEALIGNS(2,K,50)

FOR_ALL_K(INSTANTIATE_MAKEALIGNS_FOR_K,unused);

