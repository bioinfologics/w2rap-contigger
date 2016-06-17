///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// AffineAlign.  Produce visual affine alignment of two sequences.  Sadly, the two
// sequences have to match up exactly on the ends, or something bad happens.

#include "FetchReads.h"
#include "MainTools.h"
#include "PrintAlignment.h"
#include "pairwise_aligners/SmithWatAffine.h"

int main(int argc, char *argv[]) {
    RunTime( );

    BeginCommandArguments;
    CommandArgument_String_Doc(F1, "first fasta file; uses only first record");
    CommandArgument_String_Doc(F2, "second fasta file; uses only first record");
    CommandArgument_Bool_OrDefault_Doc(RC1, False,
                                       "reverse complement first entry");
    CommandArgument_Bool_OrDefault_Doc(RC2, False,
                                       "reverse complement second entry");
    CommandArgument_Int_OrDefault(MISMATCH, 3);
    CommandArgument_Int_OrDefault(GAP_OPEN, 12);
    CommandArgument_Int_OrDefault(GAP_EXTEND, 1);
    CommandArgument_Int_OrDefault_Doc(ID1, 0, "use this record on F1");
    CommandArgument_Int_OrDefault_Doc(ID2, 0, "use this record on F2");
    CommandArgument_Int_OrDefault_Doc(START2, 0,
                                      "start at this zero-based position on ID2");
    CommandArgument_Int_OrDefault_Doc(STOP2, -1,
                                      "stop at this zero-based position on ID2");
    CommandArgument_Int_OrDefault_Doc(BANDWIDTH, -1,
                                      "if set do a banded alignment using this bandwidth");
    CommandArgument_Bool_OrDefault_Doc(LIST_EVENTS, False,
                                       "list events using 1-based coordinates");
    CommandArgument_Int_OrDefault_Doc(LIST_EVENTS_OFFSET, 0,
                                      "add this to event coordinates");
    EndCommandArguments;

    vecbasevector f1, f2;
    FetchReads( f1, 0, F1 );
    FetchReads( f2, 0, F2 );

    if( RC1) {
        for( auto& entry: f1) {
            entry.ReverseComplement();
        }
    }
    if( RC2) {
        for( auto& entry: f2) {
            entry.ReverseComplement();
        }
    }

    if ( STOP2 < 0 ) STOP2 = f2[ID2].size( );
    basevector X2( f2[ID2], START2, STOP2 - START2 );

    align a;

    if ( BANDWIDTH < 0 ) {
        alignment al;
        SmithWatAffineParallel( f1[ID1], X2, al, true, true,
                                MISMATCH, GAP_OPEN, GAP_EXTEND );
        a = al;
    } else {
        int nerrors;
        int OFFSET = ( f1[ID1].isize( ) - X2.isize( ) ) / 2;
        SmithWatAffineBanded( f1[ID1], X2, OFFSET, BANDWIDTH, a, nerrors,
                              MISMATCH, GAP_OPEN, GAP_EXTEND );
    }

    if (LIST_EVENTS) {
        int LIST_EVENTS_BASE = 1;
        const basevector &query = f1[ID1], &target = X2;
        int p1 = a.pos1( ), p2 = a.pos2( );
        for ( int j = 0; j < a.Nblocks( ); j++ ) {
            int g = a.Gaps(j);
            if ( g > 0 ) {
                cout << "EVENT: deletion of " << g << " bases at ref."
                     << p2 + LIST_EVENTS_BASE + LIST_EVENTS_OFFSET << std::endl;
                p2 += g;
            }
            if ( g < 0 ) {
                String ins;
                for ( int j = 0; j < -g; j++ )
                    ins.push_back( as_base( query[ p1 + j ] ) );
                p1 -= g;
                cout << "EVENT: insertion of " << -g << " bases ("
                     << ins << ") at ref."
                     << p2 + LIST_EVENTS_BASE + LIST_EVENTS_OFFSET << std::endl;
            }
            for ( int x = 0; x < a.Lengths(j); x++ ) {
                if ( query[p1] != target[ p2] ) {
                    cout << "EVENT: substitution at ref."
                         << p2 + LIST_EVENTS_BASE + LIST_EVENTS_OFFSET
                         << ", " << as_base( target[ p2] )
                         << " changed to " << as_base( query[p1] ) << std::endl;
                }
                ++p1, ++p2;
            }
        }
    }

    PrintVisualAlignment( True, cout, f1[ID1], X2, a );
}
