///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

// SeqInfo. Find edges containing a particular sequence in the forward direction.
// Example 1.  Inside a.fin,
//             SeqInfo S=ACGGTATGATTAGCTATAC
// Example 2.  Inside assembly dir,
//             SeqInfo S=a.200:5046 DIR=a.200b

#include "Basevector.h"
#include "MainTools.h"
#include "paths/HyperBasevector.h"

int main(int argc, char *argv[]) {
    RunTime( );

    BeginCommandArguments;
    CommandArgument_String_OrDefault_Doc(DIR, ".", "looks for DIR/a.{fastb}");
    CommandArgument_String_Doc(S, "either a DNA sequence, or an expression of "
                               "the form dir:n, to find sequence n in dir/a.fastb");
    CommandArgument_Int_OrDefault_Doc(TRUNC, -1,
                                      "truncate sequence to the given length");
    EndCommandArguments;

    // Define query sequence.

    String x;
    if ( !S.Contains( ":" ) ) {
        for ( int i = 0; i < S.isize( ); i++ ) {
            if ( S[i] != 'A' && S[i] != 'C' && S[i] != 'G' && S[i] != 'T' ) {
                cout << "The character '" << S[i] << "' in S is not a DNA base."
                     << std::endl;
                Scram(1);
            }
        }
        x = S;
    } else {
        vecbasevector t;
        String dir = S.Before( ":" );
        if ( !IsDirectory(dir) ) {
            cout << "Can't find directory " << dir << "." << std::endl;
            Scram(1);
        }
        String fastb = dir + "/a.fastb";
        if ( !IsRegularFile(fastb) ) {
            cout << "Can't find file " << fastb << std::endl;
            Scram(1);
        }
        String id = S.After( ":" );
        if ( !id.IsInt( ) ) {
            cout << "Illegal id " << id << "." << std::endl;
            Scram(1);
        }
        t.ReadOne( fastb, id.Int( ) );
        x = t[0].ToString( );
    }
    if ( TRUNC >= 0 ) {
        if ( TRUNC > x.isize( ) ) {
            cout << "That truncation would make the sequence longer." << std::endl;
            Scram(1);
        }
        x.resize(TRUNC);
    }

    // Look for it.

    HyperBasevectorX hb;
    BinaryReader::readFile( DIR + "/a.hbx", &hb );
    int K = hb.K( );
    vecbasevector tigs( DIR + "/a.fastb" );
    if ( x.isize( ) >= K ) {
        vec<vec<int>> hits;
        String y = x.substr( 0, K );
        #pragma omp parallel for
        for ( int e = 0; e < (int ) tigs.size( ); e++ ) {
            String s = tigs[e].ToString( );
            int p = s.Position(y);
            if ( p < 0 ) continue; // misses cases where it occurs > once
            vec<int> E = {e};
            int epos = p;
            Bool ok = True;
            for ( int j = 0; j < x.isize( ); j++ ) {
                if ( epos == hb.EdgeObject( E.back( ) ).isize( ) ) {
                    int v = hb.ToRight( E.back( ) );
                    vec<int> exts;
                    for ( int l = 0; l < (int) hb.From(v).size( ); l++ ) {
                        if ( hb.EdgeObject( hb.IFrom( v, l ) )[K-1]
                                == as_char(x[j]) ) {
                            exts.push_back( hb.IFrom( v, l ) );
                        }
                    }
                    if ( !exts.solo( ) ) {
                        ok = False;
                        break;
                    }
                    E.push_back( exts[0] );
                    epos = K;
                } else if ( hb.EdgeObject( E.back( ) )[epos] != as_char(x[j]) ) {
                    ok = False;
                    break;
                } else epos++;
            }
            if (ok) hits.push_back(E);
        }
        Sort(hits);
        for ( int i = 0; i < hits.isize( ); i++ ) {
            if ( i > 0 ) std::cout << ";";
            cout << printSeq( hits[i] );
        }
        cout << std::endl;
    } else {
        vec<int> hits;
        #pragma omp parallel for
        for ( int e = 0; e < (int ) tigs.size( ); e++ ) {
            String s = tigs[e].ToString( );
            if ( s.Contains(x) ) {
                #pragma omp critical
                {    hits.push_back(e);    }
            }
        }
        Sort(hits);
        cout << printSeq( hits.begin( ), hits.end( ), ";" ) << std::endl;
    }
}
