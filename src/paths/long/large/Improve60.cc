///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "CoreTools.h"
#include "ParallelVecUtilities.h"
#include "Qualvector.h"
#include "feudal/PQVec.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "paths/long/large/GapToyTools.h"
#include "paths/long/large/Improve60.h"

void Improve60( HyperBasevector& hb, vec<int>& inv, ReadPathVec& paths,
                const vecbasevector& bases, const VecPQVec& qualsp, const Bool verbose ) {
    // Two passes.

    for ( int xpass = 1; xpass <= 2; xpass++ ) {
        if (verbose) std::cout << "\nstarting pass " << xpass << "\n\n";


        // Create indices.

        vec<int> to_left, to_right;
        hb.ToLeft(to_left), hb.ToRight(to_right);
        VecULongVec paths_index;
        invert( paths, paths_index, hb.EdgeObjectCount( ) );

        // Find an edge e: v --> w that looks like junk, and for which w is a dead
        // end with only e entering.  Find an alternative path f = f1,...,fn
        // starting at e, that is at least as long as e.  For the columns at which e
        // and f differ, compute the quality sums, etc.  Starting with n = 1.  And of
        // course we'll do the same thing the other way.
        //
        // Starting with special case.

        vec<int> dels;
        #pragma omp parallel for
        for ( int e = 0; e < hb.E( ); e++ ) {
            int v = to_left[e], w = to_right[e];
            if ( hb.From(v).size( ) != 2 ) continue;
            int f;
            if ( hb.IFrom(v,0) != e ) f = hb.IFrom(v,0);
            else f = hb.IFrom(v,1);
            if ( hb.Kmers(e) > hb.Kmers(f) ) continue;

            Bool OK = False;
            if ( hb.From(w).empty( ) && hb.To(w).solo( ) ) OK = True;
            // Test for simple bubble.
            if ( hb.From(v).size( ) == 2 && hb.To(w).size( ) == 2 && v != w
                    && hb.To(v).solo( ) && hb.From(w).solo( )
                    && hb.From(v)[0] == w && hb.From(v)[1] == w ) {
                OK = True;
            }
            if ( !OK ) continue;

            // if ( hb.From(w).nonempty( ) || hb.To(w).size( ) != 1 ) continue;

            vec<int> diffs;
            for ( int j = hb.K( ) - 1; j < hb.Bases(e); j++ )
                if ( hb.O(e)[j] != hb.O(f)[j] ) diffs.push_back(j);
            vec<int> delta;
            for ( int pass = 1; pass <= 2; pass++ ) {
                int s = ( pass == 1 ? e : f );
                for ( int rpass = 1; rpass <= 2; rpass++ ) {
                    int t = ( rpass == 1 ? s : inv[s] );
                    for ( int64_t i = 0; i < (int64_t) paths_index[t].size( ); i++ ) {
                        int64_t id = paths_index[t][i];
                        const ReadPath& p = paths[id];
                        int start = p.getOffset( );
                        basevector b = bases[id];
                        qualvector q;
                        qualsp[id].unpack(&q);
                        if ( rpass == 2 ) {
                            b.ReverseComplement( );
                            q.ReverseMe( );
                            start = hb.Bases(f) - start - b.isize( );
                        }
                        for ( int j = 0; j < (int) p.size( ); j++ ) {
                            if ( p[j] == t ) {
                                int qsum1 = 0, qsum2 = 0;
                                for ( int d = 0; d < diffs.isize( ); d++ ) {
                                    int l = diffs[d];
                                    if ( l-start < 0 || l-start >= b.isize( ) )
                                        continue;
                                    if ( hb.O(e)[l] != b[l-start] )
                                        qsum1 += q[l-start];
                                    if ( hb.O(f)[l] != b[l-start] )
                                        qsum2 += q[l-start];
                                }
                                if ( qsum1 != qsum2 )
                                    delta.push_back( qsum1 - qsum2 );
                            }
                            if ( rpass == 1 ) start -= hb.Kmers( p[j] );
                            else start += hb.Kmers( p[j] );
                        }
                    }
                }
            }
            ReverseSort(delta);
            if ( delta.empty( ) ) continue;
            int sum1 = 0, sum2 = 0;
            for ( int i = 0; i < delta.isize( ); i++ ) {
                if ( delta[i] > 0 && i > 0 ) sum1 += delta[i];
                else if ( delta[i] < 0 && i < delta.isize( ) - 1 ) {
                    sum2 -= delta[i];
                }
            }
            Bool del = ( sum1 >= 100 && sum2 < 100 && sum1 >= 10 * sum2 );
            if ( verbose || del ) {
                #pragma omp critical
                {
                    if (verbose) {
                        cout << "e = " << e << ", inv[e] = " << inv[e]
                        << ", delta = " << printSeq(delta);
                        if (del) std::cout << " [delete]";
                        cout << std::endl;
                    }
                    if (del) dels.push_back( e, inv[e] );
                }
            }
        }

        // Delete edges and clean up.

        ParallelUniqueSort(dels);
        hb.DeleteEdges(dels);
        Cleanup( hb, inv, paths );
    }
}
