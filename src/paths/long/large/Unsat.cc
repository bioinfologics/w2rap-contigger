///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Find and cluster unsatisfied links.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "CoreTools.h"
#include "Equiv.h"
#include "ParallelVecUtilities.h"
#include "ParseSet.h"
#include "VecUtilities.h"
#include "paths/HyperBasevector.h"
#include "paths/long/MakeKmerStuff.h"
#include "paths/long/ReadPath.h"
#include "paths/long/large/Unsat.h"
#include "system/SortInPlace.h"

vec<int> Nhood( const HyperBasevector& hb, const vec<int>& to_left,
                const vec<int>& to_right, const int e, const int radius ) {
    vec<int> x = {e};
    for ( int r = 0; r < radius; r++ ) {
        vec<int> x2 = x;
        for ( int l = 0; l < x.isize( ); l++ ) {
            int w = to_right[ x[l] ];
            for ( int j = 0; j < hb.From(w).isize( ); j++ )
                x2.push_back( hb.IFrom( w, j ) );
        }
        x = x2;
        for ( int l = 0; l < x.isize( ); l++ ) {
            int w = to_left[ x[l] ];
            for ( int j = 0; j < hb.To(w).isize( ); j++ )
                x2.push_back( hb.ITo( w, j ) );
        }
        x = x2;
    }
    UniqueSort(x);
    return x;
}

void MergeClusters( const vec< vec< std::pair<int,int> > >& x,
                    vec< vec< std::pair<int,int> > >& y, const vec< vec<int> >& n, const int N ) {
    vec< vec<int> > ind1(N), ind2(N);
    for ( int i = 0; i < x.isize( ); i++ )
        for ( int j = 0; j < x[i].isize( ); j++ ) {
            ind1[ x[i][j].first ].push_back(i);
            ind2[ x[i][j].second ].push_back(i);
        }
    #pragma omp parallel for
    for ( int i = 0; i < N; i++ ) {
        UniqueSort( ind1[i] ), UniqueSort( ind2[i] );
    }
    equiv_rel e( x.size( ) );
    #pragma omp parallel for
    for ( int i = 0; i < x.isize( ); i++ ) {
        vec<int> s1, s2, t1, t2;
        for ( int j = 0; j < x[i].isize( ); j++ ) {
            s1.push_back( x[i][j].first );
            s2.push_back( x[i][j].second );
        }
        UniqueSort(s1), UniqueSort(s2);

        vec<int> ss1, ss2;
        for ( int j = 0; j < s1.isize( ); j++ )
            ss1.append( n[ s1[j] ] );
        for ( int j = 0; j < s2.isize( ); j++ )
            ss2.append( n[ s2[j] ] );
        UniqueSort(ss1), UniqueSort(ss2);
        s1 = ss1;
        s2 = ss2;

        for ( int j = 0; j < s1.isize( ); j++ )
            t1.append( ind1[ s1[j] ] );
        for ( int j = 0; j < s2.isize( ); j++ )
            t2.append( ind2[ s2[j] ] );
        UniqueSort(t1), UniqueSort(t2);
        vec<int> t = Intersection( t1, t2 );
        #pragma omp critical
        {
            for ( int j = 1; j < t.isize( ); j++ )
                e.Join( t[0], t[j] );
        }
    }
    vec< vec< std::pair<int,int> > > z;
    vec<int> reps;
    e.OrbitReps(reps);
    for ( int j = 0; j < reps.isize( ); j++ ) {
        vec<int> o;
        e.Orbit( reps[j], o );
        vec< std::pair<int,int> > m;
        for ( int l = 0; l < o.isize( ); l++ )
            m.append( x[ o[l] ] );
        UniqueSort(m);
        z.push_back(m);
    }
    sortInPlaceParallel(z.begin(),z.end());
    y = z;
}

void PrintClusters( const vec< vec< std::pair<int,int> > >& xs,
                    std::map< std::pair<int,int>, int >& mult, const String& txt ) {
    Ofstream( out, txt );
    for ( int i = 0; i < xs.isize( ); i++ ) {
        vec< std::pair<int,int> > d = xs[i];
        out << "\n[" << i << "]\n";
        Sort(d);
        vec<int> all;
        for ( int j = 0; j < d.isize( ); j++ ) {
            int k = d.NextDiff(j);
            int e1 = d[j].first, e2 = d[j].second;
            all.push_back( e1, e2 );
            out << e1 << "," << e2 << " [" << mult[ d[j] ] << "]" << std::endl;
            j = k - 1;
        }
        UniqueSort(all);
        out << printSeq(all) << std::endl;
    }
}

void Unsat( const HyperBasevector& hb, const vec<int>& inv,
            const ReadPathVec& paths, vec< vec< std::pair<int,int> > >& xs,
            const String& work_dir, const int A2V ) {
    // Heuristics.

    const int max_depth = 15;
    const int max_verts = 50;
    const int radius = 3;
    const int merge_passes = 10;

    // Verbosity control.

    int verbosity = 1;

    // Set up data structures.

    vec<int> to_left, to_right;
    hb.ToLeft(to_left), hb.ToRight(to_right);

    // Phase 1.  Find unsatisfied links.

    std::cout << Date( ) << ": finding unsatisfieds" << std::endl;
    vec<vec< std::pair<int,int64_t> >> unsats( hb.EdgeObjectCount( ) );
    vec<Bool> u( paths.size( ) / 2, False );
    #pragma omp parallel for
    for ( int64_t i = 0; i < (int64_t) paths.size( ); i += 2 ) {
        const ReadPath &p1 = paths[i], &p2 = paths[i+1];
        if ( p1.size( ) == 0 || p2.size( ) == 0 ) continue;
        vec<int> x1, x2;
        for ( int i = 0; i < (int) p1.size( ); i++ )
            x1.push_back( p1[i] );
        for ( int i = ( (int) p2.size( ) ) - 1; i >= 0; i-- )
            x2.push_back( inv[ p2[i] ] );
        if ( Meet2( x1, x2 ) ) continue;
        int v = to_right[ x1.back( ) ], w = to_left[ x2.front( ) ];
        if ( v == w ) continue;
        Bool sat = False;
        vec<int> s = {v};
        for ( int d = 1; d <= max_depth; d++ ) {
            vec<int> s2;
            for ( int l = 0; l < s.isize( ); l++ ) {
                const int x = s[l];
                for ( int j = 0; j < hb.From(x).isize( ); j++ ) {
                    int y = hb.From(x)[j];
                    if ( y == w ) {
                        sat = True;
                        break;
                    } else s2.push_back(y);
                }
                if (sat) break;
            }
            if (sat) break;
            if ( s2.isize( ) > max_verts ) break;
            s = s2;
        }
        if (sat) continue;
        u[i/2] = True;
    }
    for ( int64_t i = 0; i < (int64_t) paths.size( ); i += 2 ) {
        if ( !u[i/2] ) continue;
        const ReadPath &p1 = paths[i], &p2 = paths[i+1];
        if ( p1.back( ) == p2.back( ) ) continue;
        unsats[ p1.back( ) ].push( inv[ p2.back( ) ], i/2 );
        unsats[ p2.back( ) ].push( inv[ p1.back( ) ], i/2 );
    }
    #pragma omp parallel for
    for ( int e = 0; e < hb.EdgeObjectCount( ); e++ )
        Sort( unsats[e] );

    // Create link multiplicity std::map.

    std::cout << Date( ) << ": creating multiplicity std::map" << std::endl;
    std::map< std::pair<int,int>, int > mult;
    for ( int e = 0; e < unsats.isize( ); e++ ) {
        UniqueSort( unsats[e] );
        for ( int i = 0; i < unsats[e].isize( ); i++ ) {
            int j;
            for ( j = i + 1; j < unsats[e].isize( ); j++ )
                if ( unsats[e][j].first != unsats[e][i].first ) break;
            mult[ std::make_pair( e, unsats[e][i].first ) ] = j - i;
            i = j - 1;
        }
    }

    // Delete duplicate links.

    std::cout << Date( ) << ": economizing links" << std::endl;
    #pragma omp parallel for
    for ( int e = 0; e < unsats.isize( ); e++ ) {
        vec<Bool> to_delete( unsats[e].size( ), False );
        for ( int j = 1; j < unsats[e].isize( ); j++ ) {
            if ( unsats[e][j] == unsats[e][j-1] ) {
                to_delete[j] = True;
            }
        }
        EraseIf( unsats[e], to_delete );
    }

    // Form neighborhoods.

    std::cout << Date( ) << ": forming neighborhoods" << std::endl;
    vec<vec<int>> n( hb.EdgeObjectCount( ) );
    #pragma omp parallel for
    for ( int e = 0; e < hb.EdgeObjectCount( ); e++ )
        n[e] = Nhood( hb, to_left, to_right, e, radius );

    // Form initial clusters.

    std::cout << Date( ) << ": forming initial clusters" << std::endl;
    xs.clear( );
    #pragma omp parallel for
    for ( int id1 = 0; id1 < hb.EdgeObjectCount( ); id1++ ) {
        for ( int m = 0; m < unsats[id1].isize( ); m++ ) {
            int id2 = unsats[id1][m].first;
            if ( m > 0 && unsats[id1][m-1].first == id2 ) continue;
            vec<int> id = { id1, id2 };
            vec< std::pair<int,int> > x;
            for ( int i = 0; i < n[ id[0] ].isize( ); i++ ) {
                int e1 = n[ id[0] ][i];
                for ( int j = 0; j < unsats[e1].isize( ); j++ ) {
                    int e2 = unsats[e1][j].first;
                    if ( BinMember( n[ id[1] ], unsats[e1][j].first ) )
                        x.push( e1, e2 );
                }
            }
            Sort(x);
            #pragma omp critical
            {
                xs.push_back(x);
            }
        }
    }
    std::cout << Date( ) << ": start sort" << std::endl; // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    double clock = WallClockTime( ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    sortInPlaceParallel(xs.begin(),xs.end());
    std::cout << TimeSince(clock) << " used sorting" << std::endl; // XXXXXXXXXXXXXXXXXXXXXXX
    Unique(xs);

    // Merge clusters.

    std::cout << Date( ) << ": merging clusters" << std::endl;
    double mclock = WallClockTime( );
    PRINT( xs.size( ) );
    for ( int p = 1; p <= merge_passes; p++ ) {
        MergeClusters( xs, xs, n, hb.EdgeObjectCount( ) );
        int n = 0;
        for ( int i = 0; i < xs.isize( ); i++ )
            if ( xs[i].size( ) > 1 || mult[ xs[i][0] ] > 1 ) n++;
        // DPRINT2( xs.size( ), n );
    }
    std::cout << TimeSince(mclock) << " used merging" << std::endl;

    // Remove giant clusters.

    const int max_cluster = 20;
    vec<Bool> delm( xs.size( ), False );
    #pragma omp parallel for
    for ( int i = 0; i < xs.isize( ); i++ ) {
        vec<int> m;
        for ( int j = 0; j < xs[i].isize( ); j++ )
            m.push_back( xs[i][j].first, xs[i][j].second );
        UniqueSort(m);
        if ( m.isize( ) > max_cluster ) delm[i] = True;
    }
    EraseIf( xs, delm );

    // Remove singleton clusters.

    vec<Bool> xdel( xs.size( ), False );
    for ( int i = 0; i < xs.isize( ); i++ ) {
        vec< std::pair<int,int> > d = xs[i];
        if ( d.solo( ) && mult[ d[0] ] == 1 ) xdel[i] = True;
    }
    EraseIf( xs, xdel );
    PrintClusters( xs, mult, work_dir + "/clusters.txt.ini" );

    // Look for cluster merges based on sequence overlaps.

    PRINT( xs.size( ) );
    for ( int opass = 1; opass <= 2; opass++ ) {
        std::cout << Date( ) << ": start overlap-based merging" << std::endl;
        int N = hb.EdgeObjectCount( );
        vec< vec<int> > ind1(N), ind2(N);
        for ( int i = 0; i < xs.isize( ); i++ )
            for ( int j = 0; j < xs[i].isize( ); j++ ) {
                ind1[ xs[i][j].first ].push_back(i);
                ind2[ xs[i][j].second ].push_back(i);
            }
        #pragma omp parallel for
        for ( int i = 0; i < N; i++ ) {
            UniqueSort( ind1[i] ), UniqueSort( ind2[i] );
        }
        vec<vec<std::pair<int,int>>> xs2(xs);
        #pragma omp parallel for
        for ( int i = 0; i < xs.isize( ); i++ ) {
            vec<int> s, m, r;

            // Let s be the set of all right hand sides of the cluster.

            for ( int j = 0; j < xs[i].isize( ); j++ )
                s.push_back( xs[i][j].second );
            UniqueSort(s);

            // Let m be the set of all clusters that left-share with xs[i].

            for ( int j = 0; j < xs[i].isize( ); j++ )
                m.append( ind1[ xs[i][j].first ] );
            UniqueSort(m);

            // Let r be the set of all right hand sides in m, excluding those in s.

            for ( int l = 0; l < m.isize( ); l++ )
                for ( int j = 0; j < xs[m[l]].isize( ); j++ )
                    r.push_back( xs[m[l]][j].second );
            UniqueSort(r);
            vec<Bool> rdel( r.size( ), False );
            for ( int j = 0; j < r.isize( ); j++ )
                if ( BinMember( s, r[j] ) ) rdel[j] = True;
            EraseIf( r, rdel );

            // Look for overlaps.

            const int maxo = 5;
            if ( r.isize( ) > maxo ) continue;
            vecbasevector all;
            all.reserve( s.size( ) + r.size( ) );
            for ( int j = 0; j < s.isize( ); j++ )
                all.push_back( hb.EdgeObject( s[j] ) );
            for ( int j = 0; j < r.isize( ); j++ )
                all.push_back( hb.EdgeObject( r[j] ) );
            const int L = 100;
            vec< triple<kmer<L>,int,int> > kmers_plus;
            MakeKmerLookup3( all, kmers_plus );
            vec<Bool> touched( r.size( ), False );
            for ( int j = 0; j < (int) kmers_plus.size( ); j++ ) {
                int k;
                for ( k = j + 1; k < (int) kmers_plus.size( ); k++ )
                    if ( kmers_plus[k].first != kmers_plus[j].first ) break;
                int m;
                for ( m = j; m < k; m++ )
                    if ( kmers_plus[m].second >= s.isize( ) ) break;
                if ( m > j ) {
                    for ( int l = m; l < k; l++ )
                        touched[ kmers_plus[l].second - s.isize( ) ] = True;
                }
                j = k - 1;
            }

            // Enlarge xs[i].

            vec<Bool> add( m.size( ), False );
            for ( int j = 0; j < m.isize( ); j++ ) {
                for ( int l = 0; l < xs[m[j]].isize( ); l++ ) {
                    int p = BinPosition( r, xs[m[j]][l].second );
                    if ( p >= 0 && touched[p] ) {
                        add[j] = True;
                        break;
                    }
                }
            }
            for ( int j = 0; j < add.isize( ); j++ )
                if ( add[j] ) xs2[i].append( xs[ m[j] ] );
            UniqueSort( xs2[i] );
        }
        xs = xs2;
        MergeClusters( xs, xs, n, hb.EdgeObjectCount( ) );
    }

    // Partially symmetrize.

    int nxs = xs.size( );
    for ( int i = 0; i < nxs; i++ ) {
        vec< std::pair<int,int> > d = xs[i], rd;
        for ( int j = 0; j < d.isize( ); j++ )
            rd.push( std::make_pair( inv[ d[j].second ], inv[ d[j].first ] ) );
        xs.push_back(rd);
    }
    MergeClusters( xs, xs, n, hb.EdgeObjectCount( ) );

    // Clean clusters.

    const int cluster_ratio = 10;
    for ( int i = 0; i < xs.isize( ); i++ ) {
        vec< std::pair<int,int> >& d = xs[i];
        vec<int> m( d.size( ) );
        for ( int j = 0; j < d.isize( ); j++ )
            m[j] = mult[ d[j] ];
        ReverseSortSync( m, d );
        for ( int j = 1; j < m.isize( ); j++ ) {
            if ( m[0] >= 1 && m[0] >= cluster_ratio * m[j] ) {
                d.resize(j);
                break;
            }
        }
    }

    // Print clusters.

    BinaryWriter::writeFile( work_dir + "/clusters.bin", xs );
    PrintClusters( xs, mult, work_dir + "/clusters.txt" );
}
