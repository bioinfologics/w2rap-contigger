///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "CoreTools.h"
#include "PairsManager.h"
#include "graph/Digraph.h"
#include "math/Functions.h"
#include "paths/Uniseq.h"
#include "paths/AssemblyCleanupTools.h"

int uniseq::Len( ) const {
    Assert( unibases_ != 0 );
    int len = -Sum( Over( ) );
    for ( int j = 0; j < N( ); j++ )
        len += Unibase( U(j) ).size( );
    return len;
}

void uniseq::ReverseMe( const vec<int>& to_rc ) {
    u_.ReverseMe( );
    overlap_.ReverseMe( );
    for ( int j = 0; j < N( ); j++ )
        u_[j] = to_rc[ u_[j] ];
}

uniseq uniseq::Reverse( const vec<int>& to_rc ) const {
    uniseq x(*this);
    x.ReverseMe( to_rc );
    return x;
}

uniseq Cat( const uniseq& s1, const uniseq& s2 ) {
    vec<int> u = s1.U( );
    u.resize( u.isize( ) - 1 );
    u.append( s2.U( ) );
    vec<int> over = s1.Over( );
    over.append( s2.Over( ) );
    return uniseq( u, over );
}

uniseq Cat( const uniseq& u1, const uniseq& u2, const uniseq& u3 ) {
    return Cat( Cat( u1, u2 ), u3 );
}

std::ostream& operator<<( std::ostream& out, const uniseq& q ) {
    Assert( q.Alive( ) );
    for ( int r = 0; r < q.N( ); r++ ) {
        out << q.U(r);
        if ( r < q.N( ) - 1 ) out << " --[" << q.Over(r) << "]--> ";
    }
    return out;
}

void uniseq::Print( std::ostream& out, const int K ) const {
    Assert( Alive( ) );
    for ( int r = 0; r < N( ); r++ ) {
        out << U(r);
        if ( r < N( ) - 1 ) {
            if ( Over(r) == K-1 ) out << "+";
            else out << "--[" << Over(r) << "]-->";
        }
    }
}

const vecbasevector* uniseq::unibases_(0);

int gapster::MinLen( const double dev_mult ) const {
    if ( Closed( ) ) {
        int min_len = 1000000000;
        {
            for ( int j = 0; j < ClosureCount( ); j++ )
                min_len = Min( min_len, Closure(j).Len( ) );
            return min_len;
        }
    } else return Sep( ) + int( floor( -dev_mult * double( Dev( ) ) ) );
}

int gapster::MaxLen( const double dev_mult ) const {
    if ( Closed( ) ) {
        int max_len = 0;
        for ( int j = 0; j < ClosureCount( ); j++ )
            max_len = Max( max_len, Closure(j).Len( ) );
        return max_len;
    } else return Sep( ) + int( ceil( dev_mult * double( Dev( ) ) ) );
}

void snark::CloseGap( const int v1, const int v2, const vec<uniseq>& closures ) {
#ifndef NDEBUG
    {
        for ( size_t i = 0; i < closures.size( ); i++ ) {
            const uniseq& c = closures[i];
            Assert( Vert(v1).U( ).back( ) == c.U( ).front( ) );
            AssertEq( c.U( ).back( ), Vert(v2).U( ).front( ) );
        }
    }
#endif
    int p = BinPosition( G_.From(v1), v2 );
    AssertGe( p, 0 );
    G_.EdgeObjectByIndexFromMutable( v1, p ).Close(closures);
}

void snark::SwallowSimpleGaps( ) {
    for ( int x = 0; x < VertN( ); x++ ) {
        if ( Vert(x).Dead( ) ) continue;
        if ( !G( ).From(x).solo( ) ) continue;
        int y = G( ).From(x)[0];
        if ( y == x ) continue;
        const gapster& g = G( ).EdgeObjectByIndexFrom( x, 0 );
        if ( g.ClosureCount( ) != 1 || !G( ).To(y).solo( ) ) continue;
        seq_[x] = Cat( Vert(x), g.Closure(0), Vert(y) );
        G_.FromMutable(x) = G( ).From(y);
        G_.FromEdgeObjMutable(x) = G( ).FromEdgeObj(y);
        G_.FromMutable(y).clear( ), G_.ToMutable(y).clear( );
        G_.FromEdgeObjMutable(y).clear( ), G_.ToEdgeObjMutable(y).clear( );
        for ( int j = 0; j < G( ).From(x).isize( ); j++ ) {
            int z = G( ).From(x)[j];
            int p = BinPosition( G( ).To(z), y );
            G_.ToMutable(z)[p] = x;
            SortSync( G_.ToMutable(z), G_.ToEdgeObjMutable(z) );
        }
        seq_[y].Kill( );
        x--;
    }

    // Sometimes x --> y --> z can be reduced to x --> z.

    for ( int y = 0; y < VertN( ); y++ ) {
        if ( Vert(y).Dead( ) ) continue;
        if ( !From(y).solo( ) || !To(y).solo( ) ) continue;
        int x = To(y)[0], z = From(y)[0];
        if ( x == y || y == z || x == z ) continue;
        const gapster& g1 = G( ).EdgeObjectByIndexTo( y, 0 );
        const gapster& g2 = G( ).EdgeObjectByIndexFrom( y, 0 );
        if ( g1.ClosureCount( ) != 1 || g2.ClosureCount( ) != 1 ) continue;
        uniseq u = Cat( g1.Closure(0), Vert(y), g2.Closure(0) );
        Gmutable( ).AddEdge( x, z, u );
        Gmutable( ).DeleteEdgeTo( y, 0 );
        Gmutable( ).DeleteEdgeFrom( y, 0 );
        seq_[y].Kill( );
        x--;
    }

    // Look for parallel and identical edges.

    vec<int> to_delete;
    for ( int x = 0; x < VertN( ); x++ ) {
        for ( int j1 = 0; j1 < From(x).isize( ); j1++ )
            for ( int j2 = j1 + 1; j2 < From(x).isize( ); j2++ ) {
                if ( From(x)[j1] != From(x)[j2] ) continue;
                if ( G( ).EdgeObjectByIndexFrom( x, j1 )
                        == G( ).EdgeObjectByIndexFrom( x, j2 ) ) {
                    to_delete.push(
                        G( ).EdgeObjectIndexByIndexFrom( x, j2 ) );
                }
            }
    }
    UniqueSort(to_delete);
    Gmutable( ).DeleteEdges(to_delete);
}

const vecbasevector* snark::unibases_(0);
// const vec<int>* snark::to_rc_(0);

void FindPartnersInGap(

    // vertices to look between

    const int x1, const int x2,

    // distance to scan on left and right

    const int flank,

    // assembly

    const snark& S,

    // reads and pairs

    const vecbasevector& bases, const PairsManager& pairs,

    // 1. read --> ( unipath, pos, fw? )
    // 2. unipath --> ( read_id, pos, fw? )

    const vec< vec< triple<int,int,Bool> > >& placements_by_read,
    const vec< vec< triple<int64_t,int,Bool> > >& placements_by_unipath,

    // ( unipath, start, stop )

    vec< triple<int,int,int> >& between ) {
    int len1 = S.Vert(x1).Len( );
    int pos = 0;
    for ( int i = 0; i < S.Vert(x1).N( ); i++ ) {
        int u1 = S.Vert(x1).U(i);
        for ( int j = 0; j < placements_by_unipath[u1].isize( ); j++ ) {
            if ( !placements_by_unipath[u1][j].third ) continue;
            int start1 = pos + placements_by_unipath[u1][j].second;
            if ( len1 - start1 > flank ) continue;
            int64_t id1 = placements_by_unipath[u1][j].first;
            int64_t id2 = pairs.getPartnerID(id1);
            for ( int l = 0; l < placements_by_read[id2].isize( ); l++ ) {
                if ( placements_by_read[id2][l].third ) continue;
                int u2 = placements_by_read[id2][l].first;
                int start2 = placements_by_read[id2][l].second;
                int stop2 = start2 + bases[id2].isize( );
                between.push( u2, start2, stop2 );
            }
        }
        if ( i < S.Vert(x1).N( ) - 1 )
            pos += S.Unibase(u1).isize( ) - S.Vert(x1).Over(i);
    }
    pos = 0;
    for ( int i = 0; i < S.Vert(x2).N( ); i++ ) {
        int u2 = S.Vert(x2).U(i);
        for ( int j = 0; j < placements_by_unipath[u2].isize( ); j++ ) {
            if ( placements_by_unipath[u2][j].third ) continue;
            int start2 = pos + placements_by_unipath[u2][j].second;
            int64_t id2 = placements_by_unipath[u2][j].first;
            int stop2 = start2 + bases[id2].isize( );
            if ( stop2 > flank ) continue;
            int64_t id1 = pairs.getPartnerID(id2);
            for ( int l = 0; l < placements_by_read[id1].isize( ); l++ ) {
                if ( !placements_by_read[id1][l].third ) continue;
                int u1 = placements_by_read[id1][l].first;
                int start1 = placements_by_read[id1][l].second;
                int stop1 = start1 + bases[id1].isize( );
                between.push( u1, start1, stop1 );
            }
        }
    }
    Sort(between);
}

void GetPairPlacements(

    // assembly and ancillary data:

    const snark& S, const vec<int>& to_right,
    const vec<int>& min_len, const vec<int>& max_len,

    // input pair:

    const int64_t pid,

    // reads and pairs:

    const vecbasevector& bases, const PairsManager& pairs,

    // read --> ( unipath, pos, fw? ):

    const vec< vec< triple<int,int,Bool> > >& placements_by_read,

    // edge_id --> ( u, pos, apos ):

    const vec< vec< triple<int,int,int> > >& u_pos_apos,

    // placements of the reads -- pos_on[0] provides placements of the first read
    // and pos_on[1] provided placements of the second read:

    vec< vec<placement_on> >& pos_on,

    // placements of the pairs -- the first and second entries are indices in
    // pos_on[0] and pos_on[1], respectively; the third entry is set if the second
    // read is fw rather than the first:

    vec< triple<int,int,Bool> >& placements ) {
    // Define some basic stuff about the pair.

    int64_t id1 = pairs.ID1(pid), id2 = pairs.ID2(pid);
    int mean_sep = pairs.sep(pid), dev = pairs.sd(pid);
    int read_len1 = bases[id1].size( ), read_len2 = bases[id2].size( );
    vec<int> read_len;
    read_len.push_back( read_len1, read_len2 );

    // Find all the places that the reads could be.

    pos_on.clear_and_resize(2);
    for ( int pass = 0; pass < 2; pass++ ) {
        int64_t id = ( pass == 0 ? id1 : id2 );
        vec<int> uni1;
        for ( int j = 0; j < placements_by_read[id].isize( ); j++ )
            uni1.push_back( placements_by_read[id][j].first );
        UniqueSort(uni1);
        for ( int x = 0; x < S.VertN( ); x++ ) {
            vec< std::pair<int,Bool> > pos_on_x;
            const uniseq& r = S.Vert(x);
            int len = r.Len( ), pos = 0;
            for ( int j = 0; j < r.N( ); j++ ) {
                int u = r.U(j);
                if ( j > 0 ) pos += S.Unibase( r.U(j-1) ).isize( ) - r.Over(j-1);
                if ( !BinMember( uni1, u ) ) continue;
                for ( int l = 0; l < placements_by_read[id].isize( ); l++ ) {
                    if ( placements_by_read[id][l].first != u ) continue;
                    int xpos = pos + placements_by_read[id][l].second;
                    Bool xfw = placements_by_read[id][l].third;
                    pos_on_x.push( xpos, xfw );
                }
            }
            UniqueSort(pos_on_x);
            for ( int j = 0; j < pos_on_x.isize( ); j++ ) {
                pos_on[pass].push( x, pos_on_x[j].first,
                                   len - pos_on_x[j].first - read_len[pass],
                                   pos_on_x[j].second );
            }
        }
        for ( int gi = 0; gi < S.EdgeN( ); gi++ ) {
            const gapster& g = S.Edge(gi);
            if ( g.Open( ) ) continue;
            vec< triple<int,int,Bool> > pos_on_e;
            for ( int j = 0; j < u_pos_apos[gi].isize( ); j++ ) {
                int u = u_pos_apos[gi][j].first;
                if ( !BinMember( uni1, u ) ) continue;
                int pos = u_pos_apos[gi][j].second;
                int apos = u_pos_apos[gi][j].third;
                for ( int l = 0; l < placements_by_read[id].isize( ); l++ ) {
                    if ( placements_by_read[id][l].first != u ) continue;
                    int xpos = pos + placements_by_read[id][l].second;
                    int axpos = apos - placements_by_read[id][l].second
                                - read_len[pass];
                    Bool xfw = placements_by_read[id][l].third;
                    pos_on_e.push( xpos, axpos, xfw );
                }
            }
            UniqueSort(pos_on_e);
            for ( int j = 0; j < pos_on_e.isize( ); j++ ) {
                pos_on[pass].push( S.VertN( ) + gi, pos_on_e[j].first,
                                   pos_on_e[j].second, pos_on_e[j].third );
            }
        }
    }

    // Now look for pair placements.

    const int dev_mult = 3;
    placements.clear( );
    for ( int j1 = 0; j1 < pos_on[0].isize( ); j1++ )
        for ( int j2 = 0; j2 < pos_on[1].isize( ); j2++ ) {
            placement_on p1 = pos_on[0][j1], p2 = pos_on[1][j2];
            if ( p1.fw == p2.fw ) continue;
            int len1 = read_len1, len2 = read_len2;
            Bool swapped = False;
            if ( !p1.fw ) {
                std::swap( p1, p2 );
                std::swap( len1, len2 );
                swapped = True;
            }
            if ( !p1.fw ) continue;
            int min_sep = mean_sep - dev_mult * dev;
            int max_sep = mean_sep + dev_mult * dev;
            int sep0 = p1.apos + p2.pos;
            int p1id = p1.id, p2id = p2.id;
            int p1pos = p1.pos, p2pos = p2.pos;
            if ( p1id < S.VertN( ) && p2id < S.VertN( ) && p1id == p2id ) {
                int sep = p2pos - p1pos - len1;
                if ( Abs( sep - mean_sep ) <= dev_mult * dev ) {
                    placements.push( j1, j2, swapped );
                    continue;
                }
            }
            if ( p1id >= S.VertN( ) && p2id < S.VertN( ) ) {
                p1id = to_right[ p1id - S.VertN( ) ];
                if ( p1id == p2id ) {
                    if ( Abs( sep0 - mean_sep ) <= dev_mult * dev ) {
                        placements.push( j1, j2, swapped );
                        continue;
                    }
                }
                sep0 += S.Vert(p1id).Len( );
            }
            if ( p1id >= S.VertN( ) && p2id >= S.VertN( ) ) {
                if ( p1id == p2id ) {
                    int sep = p2pos - p1pos - len1;
                    if ( Abs( sep - mean_sep ) <= dev_mult * dev ) {
                        placements.push( j1, j2, swapped );
                        continue;
                    }
                }
                p1id = to_right[ p1id - S.VertN( ) ];
                sep0 += S.Vert(p1id).Len( );
            }
            if ( sep0 > max_sep ) continue;
            vec< triple<int,int,int> > pp;
            pp.push( p1id, sep0, sep0 );
            while( pp.nonempty( ) ) {
                int x = pp.back( ).first;
                int sep1 = pp.back( ).second, sep2 = pp.back( ).third;
                pp.pop_back( );
                for ( int m = 0; m < S.G( ).From(x).isize( ); m++ ) {
                    int y = S.G( ).From(x)[m];
                    int sep1_new(sep1), sep2_new(sep2);
                    int z = S.G( ).EdgeObjectIndexByIndexFrom( x, m );
                    const gapster& g = S.G( ).EdgeObjectByIndexFrom( x, m );
                    if ( z == p2.id - S.VertN( ) && sep2_new >= min_sep ) {
                        placements.push( j1, j2, swapped );
                        pp.clear( );
                        break;
                    }
                    sep1_new +=
                        min_len[z] - S.Unibase( S.Vert(x).U( ).back( ) ).isize( )
                        - S.Unibase( S.Vert(y).U( ).front( ) ).isize( );
                    sep2_new +=
                        max_len[z] - S.Unibase( S.Vert(x).U( ).back( ) ).isize( )
                        - S.Unibase( S.Vert(y).U( ).front( ) ).isize( );
                    if ( sep1_new > max_sep ) continue;
                    if ( y == p2id && sep2_new >= min_sep ) {
                        placements.push( j1, j2, swapped );
                        pp.clear( );
                        break;
                    }
                    sep1_new += S.Vert(y).Len( ), sep2_new += S.Vert(y).Len( );
                    const int insane = 100000;
                    if ( Abs(sep1_new) > insane || Abs(sep2_new) > insane ) continue;
                    if ( sep1_new <= max_sep )
                        pp.push( y, sep1_new, sep2_new );
                }
            }
        }
}

void snark::BringOutTheDead( ) {
    vec<int> to_remove;
    vec<Bool> to_remove2( VertN( ), False );
    for ( int i = 0; i < VertN( ); i++ ) {
        if ( Vert(i).Dead( ) ) {
            to_remove.push_back(i);
            to_remove2[i] = True;
        }
    }
    G_.RemoveEdgelessVertices(to_remove);
    EraseIf( seq_, to_remove2 );
    Gmutable( ).RemoveDeadEdgeObjects( );
}

void placement_on::Print( std::ostream& out, const snark& S, const vec<int>& to_left,
                          const vec<int>& to_right ) {
    out << "position " << pos << "/" << apos << " on ";
    if ( id < S.VertN( ) ) {
        out << "vert " << S.Vert(id).U( ).front( )
            << ".." << S.Vert(id).U( ).back( );
    } else {
        int e = id - S.VertN( );
        out << "edges between verts " << S.Vert( to_left[e] ).U( ).front( )
            << ".." << S.Vert( to_left[e] ).U( ).back( ) << " and "
            << S.Vert( to_right[e] ).U( ).front( ) << ".."
            << S.Vert( to_right[e] ).U( ).back( );
    }
}

Bool uniseq::Contains( const uniseq& x ) const {
    for ( int p = 0; p < N( ); p++ ) {
        if ( U( ).Contains( x.U( ), p ) && Over( ).Contains( x.Over( ), p ) )
            return True;
    }
    return False;
}

void Print( std::ostream& out, const vec< vec<uniseq> >& ul, const int K ) {
    for ( int j = 0; j < ul.isize( ); j++ ) {
        if ( j > 0 ) out << " --> ";
        const vec<uniseq>& x = ul[j];
        out << "{";
        for ( int k = 0; k < x.isize( ); k++ ) {
            if ( k > 0 ) out << ",";
            x[k].Print( out, K );
        }
        out << "}";
    }
}

vec<int> Common( const vec<uniseq>& v ) {
    vec< vec<int> > c( v.size( ) );
    for ( int i = 0; i < v.isize( ); i++ ) {
        c[i] = v[i].U( );
        UniqueSort( c[i] );
    }
    vec<int> intersection;
    Intersection( c, intersection );
    return intersection;
}

basevector uniseq::Bases( ) const {
    Assert( Alive( ) );
    basevector b = Unibase( U(0) );
    for ( int j = 1; j < N( ); j++ ) {
        b.resize( b.isize( ) - Over(j-1) );
        b = Cat( b, Unibase( U(j) ) );
    }
    return b;
}

int64_t snark::EstimatedGenomeSize( ) const {
    int64_t g = 0;
    for ( int v = 0; v < VertN( ); v++ )
        g += Vert(v).Len( );
    for ( int e = 0; e < EdgeN( ); e++ ) {
        g += Edge(e).MidLen( );
        if ( Edge(e).Closed( ) ) {
            g -= Unibase( Edge(e).Closure(0).U( ).front( ) ).isize( );
            g -= Unibase( Edge(e).Closure(0).U( ).back( ) ).isize( );
        }
    }
    return g;
}

void snark::HandleSimpleInvertedRepeats( ) {
    vec<int> multi;
    for ( int x = 0; x < VertN( ); x++ ) {
        if ( From(x).empty( ) && To(x).empty( ) ) continue;
        if ( !From(x).solo( ) || !To(x).solo( ) ) multi.push_back(x);
    }
    if ( multi.size( ) == 2
            && Vert( multi[1] ) == Vert( multi[0] ).Reverse( ToRc( ) ) ) {
        int x = multi[0], rcx = multi[1];
        if ( From(x).size( ) == 2 && To(x).size( ) == 2 ) {
            vec< vec<int> > chains(2);
            Bool fail = False;
            for ( int pass = 0; pass < 2; pass++ ) {
                int y = ( pass == 0 ? x : rcx );
                while(1) {
                    chains[pass].push_back(y);
                    if ( y == ( pass == 0 ? rcx : x ) ) break;
                    if ( chains[pass].size( ) > 1 && From(y).size( ) != 1 ) {
                        fail = True;
                        break;
                    }
                    y = From(y)[0];
                }
            }
            if ( !fail ) {
                if ( chains[1].size( ) > chains[0].size( ) )
                    swap( chains[0], chains[1] );
                for ( int j = 0; j < chains[0].isize( ) - 1; j++ ) {
                    int a = chains[0][j], b = chains[0][j+1];
                    for ( int l = 0; l < From(a).isize( ); l++ ) {
                        if ( From(a)[l] == b ) {
                            Gmutable( ).DeleteEdgeFrom( a, l );
                        }
                    }
                }
                for ( int j = 1; j < chains[0].isize( ) - 1; j++ )
                    VertMutable( chains[0][j] ).Kill( );
            }
        }
    }
    BringOutTheDead( );
}

void snark::RemoveSubsumedStuff( ) {
    for ( int x = 0; x < VertN( ); x++ ) {
        if ( From(x).isize( ) > 0 || To(x).isize( ) > 0 ) continue;
        uniseq un = Vert(x);
        for ( int pass = 1; pass <= 2; pass++ ) {
            if ( Vert(x).Dead( ) ) break;
            for ( int y = 0; y < VertN( ); y++ ) {
                if ( Vert(y).Contains(un) && Vert(y) != un ) {
                    VertMutable(x).Kill( );
                    break;
                }
            }
            un.ReverseMe( ToRc( ) );
        }
        if ( Vert(x).Dead( ) || Vert(x).N( ) != 1 ) continue;
        for ( int e = 0; e < EdgeN( ); e++ ) {
            if ( !Edge(e).Closed( ) ) continue;
            vec<int> common = Common( Edge(e).Closures( ) );
            int u = Vert(x).U(0);
            if ( BinMember( common, u ) || BinMember( common, ToRc(u) ) ) {
                VertMutable(x).Kill( );
                break;
            }
        }
    }
    BringOutTheDead( );

    // Look for edges that fold in.

    for ( int x = 0; x < VertN( ); x++ ) {
        if ( From(x).size( ) != 1 || To(x).size( ) != 0 ) continue;
        const gapster& gx = G( ).EdgeObjectByIndexFrom( x, 0 );
        uniseq p = Vert(x);
        int y = From(x)[0];
        if ( !To(y).size( ) == 2 ) continue;
        int e = -1;
        for ( int j = 0; j < To(y).isize( ); j++ )
            if ( To(y)[j] != x ) e = j;
        if ( e < 0 ) continue;
        const gapster& gy = G( ).EdgeObjectByIndexTo( y, e );
        if ( !gy.Closed( ) ) continue;
        Bool contained = True;
        for ( int j = 0; j < gy.ClosureCount( ); j++ ) {
            const uniseq& v = gy.Closure(j);
            if ( !v.Contains(p) ) contained = False;
        }
        if (contained) {
            Gmutable( ).DeleteEdgeFrom( x, 0 );
            VertMutable(x).Kill( );
        }
    }
    BringOutTheDead( );
    for ( int x = 0; x < VertN( ); x++ ) {
        if ( To(x).size( ) != 1 || From(x).size( ) != 0 ) continue;
        const gapster& gx = G( ).EdgeObjectByIndexTo( x, 0 );
        uniseq p = Vert(x);
        int y = To(x)[0];
        if ( !From(y).size( ) == 2 ) continue;
        int e = -1;
        for ( int j = 0; j < From(y).isize( ); j++ )
            if ( From(y)[j] != x ) e = j;
        if ( e < 0 ) continue;
        const gapster& gy = G( ).EdgeObjectByIndexFrom( y, e );
        if ( !gy.Closed( ) ) continue;
        Bool contained = True;
        for ( int j = 0; j < gy.ClosureCount( ); j++ ) {
            const uniseq& v = gy.Closure(j);
            if ( !v.Contains(p) ) contained = False;
        }
        if (contained) {
            Gmutable( ).DeleteEdgeTo( x, 0 );
            VertMutable(x).Kill( );
        }
    }
    BringOutTheDead( );

    // Look for vertex-edge-vertex that is subsumed.

    vec<int> to_delete;
    for ( int x = 0; x < VertN( ); x++ ) {
        if ( From(x).size( ) != 1 || To(x).size( ) != 0 ) continue;
        const gapster& gx = G( ).EdgeObjectByIndexFrom( x, 0 );
        if ( !gx.Closed( ) ) continue;
        int y = From(x)[0];
        if ( To(y).size( ) != 1 || From(y).size( ) != 0 ) continue;
        const uniseq &p1 = Vert(x), &p3 = Vert(y);
        Bool found = False;
        for ( int j = 0; j < gx.ClosureCount( ); j++ ) {
            const uniseq& p2 = gx.Closure(j);
            vec<int> u;
            for ( int l = 0; l < p1.N( ) - 1; l++ )
                u.push_back( p1.U(l) );
            for ( int l = 0; l < p2.N( ); l++ )
                u.push_back( p2.U(l) );
            for ( int l = 1; l < p3.N( ); l++ )
                u.push_back( p3.U(l) );
            for ( int z = 0; z < VertN( ); z++ ) {
                if ( z == x || z == y ) continue;
                const vec<int>& v = Vert(z).U( );
                vec<int> vr(v);
                vr.ReverseMe( );
                for ( int l = 0; l < vr.isize( ); l++ )
                    vr[l] = ToRc( vr[l] );
                if ( v.Contains(u) || vr.Contains(u) ) {
                    found = True;
                    break;
                }
            }
            if (found) break;
        }
        if (found) to_delete.push_back( x, y );
    }
    UniqueSort(to_delete);
    for ( int j = 0; j < to_delete.isize( ); j++ ) {
        int x = to_delete[j];
        if ( From(x).solo( ) ) Gmutable( ).DeleteEdgeFrom( x, 0 );
        VertMutable(x).Kill( );
    }
    BringOutTheDead( );
}

// ComputeScaffolds
class sepdevkind {
  public:
    sepdevkind() {};
    sepdevkind( int sep, int dev, Bool open ) {
        sep_  = sep;
        dev_  = dev;
        open_ = open;
    }

    Bool Open() const {
        return open_;
    }
    Bool Closed() const {
        return ! open_;
    }
    int  Sep()  const {
        return sep_;
    }
    int  Dev()  const {
        return dev_;
    }

  private:
    int sep_;
    int dev_;
    Bool open_;
};

void snark::ComputeScaffolds( vec<superb>& superbs, VecEFasta& econtigs,
                              digraphE<sepdev>& SG ) const {

    superbs.clear();
    econtigs.clear();
    SG.Clear();

    if ( G_.N() == 0 ) {
        std::cout << Date() << " : Nothing to scaffold. Graph is empty." << std::endl;
        return;
    }
    // convert to digraphE VG having vertices assosiated with sequence (gapster)
    // and edges with separation of sequence objects (sepdev).

    std::cout << Date() << ": computing scaffolds for snark with " << G_.N()
              << " vertices and " << G_.EdgeObjectCount() << " edges" << std::endl;

    vec<int> ulens( unibases_->size(), 0);
    for ( size_t ui = 0; ui < ulens.size(); ui++ )
        ulens[ui] = (*unibases_)[ui].isize();

    digraphE<sepdevkind> VG;
    vec<gapster> vsters;
    {
        std::cout << Date() << ": converting to vertex graph" << std::endl;
        vec<int> ei2newvert( G_.EdgeObjectCount(), -1 );
        int ci = -1;
        int nGapsOpen=0;
        for ( int ei = 0; ei < G_.EdgeObjectCount(); ei++ ) {
            if ( G_.Edges()[ei].Closed() ) {
                ForceAssertGt( G_.Edges()[ei].Closures().size(), 0u );
                ci++;
                ei2newvert[ei] = G_.N() + ci;
            } else {
                nGapsOpen++;
            }
        }

        int nNewVerts = ( ei2newvert.size() == 0 || Max(ei2newvert) < 0 ) ? G_.N() : Max( ei2newvert ) + 1;
        vec< vec<int> > from(nNewVerts), to(nNewVerts), to_edge_obj(nNewVerts), from_edge_obj(nNewVerts);
        vsters.resize( nNewVerts );
        for ( int i = 0; i < G_.N(); i++ )
            vsters[i] = gapster( seq_[i] );
        for ( int ei = 0; ei < G_.EdgeObjectCount(); ei++ ) {
            if ( ei2newvert[ei] >= 0 ) {
                int nvi = ei2newvert[ei];
                vsters[nvi] = G_.Edges()[ei];
                ForceAssertGt( vsters[nvi].Closures().size(), 0u);
            }
        }

        vec<Bool> edgesUsed( G_.Edges().size(), False );
        vec<sepdevkind> newEdges;
        for ( int vi = 0; vi < G_.N(); vi++ ) {
            for ( int j = 0; j < G_.From(vi).isize(); j++ ) {
                int ei = G_.EdgeObjectIndexByIndexFrom(vi,j);
                edgesUsed[ei]=True;
                int nvi, sep, dev;
                Bool open;
                if ( ei2newvert[ei] >= 0 ) { // edge has a sequence
                    nvi =  ei2newvert[ei];
                    dev = 0;
                    vec<uniseq> closures = vsters.at(nvi).Closures();
                    ForceAssertGt( closures.size(), 0u );
                    uniseq useq = closures.at(0);
                    vec<int> us = useq.U();
                    int uid = us.front();
                    sep = -ulens[uid];
                    open = False;
                } else { // edge is a gap
                    nvi = G_.From(vi)[j];
                    sep = G_.Edges()[ei].Sep();
                    dev = G_.Edges()[ei].Dev();
                    open = True;
                }
                newEdges.push_back( sepdevkind( sep, dev, open ) );
                from[vi].push_back(nvi);
                to[nvi].push_back(vi);
                from_edge_obj[vi].push_back( newEdges.isize() -1 );
                to_edge_obj[nvi].push_back( newEdges.isize() -1 );
            }
        }

        for ( int vi = 0; vi < G_.N(); vi++ ) {
            for ( int j = 0; j < G_.To(vi).isize(); j++ ) {
                int ei = G_.EdgeObjectIndexByIndexTo(vi,j);
                edgesUsed[ei] = True;
                int pvi, sep, dev;
                Bool open;
                if ( ei2newvert[ei] >= 0 ) { // edge has a sequence
                    pvi =  ei2newvert[ei];
                    dev = 0;
                    vec<uniseq> closures = vsters.at(pvi).Closures();
                    ForceAssertGt( closures.size(), 0u );
                    uniseq useq = closures.at(0);
                    vec<int> us = useq.U();
                    int uid = us.back();
                    sep = -ulens[uid];
                    open = False;
                    newEdges.push_back( sepdevkind( sep, dev, open ) );
                    from[pvi].push_back(vi);
                    to[vi].push_back(pvi);
                    from_edge_obj[pvi].push_back( newEdges.isize() -1 );
                    to_edge_obj[vi].push_back( newEdges.isize() -1 );
                } else { // edge is a gap and should already exist in the graph. We need to find its index
                    pvi = G_.To(vi)[j];
                    int eix = -1;
                    for ( int k = 0; k < from.at(pvi).isize(); k++ ) {
                        if ( from[pvi][k] == vi ) {
                            eix = from_edge_obj.at(pvi).at(k);
                            break;
                        }
                    }
                    ForceAssertGe( eix, 0 );
                }
            }
        }
        ForceAssertEq( Sum(edgesUsed), edgesUsed.isize() );
        for ( size_t vi = 0; vi < from.size(); vi++ ) {
            SortSync( from[vi], from_edge_obj[vi] );
            SortSync( to[vi], to_edge_obj[vi] );
        }
        VG.Initialize( from, to, newEdges, to_edge_obj, from_edge_obj );
    }

    vec< std::pair<int,int> > vcuts( vsters.size() );
    for ( size_t vi = 0; vi < vsters.size(); vi++ ) {
        int luid = vsters[ vi ].Closures()[0].U().front();
        int ruid = vsters[ vi ].Closures()[0].U().back();
        vcuts[vi] = std::pair<int,int>( ulens[luid], ulens[ruid] );
    }


    // order vertices so that we preferentially start with each scaffolding process
    // with a vertex having well defined overlaps with neighbors. This should help keeping
    // open gaps in scaffolds

    std::cout << Date() << ": preparing for scaffolding, ordering vertices" << std::endl;
    vec<int> closedNeighborsCt( VG.N(), 0 );
    for ( int v = 0; v < VG.N(); v++ ) {
        if ( VG.From(v).size() == 1 && VG.Edges()[ VG.EdgeObjectIndexByIndexFrom(v,0) ].Sep() < 0 )
            closedNeighborsCt[v]++;
        if ( VG.To(v).size() == 1 && VG.Edges()[ VG.EdgeObjectIndexByIndexTo(v,0) ].Sep() < 0 )
            closedNeighborsCt[v]++;
    }
    vec<int> orderedVerts( VG.N(), vec<int>::IDENTITY );
    ReverseSortSync( closedNeighborsCt, orderedVerts );


    // start scaffolding by finding unambiguous lines in VG graph

    std::cout << Date() << ": building vertex scaffolds" << std::endl;
    vec<Bool> vertsUsed( VG.N(), False );
    int nused = 0;
    int nall = VG.N();
    vec< vec<int> > vscaffolds;

    for ( size_t oi = 0; oi < orderedVerts.size(); oi++ ) {
        int vs = orderedVerts[oi];
        if ( vertsUsed[vs] ) continue;
        vec<int> vline;
        vline.push_back( vs );
        vertsUsed[vs] = True;
        int vc = vs;
        while(1) {
            if ( VG.From(vc).size() != 1 ) break;

            // check if next vertex should be added
            int vn = VG.From(vc)[0];
            if ( vertsUsed[vn] ) break;

            if ( VG.To(vn).size() == 1 ) {
                vline.push_back(vn);
                vertsUsed[vn] = True;
                vc = vn;  // update vc
            } else break;
        }

        // now go in the opposite direction
        vc = vs;
        while(1) {
            if ( VG.To(vc).size() != 1 ) break;

            // check if next vertex should be added
            int vp = VG.To(vc)[0];
            if ( vertsUsed[vp] ) break;

            if ( VG.From(vp).size() == 1 ) {
                vline.push_front(vp);
                vertsUsed[vp] = True;
                vc = vp;  // update vc
            } else break;
        }

        vscaffolds.push_back( vline );
    }
    ForceAssertEq( VG.N(), Sum( vertsUsed ) ); // make sure all vertices were used


    // compute scaffold graph

    std::cout << Date() << ": computing scaffold graph" << std::endl;
    int Ns = vscaffolds.isize();
    vec< vec<int> > sfrom(Ns), sto(Ns), sto_edge_obj(Ns), sfrom_edge_obj(Ns);
    vec<sepdev> sedges;
    for ( int s1i = 0; s1i < Ns; s1i++ ) {
        int ve = vscaffolds[s1i].back();
        for ( int s2i = 0; s2i < Ns; s2i++ ) {
            int vb  = vscaffolds[s2i].front();
            if ( VG.HasEdge(ve,vb ) ) {
                int ei  = VG.EdgesBetween(ve,vb)[0];
                int sep = VG.Edges()[ei].Sep();
                int dev = VG.Edges()[ei].Dev();

                sepdev sedge( sep, dev );
                sedges.push_back( sedge );
                sfrom[s1i].push_back( s2i );
                sfrom_edge_obj[s1i].push_back( sedges.size() -1 );
                sto[s2i].push_back( s1i );
                sto_edge_obj[s2i].push_back( sedges.size() -1 );
            }
        }
    }
    for ( size_t vi = 0; vi < sfrom.size(); vi++ ) {
        SortSync( sfrom[vi], sfrom_edge_obj[vi] );
        SortSync( sto[vi], sto_edge_obj[vi] );
    }
    SG.Initialize( sfrom, sto, sedges, sto_edge_obj, sfrom_edge_obj );


    // compute vertex contigs

    vec<superb> supervs( vscaffolds.size() );
    vec< vec<int> > vcontigs;
    for ( int si = 0; si < vscaffolds.isize(); si++ ) {
        const vec<int>& vline = vscaffolds[si];
        superb sv;
        vec<int> vcontig;
        vcontig.push_back( vline[0] );
        for ( int vi = 1; vi < vline.isize(); vi++ ) {
            int vp = vline[vi-1];
            int vc = vline[vi];
            int ei = VG.EdgesBetween(vp,vc)[0];
            int sep = VG.Edges()[ei].Sep();
            int dev = VG.Edges()[ei].Dev();
            int open = VG.Edges()[ei].Open();
            if ( !open ) {
                vcontig.push_back( vline[vi] );
            } else {
                vcontigs.push_back(vcontig);
                sv.AppendTig( vcontigs.size() -1, vcontigs.back().isize(), sep, dev );
                vcontig.clear();
                vcontig.push_back( vline[vi] );
            }
        }
        if ( vcontig.size() > 0 ) {
            vcontigs.push_back(vcontig);
            sv.AppendTig( vcontigs.size() -1, vcontigs.back().isize() );
            vcontig.clear();
        }
        supervs[si] = sv;
    }

    vec<int> vtig2scaff( vcontigs.size(), -1 );
    for ( int svi = 0; svi < supervs.isize(); svi++ ) {
        for ( int ti = 0; ti < supervs[svi].Ntigs(); ti++ )
            vtig2scaff[ supervs[svi].Tig(ti) ] = svi;
    }
    ForceAssertGe( Min(vtig2scaff), 0 );


    // convert vcontigs to sequence

    econtigs.clear();
    econtigs.resize( vcontigs.size() );
    for ( size_t eci = 0; eci < vcontigs.size(); eci++ ) {
        vec<int>& verts = vcontigs[eci];
        int firstcut = verts.front() < G_.N() ? 0 : vcuts[ verts.front() ].first;
        int lastcut = verts.back() < G_.N() ? 0 : vcuts[ verts.back() ].second;
        // compute efastas

        efasta econtig;

        for ( size_t vi = 0; vi < verts.size(); vi++ ) {
            vec<uniseq> cs = vsters[ verts[vi] ].Closures();
            int lcut = vi == 0u ? firstcut : vcuts[ verts[vi] ].first;
            vec<BaseVec> seqs;
            for ( size_t csi = 0; csi < cs.size(); csi++ ) {
                int inisize =  cs[csi].Bases().isize();
                ForceAssertGe( inisize, lcut );
                BaseVec seq( cs[csi].Bases(), lcut, inisize - lcut );
                seqs.push_back( seq );
            }
            efasta eseq( seqs );
            econtig += eseq;
        }
        ValidateEfastaRecord( econtig );
        // check if there is enough room for the cut on the right

        vec<int> gaps;
        if ( lastcut > 0 ) {

            if ( econtig.MinLength() < lastcut ) {
                if ( econtig.AmbEventCount() == 0 ) {
                    gaps.push_back( econtig.isize() - lastcut );
                    econtig.clear();
                } else {
                    vec< vec<basevector> > G;
                    vec< vec<int> > gpaths;
                    if ( ! econtig.ExpandTo( G, gpaths ) )
                        FatalErr("Efasta expansion exploded");

                    vec<BaseVec> bpaths;
                    gaps.clear();
                    if ( ! econtig.ExpandTo( bpaths ) )
                        FatalErr("Efasta expansion exploded");
                    PRINT2( gpaths.size(), bpaths.size() );
                    ForceAssertEq( bpaths.size(), gpaths.size() );

                    vec<int> glens( gpaths.size() );
                    for ( size_t gi = 0; gi < gpaths.size(); gi++ )
                        for ( size_t j = 0; j < gpaths[gi].size(); j++ )
                            glens[gi] += G[gi][j].isize();

                    vec< vec<Bool> > to_remove( G.size() );
                    for ( size_t gi = 0; gi < G.size(); gi++ ) {
                        to_remove[gi].resize( G[gi].size() );
                        for ( size_t j = 0; j < G[gi].size(); j++ )
                            to_remove[gi][j] = True;
                    }
                    for ( size_t gi = 0; gi < gpaths.size(); gi++ ) {
                        if ( glens[gi] < lastcut ) {
                            gaps.push_back( glens[gi] - lastcut );
                        } else {
                            for ( size_t j = 0; j < gpaths[gi].size(); j++ )
                                to_remove[gi][j] = False;
                        }
                    }

                    for ( size_t gi = 0; gi < G.size(); gi++ )
                        EraseIf( G[gi], to_remove[gi] );

                    efasta econtig_new;
                    econtig_new.MakeFromGraph( G );
                    if ( econtig_new.MinLength() >= lastcut ) {
                        econtig = econtig_new;
                    } else {
                        std::cout << "Unable to remove short paths from efasta efficiently. ";
                        std::cout << "Going to allow long 'unzipped' alternatives." << std::endl;

                        vec<BaseVec> bpaths;
                        gaps.clear();
                        if ( ! econtig.ExpandTo( bpaths ) )
                            FatalErr("Efasta expansion exploded");
                        ForceAssertEq( bpaths.size(), gpaths.size() );
                        vec<BaseVec> seqs;
                        for ( size_t bi = 0; bi < bpaths.size(); bi++ ) {
                            if ( bpaths[bi].isize() >= lastcut ) {
                                BaseVec seq( bpaths[bi], 0, bpaths[bi].isize() - lastcut );
                                seqs.push_back(seq);
                            } else {
                                gaps.push_back( bpaths[bi].isize() - lastcut );
                            }
                        }

                        if ( seqs.size() > 0 ) {
                            efasta eseq( seqs );
                            econtig = eseq;
                        } else {
                            econtig.clear();
                        }
                    }
                }
            }

            if ( econtig.size() > 0 ) {
                size_t braRightLoc = econtig.rfind("}");
                if ( braRightLoc == String::npos ||
                        (ssize_t)econtig.size() - (ssize_t)braRightLoc - (ssize_t)1 >= (ssize_t)lastcut ) {
                    econtig.erase( econtig.size() - lastcut, lastcut );
                } else {
                    Bool Found = False;
                    for ( int ie = econtig.isize() -1; ie >= 0; ie-- ) {
                        if ( ! Found && ( econtig[ie] == '{' || ie == 0 || ! econtig.IndexInAmbiguity(ie) ) ) {
                            String secontigR( econtig, ie, econtig.isize() - ie );
                            efasta econtigR( secontigR );
                            if ( econtigR.MinLength() >= lastcut ) {
                                Found = True;
                                vec<BaseVec> bpaths;
                                if ( ! econtigR.ExpandTo( bpaths ) )
                                    FatalErr("Efasta expansion exploded");;
                                vec<BaseVec> seqs;
                                for ( size_t bi = 0; bi < bpaths.size(); bi++ ) {
                                    ForceAssertGe( bpaths[bi].isize(), lastcut );
                                    BaseVec seq( bpaths[bi], 0, bpaths[bi].isize() - lastcut );
                                    seqs.push_back(seq);
                                }
                                efasta eseq( seqs );
                                String secontigL( econtig, 0, ie );
                                efasta econtigL( secontigL );
                                efasta echeck = econtigL + econtigR;
                                ForceAssertEq( econtig, echeck );
                                econtig = efasta(econtigL) + eseq;
                            }
                        }
                    }
                    ForceAssert( Found );
                }
            }
        }

        if ( econtig.size() > 0 )
            ValidateEfastaRecord(econtig);
        econtigs[eci] = econtig;

        // update scaffold graph

        if ( gaps.isize() > 0 ) {
            int svi = vtig2scaff[eci];
            ForceAssertEq( supervs[svi].Ntigs(), 1 );
            int mean = Mean(gaps);
            int dev = StdDev(gaps,mean);
            for ( int ii = 0; ii < SG.To(svi).isize(); ii++ ) {
                int vp = SG.To(svi)[ii];
                for ( int jj = 0; jj < SG.From(svi).isize(); jj++ ) {
                    int vn = SG.From(svi)[jj];
                    SG.AddEdge( vp, vn, sepdev(mean,dev) );
                }
            }
        }

        int svi = vtig2scaff[eci];
        int tigFirst = supervs[svi].Tig(0);
        if ( (int)eci == tigFirst ) {
            for ( int ii = 0; ii < SG.ToEdgeObj(svi).isize(); ii++ ) {
                int ei = SG.ToEdgeObj(svi)[ii];
                SG.EdgesMutable()[ei].AddToSep(firstcut);
            }
        }
        int tigLast = supervs[svi].Tig( supervs[svi].Ntigs() -1 );
        if ( (int)eci == tigLast ) {
            for ( int ii = 0; ii < SG.FromEdgeObj(svi).isize(); ii++ ) {
                int ei = SG.FromEdgeObj(svi)[ii];
                SG.EdgesMutable()[ei].AddToSep(lastcut);
            }
        }
    }


    superbs.clear();
    superbs = supervs;
    longlong eLenSum = 0;
    for ( size_t si = 0; si < superbs.size(); si++ ) {
        for ( int ti = 0; ti < superbs[si].Ntigs(); ti++ ) {
            int tig = superbs[si].Tig(ti);
            superbs[si].SetLen( ti, econtigs[tig].Length1() );
            eLenSum += econtigs[tig].size();
        }
    }

    PRINT2( superbs.size(), econtigs.size() );
    int fullLength = 0, reducedLength = 0;
    for ( size_t si = 0; si < superbs.size(); si++ ) {
        PRINT2( si, superbs[si].Ntigs() );
        superbs[si].Print(std::cout,ToString(si));
        fullLength += superbs[si].FullLength();
        reducedLength += superbs[si].ReducedLength();
    }
    PRINT2( fullLength, reducedLength );
    PRINT( eLenSum );
    std::cout << Date() << " : Done with scaffolding" << std::endl;
    return;
}

void uniseq::TrimLeft( const int n ) {
    AssertGe( u_.isize(), n );
    u_ = SubOf( u_, n, u_.isize( ) - n );
    if ( overlap_.isize() >= n )
        overlap_ = SubOf( overlap_, n, overlap_.isize( ) - n );
    else overlap_.clear();
}

void uniseq::TrimRight( const int n ) {
    AssertGe( u_.isize(), n );
    u_ = SubOf( u_, 0, u_.isize() - n );
    if ( overlap_.isize() >= n )
        overlap_ = SubOf( overlap_, 0, overlap_.isize() - n );
    else overlap_.clear();
}

void uniseq::TrimEnds( const int n, const int m ) {
    AssertGe( u_.isize(), n+m );
    u_ = SubOf( u_, n, u_.isize() - m - n );
    if ( overlap_.isize() >= n+m )
        overlap_ = SubOf( overlap_, n, overlap_.isize() - m - n );
    else overlap_.clear();
}

void uniseq::writeBinary( BinaryWriter& writer ) const {
    writer.write(u_);
    writer.write(overlap_);
}

void uniseq::readBinary( BinaryReader& reader ) {
    reader.read(&u_);
    reader.read(&overlap_);
}

void gapster::writeBinary( BinaryWriter& writer ) const {
    writer.write(open_);
    writer.write(sep_);
    writer.write(dev_);
    writer.write(closures_);
}

void gapster::readBinary( BinaryReader& reader ) {
    reader.read(&open_);
    reader.read(&sep_);
    reader.read(&dev_);
    reader.read(&closures_);
}

void snark::writeBinary( BinaryWriter& writer ) const {
    writer.write(G_);
    writer.write(seq_);
}

void snark::readBinary( BinaryReader& reader ) {
    reader.read(&G_);
    reader.read(&seq_);
}

void snark::AssignVertexColors( vec<int>& color ) const {
    color.resize_and_set( VertN( ), 0 );
    vec<int> to_rc( VertN( ), -1 );
    for ( int v1 = 0; v1 < VertN( ); v1++ )
        for ( int v2 = 0; v2 < VertN( ); v2++ ) {
            if ( v1 == v2 ) continue;
            uniseq x1 = Vert(v1), x2 = Vert(v2).Reverse( ToRc( ) );
            if ( x1 == x2 ) to_rc[v1] = v2;
        }
    while(1) {
        Bool progress = False;
        for ( int v1 = 0; v1 < VertN( ); v1++ ) {
            if ( color[v1] == 0 && to_rc[v1] >= 0 && color[ to_rc[v1] ] == 0 ) {
                color[v1] = 1, color[ to_rc[v1] ] = 2;
                progress = True;
                break;
            }
        }
        while(1) {
            Bool progress2 = False;
            for ( int v = 0; v < VertN( ); v++ ) {
                if ( color[v] != 1 ) continue;
                vec<int> n;
                for ( int j = 0; j < From(v).isize( ); j++ ) {
                    int w = From(v)[j];
                    if ( color[w] == 0 && to_rc[w] != 0 )
                        n.push_back(w);
                }
                for ( int j = 0; j < To(v).isize( ); j++ ) {
                    int w = To(v)[j];
                    if ( color[w] == 0 && to_rc[w] != 0 )
                        n.push_back(w);
                }
                UniqueSort(n);
                for ( int j = 0; j < n.isize( ); j++ ) {
                    int w = n[j];
                    int z = to_rc[w];
                    if ( color[w] > 0 || z < 0 || color[z] > 0 ) continue;
                    Bool ok = False;
                    for ( int j = 0; j < From(z).isize( ); j++ ) {
                        int t = From(z)[j];
                        if ( color[t] == 2 ) ok = True;
                    }
                    for ( int j = 0; j < To(z).isize( ); j++ ) {
                        int t = To(z)[j];
                        if ( color[t] == 2 ) ok = True;
                    }
                    if ( !ok ) continue;
                    color[w] = 1, color[z] = 2;
                    progress = True, progress2 = True;
                }
            }
            if ( !progress2 ) break;
        }
        if ( !progress ) break;
    }
}

void snark::DeleteComplements( ) {
    vec<int> color;
    AssignVertexColors(color);
    vec<String> vertex_colors( VertN( ) );

    // Test to see if blue is connected, red is connected, and blue is
    // connected to red.

    vec<int> blues, reds;
    for ( int v = 0; v < VertN( ); v++ ) {
        if ( color[v] == 1 ) blues.push_back(v);
        if ( color[v] == 2 ) reds.push_back(v);
    }
    equiv_rel e_blues( blues.size( ) ), e_reds( reds.size( ) );
    Bool tied = False;
    for ( int i = 0; i < blues.isize( ); i++ ) {
        for ( int j = 0; j < From( blues[i] ).isize( ); j++ ) {
            int p = BinPosition( blues, From( blues[i] )[j] );
            if ( p >= 0 ) e_blues.Join( i, p );
            int q = BinPosition( reds, From( blues[i] )[j] );
            if ( q >= 0 ) tied = True;
        }
    }
    for ( int i = 0; i < reds.isize( ); i++ ) {
        for ( int j = 0; j < From( reds[i] ).isize( ); j++ ) {
            int p = BinPosition( reds, From( reds[i] )[j] );
            if ( p >= 0 ) e_reds.Join( i, p );
            int q = BinPosition( blues, From( reds[i] )[j] );
            if ( q >= 0 ) tied = True;
        }
    }
    if ( tied && e_blues.OrbitCount( ) == 1 && e_reds.OrbitCount( ) == 1 ) {
        // Delete red.

        for ( int r = 0; r < reds.isize( ); r++ ) {
            int v = reds[r];
            Gmutable( ).DeleteEdgesAtVertex(v);
            VertMutable(v).Kill( );
        }
        BringOutTheDead( );
    }
}

#include "graph/DigraphTemplate.h"
template digraphE<gapster>::digraphE(const vec<vec<int> >&, const vec<vec<int> >&, const vec<gapster>&, const vec<vec<int> >&, const vec<vec<int> >&, const Bool);
template int digraphE<gapster>::AddEdge(const int, const int, const gapster&);
template void digraphE<gapster>::DeleteEdgeFrom(int, int);
template void digraphE<gapster>::DeleteEdges(const vec<int>&);
template void digraphE<gapster>::DeleteEdgesAtVertex(int);
template void digraphE<gapster>::DeleteEdgeTo(int, int);
template const gapster& digraphE<gapster>::EdgeObject(int) const;
template gapster const& digraphE<gapster>::EdgeObjectByIndexFrom(int, int) const;
template gapster& digraphE<gapster>::EdgeObjectByIndexFromMutable(int, int);
template gapster const& digraphE<gapster>::EdgeObjectByIndexTo(int, int) const;
template gapster& digraphE<gapster>::EdgeObjectByIndexToMutable(int, int);
template int digraphE<gapster>::EdgeObjectCount() const;
template int digraphE<gapster>::EdgeObjectIndexByIndexFrom(int, int) const;
template int digraphE<gapster>::EdgeObjectIndexByIndexTo(int, int) const;
template gapster& digraphE<gapster>::EdgeObjectMutable(int);
template vec<gapster> const& digraphE<gapster>::Edges() const;
template vec<int> const& digraphE<gapster>::FromEdgeObj(int) const;
template vec<int>& digraphE<gapster>::FromEdgeObjMutable(int);
template void digraphE<gapster>::Initialize(vec<vec<int> > const&, vec<vec<int> > const&, vec<gapster> const&, vec<vec<int> > const&, vec<vec<int> > const&, const Bool);
template int digraphE<gapster>::InputFromOutputTo(int, int) const;
template int digraphE<gapster>::InputToOutputFrom(int, int) const;
template vec<int> digraphE<gapster>::RemoveDeadEdgeObjects();
template void digraphE<gapster>::RemoveEdgelessVertices(const vec<int>&);
template vec<int>& digraphE<gapster>::ToEdgeObjMutable(int);
template void digraphE<gapster>::ToLeft(vec<int>&) const;
template void digraphE<gapster>::ToRight(vec<int>&) const;
template void digraphE<gapster>::Used(vec<unsigned char>&) const;
template void digraphE<gapster>::readBinary(BinaryReader&);
template void digraphE<gapster>::writeBinary(BinaryWriter&) const;

template int digraphE<sepdevkind>::EdgeObjectIndexByIndexTo(int, int) const;
template vec<int> digraphE<sepdevkind>::EdgesBetween(const int, const int) const;
template void digraphE<sepdevkind>::Initialize(const vec<vec<int> >&, const vec<vec<int> >&, const vec<sepdevkind>&, const vec<vec<int> >&, const vec<vec<int> >&, const Bool);
