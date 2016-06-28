///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "paths/long/DisplayTools.h"
#include "paths/long/large/GapToyTools.h"
#include "paths/long/large/tools/NhoodInfoStuff.h"

void TestDot( ) {
    String dot = LineOfOutput( "dot -V", true, true );
    Bool ok = True;
    if ( !dot.Contains( "version " )
            || !dot.After( "version " ).Contains( " " ) ) {
        ok = False;
    }
    if (ok) {
        String ver = dot.Between( "version ", " " );
        if ( !ver.Contains( "." ) || !ver.After( "." ).Contains( "." ) )
            ok = False;
        if (ok) {
            String v1 = ver.Before( "." );
            String v2 = ver.Between( ".", "." );
            if ( !v1.IsInt( ) || !v2.IsInt( ) )
                ok = False;
            else {
                int i1 = v1.Int( ), i2 = v2.Int( );
                if ( i1 > 2 || ( i1 == 2 && i2 >= 28 ) ) ok = True;
                else {
                    ok = False;
                    std::cout << "see dot version " << ver
                              << std::endl;
                }
            }
        }
    }
    if ( !ok ) {
        std::cout << "For PDF and PNG options, you need dot version "
                  << "2.28 or later.\nIt appears that either your version of "
                  << "dot is too old\nor else you don't have it in your path "
                  << "at all." << std::endl;
        Scram(1);
    }
}

void CreateEdgeLabels( const HyperBasevectorX& hb, const vec<int>& inv,
                       const vec<vec<vec<vec<int>>>>& lines, const vec<int>& tol,
                       const vec<int>& llens, const vec<covcount>& cov, const vec<vec<covcount>>& covs,
                       vec< vec< std::pair<int,int> > >& hits, const vec<String>& subsam_names,
                       const vec<vec<int>>& count, const vec<String>& genome_names,
                       const vec<int>& used, const nhood_info_state& state, vec<String>& edge_names ) {
    edge_names.resize( hb.E( ) );
    for ( int e = 0; e < hb.E( ); e++ ) {
        if ( !BinMember( used, e ) ) continue;
        if ( hb.Bases(e) == 0 ) continue;
        int pe = e;
        if (state.REL) pe = 1 + BinPosition( used, e );
        edge_names[e] = ToString(pe);
        if (state.SHOW_INV) {
            edge_names[e] += "<" + ToString( inv[e] );
            if (state.LINENO) edge_names[e] += " L" + ToString( tol[inv[e]] );
            edge_names[e] += ">";
        }
        if (state.LINENO) edge_names[e] += " L" + ToString( tol[e] );
        if ( state.SHOW_CN ) {
            if ( covs.nonempty( ) ) {
                Bool defined = False;
                int ns = covs.size( );
                for ( int ss = 0; ss < ns; ss++ )
                    if ( covs[ss][e].Def( ) ) defined = True;
                if (defined) {
                    edge_names[e] += " [";
                    for ( int ss = 0; ss < ns; ss++ ) {
                        if ( ss > 0 ) edge_names[e] += ";";
                        if ( covs[ss][e].Def( ) ) {
                            edge_names[e] +=
                                ToString( covs[ss][e].Cov( ), 2 ) + "x";
                        } else edge_names[e] += "?x";
                    }
                    edge_names[e] += "]";
                }
            } else if ( cov.nonempty( ) ) { // for backward compatibility
                if ( cov[e].Def( ) ) {
                    edge_names[e] +=
                        " [" + ToString( cov[e].Cov( ), 2 ) + "x]";
                }
            }
        }
        if ( state.SHOW_ALIGN && hits.nonempty( ) )
            edge_names[e] += PrintHits( e, hits, hb, inv, genome_names );
        if ( state.COUNT && count.nonempty( ) ) {
            edge_names[e] += " ";
            for ( int j = 0; j < subsam_names.isize( ); j++ ) {
                if ( j > 0 ) edge_names[e] += "/";
                edge_names[e] += subsam_names[j];
            }
            edge_names[e] += "=";
            for ( int j = 0; j < subsam_names.isize( ); j++ ) {
                if ( j > 0 ) edge_names[e] += "/";
                edge_names[e] += ToString( count[j][e] );
            }
        }
        if (state.BNG) {
            String cut = "GCTCTTC", rcut;
            StringReverseComplement( cut, rcut );
            vec<int> sites;
            String E = hb.EdgeObject(e).ToString( );
            for ( int j = 0; j < E.isize( ); j++ ) {
                if ( E.Contains( cut, j ) ) sites.push_back(j);
                if ( E.Contains( rcut, j ) ) sites.push_back(j);
            }
            if ( sites.nonempty( ) ) {
                std::ostringstream out;
                out << " BNG={" << printSeq(sites) << "}";
                edge_names[e] += out.str( );
            }
        }
    }
}

void ColorEdges( const HyperBasevectorX& hb, const vec<vec<int>>& count,
                 const vec<int>& seeds, const vec<Bool>& invisible,
                 const nhood_info_state& state, vec<String>& edge_color, vec<int>& pen_widths ) {
    edge_color.resize( hb.E( ) );
    for ( int e = 0; e < hb.E( ); e++ ) {
        if ( !invisible[e] ) {
            if ( hb.Bases(e) == 0 ) edge_color[e] = "brown";
            if ( state.PURPLE.size( ) > 0 ) {
                vec<Bool> x( count.size( ), False );
                for ( int i = 0; i < count.isize( ); i++ )
                    if ( state.PURPLE[i] == '1' ) x[i] = True;
                Bool OK = True;
                for ( int i = 0; i < count.isize( ); i++ ) {
                    if ( x[i] && count[i][e] == 0 ) OK = False;
                    if ( !x[i] && count[i][e] > 0 ) OK = False;
                }
                if (OK) {
                    edge_color[e] = "\"#800080\"";
                    pen_widths[e] = 10;
                }
            }
        }
    }
    for ( int i = 0; i < seeds.isize( ); i++ )
        if (state.GREEN) edge_color[ seeds[i] ] = "\"#00FF00\"";
}

Bool DefineSeeds( const HyperBasevectorX& hb, const vec<int>& inv,
                  const vec< triple<kmer<20>,int,int> >& kmers_plus,
                  const vec<vec<vec<vec<int>>>>& lines,
                  const vec<int>& tol, const vec<String>& genome_names,
                  const vec< std::pair<int,ho_interval> >& ambint, Bool& ambflag,
                  const vec< vec< std::pair<int,int> > >& hits, const nhood_info_state& state,
                  const int RANDOM_SEED, const String& SEEDS_MINUS, vec<int>& seeds,
                  const int max_seeds, std::ostream& tout ) {
    String SEEDSX(state.SEEDS);
    if ( state.SEEDS.Contains( "," ) && state.SEEDS[0] != '{' )
        SEEDSX = "{" + state.SEEDS + "}";
    if ( !ParseSeeds( hb, inv, kmers_plus, lines, genome_names, ambint, ambflag,
                      hits, SEEDSX, RANDOM_SEED, SEEDS_MINUS, seeds, max_seeds, tout ) ) {
        return False;
    }
    if (state.EXT) {
        vec<int> seeds2;
        for ( int s = 0; s < seeds.isize( ); s++ ) {
            const vec<vec<vec<int>>>& L = lines[ tol[ seeds[s] ] ];
            for ( int i = 0; i < L.isize( ); i++ )
                for ( int j = 0; j < L[i].isize( ); j++ )
                    for ( int k = 0; k < L[i][j].isize( ); k++ )
                        seeds2.push_back( L[i][j][k] );
        }
        seeds = seeds2;
        UniqueSort(seeds);
    }
    return True;
}

void MakeDot( const HyperBasevectorX& hb, const vec<int>& inv,
              const vec<vec<vec<vec<int>>>>& lines, const vec<int>& tol,
              const vec<int>& llens, const vec<covcount>& cov, const vec<vec<covcount>>& covs,
              vec< vec< std::pair<int,int> > >& hits, const vec<String>& genome_names,
              const vec<int>& used, const nhood_info_state& state, const vec<int>& seeds,
              const vec<Bool>& invisible, const vec<String>& subsam_names,
              const vec<vec<int>>& count, const Bool LABEL_CONTIGS, const String& dot_file ) {
    vec<String> edge_color, edge_names;
    vec<int> pen_widths( hb.E( ), 0 );
    ColorEdges( hb, count, seeds, invisible, state, edge_color, pen_widths );
    CreateEdgeLabels( hb, inv, lines, tol, llens, cov, covs, hits,
                      subsam_names, count, genome_names, used, state, edge_names );
    vec<Bool> hide( hb.E( ), False );
    const Bool DOT_LABEL_VERTICES = False;
    vec<double> lengths( hb.E( ) );
    for ( int i = 0; i < hb.E( ); i++ )
        lengths[i] = hb.Kmers(i);
    Ofstream( dout, dot_file );
    hb.PrettyDOT( dout, lengths, HyperBasevectorX::edge_label_info(
                      HyperBasevectorX::edge_label_info::DIRECT, &edge_names ),
                  LABEL_CONTIGS, DOT_LABEL_VERTICES, NULL, NULL, NULL, &hide,
                  &invisible, &edge_color, &pen_widths, ( state.NEATO ? "neato" : "" ),
                  hb.K( ) + 50, state.FONTSIZE, state.SCALE, state.ASPECT );
}

void MakeFasta( const HyperBasevectorX& hb, const vec<int>& used,
                const String& fasta_file ) {
    Ofstream( fout, fasta_file );
    for ( int i = 0; i < used.isize( ); i++ ) {
        int e = used[i];
        hb.EdgeObject(e).Print( fout,
                                ToString(e) + ": " + ToString( hb.ToLeft(e) ) + " --> "
                                + ToString( hb.ToRight(e) ) );
    }
}
