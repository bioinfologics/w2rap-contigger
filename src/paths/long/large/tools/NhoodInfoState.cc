///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "PrintAlignment.h"
#include "TokenizeString.h"
#include "pairwise_aligners/SmithWatAffine.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "paths/long/large/tools/NhoodInfoState.h"

Bool nhood_info_state::SetState( const String & line, std::ostream& out,
                                 const Bool SERVER, const HyperBasevectorX& hb,
                                 const NhoodInfoEngine& engine ) {
    String breaker = "============"
                     "======================================================";
    if ( line == "H" ) {
        out << "Available commands:\n";
        out << "- H: help (see also the file doc/NhoodInfo)" << std::endl;
        out << "- H CORE: core help" << std::endl;
        out << "- H SEEDS: detailed help for SEEDS" << std::endl;
        out << "- RESET: reset to default arguments" << std::endl;
        out << "- CLEAR: clear screen" << std::endl;
        out << "- BASES=n: show the base sequence for edge n" << std::endl;
        if ( !SERVER )
            out << "- READ=id: show where read is placed on the assembly" << std::endl;
        out << "- ALIGN a,b: align edges a and b, which are assumed to comprise "
            << "a\n"
            << "  bubble." << std::endl;
        out << "For both BASES and ALIGN, if relative numbering is in use\n(via "
            << "REL option), then the edge numbers are relative." << std::endl;
        if ( !SERVER ) {
            out << breaker << std::endl;
            out << "Argument values persist until changed." << std::endl;
            out << "Ctrl-C to exit." << std::endl;
        }
        return False;
    }
    if ( line == "H CORE" ) {
        out << "Core command type: a blank-separated list of arg=value "
            << "constructs,\n"
            << "where arg is:" << std::endl;
        out << "- SEEDS or S: seeds" << std::endl;
        if ( !SERVER ) out << "- +s: add s to SEEDS" << std::endl;
        out << "- DEPTH or D: depth (nonnegative integer)" << std::endl;
        out << "- EXT OR E: default False, or True "
            << "to extend to end of lines" << std::endl;
        out << "- COUNT: default False, or True to show read ";
        out      << "count for each edge" << std::endl;
        out << "- GREEN or G: default False, or True to make seeds green" << std::endl;
        out << "- PURPLE: for multi-sample assemblies, if you supply a string "
            << "of zeros and ones,\n"
            << "  then edges exhibiting the corresponding pattern of presence "
            << "or absence will be\n"
            << "  shown as purple and bold; for example in a two-sample "
            << "assembly, PURPLE=10 will\n"
            << "  flag edges present in only the first sample" << std::endl;
        out << "- BNG: default False, or True to show BioNano Genomics "
            << "restriction sites" << std::endl;
        out << "- NEATO or N: default False, or True to use neato layout" << std::endl;
        out << "- FONTSIZE or FS: font size (positive number)" << std::endl;
        out << "- SCALE: graph feature scale (positive number)" << std::endl;
        out << "- REL: default False, or True to "
            << "use relative numbering for edges" << std::endl;
        out << "- LINENO: show line numbers" << std::endl;
        out << "- SHOW_INV: default False, or True to show inverse edge ids"
            << " (ignores REL)" << std::endl;
        out << "- SHOW_CN: default True, to show predicted edge copy number"
            << std::endl;
        out << "- SHOW_ALIGN: default True, to show edge alignments" << std::endl;
        out << "- ALTREF: default False, to use alternate reference for\n"
            << "  alignments if available" << std::endl;
        out << "- ASPECT: default -1, if set provides aspect ratio" << std::endl;
        out << "- SVG: default False, if True generate svg file as output"
            << std::endl;
        return False;
    }

    if ( line == "H SEEDS" ) {
        out << "usage: SEEDS=... where ... is a comma-separated list of the "
            << "following types:\n";
        if ( !SERVER ) out << "- all\n";
        out << "- g:a, to get all edges mapping to the given base on grch38\n"
            << "  where g = 1,...,22,X,Y,M ";
        if ( !SERVER ) out << "or an 'extra' record as in genome.names";
        out << "\n"
            << "- g:a-b, to get all edges mapping to the given "
            << "interval on grch38\n"
            << "- a single edge id x\n"
            << "- x..y, means all edges between edges x and y\n";
        if ( !SERVER ) out << "- Ln, means line n\n";
        out << "- random:n, to get n random edges\n"
            << "- random:n:l, to get n random edges of size >= l kmers each\n"
            << "- trandom:n, to get n clock-seeded random edges\n"
            << "- trandom:n:l, to get n clock-seeded random edges of "
            << "size >= l kmers\n"
            << "- a sequence of >= 20 bases ACGT, to represent all edges that "
            << "exactly match it" << std::endl;
        if ( !SERVER )
            out << "  (note slow unless SEQ_LOOKUP=True specified)" << std::endl;
        return False;
    }

    if ( line == "CLEAR" ) {
        out << "[2J" << std::endl;
        return False;
    }

    // Process READ.

    if ( line.Contains( "READ=", 0 ) && line.After( "READ=" ).IsInt( ) ) {
        int64_t rid = line.After( "READ=" ).Int( );
        String paths = engine.Dir( ) + "/a.paths";
        if ( !IsRegularFile(paths) ) {
            out << "I can't run this command without having the file a.paths."
                << std::endl;
            return False;
        }
        int64_t N = MastervecFileObjectCount(paths);
        if ( rid < 0 || rid >= N ) {
            out << "That read id doesn't make sense." << std::endl;
            return False;
        }
        ReadPathVec p;
        p.ReadOne( paths, rid );
        if ( p[0].size( ) == 0 ) out << "unplaced" << std::endl;
        else {
            int offset = p[0].getOffset( );
            out << "placed fw starting at position " << offset
                << " on edge " << p[0][0] << " (len=" << hb.Bases(p[0][0])
                << ")" << std::endl;
            out << "full path = " << printSeq( p[0] ) << std::endl;
        }
        return False;
    }

    // Process BASES.

    if ( line.Contains( "BASES=", 0 ) && line.After( "BASES=" ).IsInt( ) ) {
        int e = line.After( "BASES=" ).Int( );
        if (REL) {
            if ( e <= 0 || e > used.isize( ) ) {
                out << "You're in REL mode, and those edge numbers don't "
                    << "make sense." << std::endl;
                return False;
            }
            e = used[e-1];
        }
        if ( e < 0 || e > hb.E( ) )
            out << "Your edge id doesn't make sense." << std::endl;
        else {
            const basevector& b = hb.EdgeObject(e);
            out << hb.Bases(e) << " bases, " << hb.Kmers(e) << " kmers" << std::endl;
            for ( int i = 0; i < b.isize( ); i++ ) {
                if ( i > 0 && i % 80 == 0 ) out << std::endl;
                out << as_base( b[i] );
            }
            out << std::endl;
        }
        return False;
    }

    // Process ALIGN.

    if ( line.Contains( "ALIGN ", 0 ) && line.Contains( "," )
            && line.Between( " ", "," ).IsInt( ) && line.After( "," ).IsInt( ) ) {
        int e1 = line.Between( " ", "," ).Int( );
        int e2 = line.After( "," ).Int( );
        if ( e1 < 0 || e2 < 0 ) {
            out << "Negative edge numbers don't make sense." << std::endl;
            return False;
        }
        if (REL) {
            if ( e1 == 0 || e2 == 0 || e1 > used.isize( ) || e2 > used.isize( ) ) {
                out << "You're in REL mode, and those edge numbers don't "
                    << "make sense." << std::endl;
                return False;
            }
            e1 = used[e1-1], e2 = used[e2-1];
        }
        if ( e1 >= hb.E( ) || e2 >= hb.E( ) ) {
            out << "Edge numbers are too large." << std::endl;
            return False;
        }
        if ( hb.ToLeft(e1) != hb.ToLeft(e2) || hb.ToRight(e1) != hb.ToRight(e2) ) {
            out << "Edges have to comprise a bubble." << std::endl;
            return False;
        }
        const basevector &E1 = hb.EdgeObject(e1), &E2 = hb.EdgeObject(e2);
        alignment al;
        SmithWatAffine( E1, E2, al );
        align a = al;
        std::ostringstream xout;
        PrintVisualAlignment( True, xout, E1, E2, a );
        String s = xout.str( );
        s.GlobalReplaceBy( "\n\n\n", "\n\n" );
        out << s;
        return False;
    }

    // Proceed.

    vec<String> x;
    Tokenize( line, ' ', x );
    for ( int i = 0; i < x.isize( ); i++ ) {
        String s = x[i];
        if ( s == "RESET" ) {
            Initialize( );
            return False;
        }

        else if ( s.Contains( "SEEDS=", 0 ) || s.Contains( "S=", 0 ) ) {
            String SEEDS_NEW;
            if ( s.Contains( "SEEDS=", 0 ) ) SEEDS_NEW = s.After( "SEEDS=" );
            else SEEDS_NEW = s.After( "S=" );
            if ( SEEDS_NEW == SEEDS ) {
                out << "Sure, I'll do that, but you're setting the seed "
                    << "value to what it just was." << std::endl;
            }
            SEEDS = SEEDS_NEW;
        }

        else if ( s.Contains( "+", 0 ) ) SEEDS += "," + s.After( "+" );
        else if ( s.Contains( "DEPTH=", 0 ) && s.After( "DEPTH=" ).IsInt( )
                  && s.After( "DEPTH=" ).Int( ) >= 0 ) {
            DEPTH = s.After( "DEPTH=" ).Int( );
        } else if ( s.Contains( "D=", 0 ) && s.After( "D=" ).IsInt( )
                    && s.After( "D=" ).Int( ) >= 0 ) {
            DEPTH = s.After( "D=" ).Int( );
        }

        else if ( s == "EXT=True" || s == "E=True" ) {
            if ( !engine.HasLines( ) ) {
                out << "Lines were not computed for this assembly directory,\n"
                    << "so you can't set EXT=True or E=True." << std::endl;
                return False;
            }
            EXT = True;
        } else if ( s == "EXT=False" || s == "E=False" ) EXT = False;

        else if ( s == "COUNT=True" ) {
            if ( !engine.HasCounts( ) ) {
                out << "COUNT=True only works if there is a file a.counts"
                    << std::endl;
                return False;
            }
            COUNT = True;
        } else if ( s == "COUNT=False" ) COUNT = False;
        else if ( s == "GREEN=True" ) GREEN = True;
        else if ( s == "GREEN=False" ) GREEN = False;
        else if ( s == "G=True" ) GREEN = True;
        else if ( s == "G=False" ) GREEN = False;
        else if ( s == "BNG=True" ) BNG = True;
        else if ( s == "BNG=False" ) BNG = False;

        else if ( s.Contains( "PURPLE=", 0 ) ) {
            if ( !engine.HasCounts( ) ) {
                out << "PURPLE only works if there is a file a.counts" << std::endl;
                return False;
            }
            String y = s.After( "PURPLE=" );
            Bool OK = True;
            for ( int i = 0; i < y.isize( ); i++ )
                if ( y[i] != '0' && y[i] != '1' ) OK = False;
            if ( y.size( ) != engine.Count( ).size( ) ) OK = False;
            if ( !OK ) {
                out << "Your PURPLE value doesn't make sense." << std::endl;
                return False;
            }
            PURPLE = y;
        }

        else if ( s == "NEATO=True" ) NEATO = True;
        else if ( s == "NEATO=False" ) NEATO = False;
        else if ( s == "SVG=True" ) SVG = True;
        else if ( s == "SVG=False" ) SVG = False;
        else if ( s == "N=True" ) NEATO = True;
        else if ( s == "N=False" ) NEATO = False;
        else if ( s == "REL=True" ) REL = True;
        else if ( s == "REL=False" ) REL = False;
        else if ( s == "LINENO=True" ) LINENO = True;
        else if ( s == "LINENO=False" ) LINENO = False;
        else if ( s == "SHOW_INV=True" ) SHOW_INV = True;
        else if ( s == "SHOW_INV=False" ) SHOW_INV = False;
        else if ( s == "SHOW_CN=True" ) SHOW_CN = True;
        else if ( s == "SHOW_CN=False" ) SHOW_CN = False;
        else if ( s == "SHOW_ALIGN=True" ) SHOW_ALIGN = True;
        else if ( s == "SHOW_ALIGN=False" ) SHOW_ALIGN = False;
        else if ( s == "ALTREF=True" ) ALTREF = True;
        else if ( s == "ALTREF=False" ) ALTREF = False;
        else if ( s.Contains( "FONTSIZE=", 0 )
                  && s.After( "FONTSIZE=" ).IsDouble( )
                  && s.After( "FONTSIZE=" ).Double( ) > 0 ) {
            FONTSIZE = s.After( "FONTSIZE=" ).Double( );
        } else if ( s.Contains( "ASPECT=", 0 )
                    && s.After( "ASPECT=" ).IsDouble( )
                    && s.After( "ASPECT=" ).Double( ) > 0 ) {
            ASPECT = s.After( "ASPECT=" ).Double( );
        } else if ( s.Contains( "FS=", 0 ) && s.After( "FS=" ).IsDouble( )
                    && s.After( "FS=" ).Double( ) > 0 ) {
            FONTSIZE = s.After( "FS=" ).Double( );
        } else if ( s.Contains( "SCALE=", 0 )
                    && s.After( "SCALE=" ).IsDouble( )
                    && s.After( "SCALE=" ).Double( ) > 0 ) {
            SCALE = s.After( "SCALE=" ).Double( );
        } else {
            out << "I don't grok " << s << "." << std::endl;
            return False;
        }
    }
    return True;
}

// Note that the values below are mirrored in NhoodInfoCore.cc.

void nhood_info_state::Initialize( ) {
    SEEDS = "";
    DEPTH = 5;
    EXT = False;
    COUNT = False;
    GREEN = False;
    BNG = False;
    PURPLE = "";
    NEATO = False;
    REL = False;
    LINENO = False;
    SHOW_INV = False;
    SHOW_CN = True;
    SHOW_ALIGN = True;
    ALTREF = False;
    FONTSIZE = 20;
    SCALE = 1.3;
    COV2 = False;
    ASPECT = -1;
    SVG = False;
}
