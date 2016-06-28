///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "FastIfstream.h"
#include "FeudalMimic.h"
#include "math/Functions.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "paths/long/large/DiscoStats.h"
#include "paths/long/large/FinalFiles.h"
#include "paths/long/large/GapToyTools.h"
#include "paths/long/large/Lines.h"
#include "random/Random.h"

// Build final assembly files, starting from the results of scaffolding.

void FinalFiles( const HyperBasevector& hb, const vec<int>& inv, const
                 ReadPathVec& paths, const vec<String>& subsam_names,
                 const vec<int64_t>& subsam_starts, const String& work_dir,
                 const String& final_dir, const int MAX_CELL_PATHS, const int MAX_DEPTH,
                 const Bool ALIGN_TO_GENOME,
                 const String& EVALUATE, const Bool& EVALUATE_VERBOSE,
                 const String& X, std::map<String,GapToyResults>& res, const String& SAMPLE,
                 const String& species, const vec<int>& fosmids, const vecbasevector& G,
                 const Bool SAVE_FASTA ) {
    // Write some assembly files.

    BinaryWriter::writeFile( final_dir + "/a.hbv", hb );
    BinaryWriter::writeFile( final_dir + "/a.hbx", HyperBasevectorX(hb) );
    BinaryWriter::writeFile( final_dir + "/a.inv", inv );
    TestInvolution( hb, inv );

    // Align to genome.

    vec< vec< std::pair<int,int> > > hitsx;
    if ( ALIGN_TO_GENOME && IsRegularFile( work_dir + "/genome.fastb" ) ) {
        vecbasevector genome( work_dir + "/genome.fastb" );
        AlignToGenome( hb, inv, genome, hitsx );
        BinaryWriter::writeFile( final_dir + "/a.aligns", hitsx );
    } else Remove( final_dir + "/a.aligns" );
    if ( ALIGN_TO_GENOME && IsRegularFile( work_dir + "/genome.fastb_alt" ) ) {
        vec< vec< std::pair<int,int> > > hitsx_alt;
        vecbasevector genome( work_dir + "/genome.fastb_alt" );
        AlignToGenome( hb, inv, genome, hitsx_alt );
        BinaryWriter::writeFile( final_dir + "/a.aligns_alt", hitsx_alt );
    } else Remove( final_dir + "/a.aligns_alt" );

    // Make scaffold lines.

    vec<vec<vec<vec<int>>>> linesx;
    // std::cout << "finding scaffold lines, mem usage = "
    //      << ToStringAddCommas( MemUsageBytes( ) ) << std::endl;
    FindLines( hb, inv, linesx, MAX_CELL_PATHS, MAX_DEPTH );
    SortLines( linesx, hb, inv );
    BinaryWriter::writeFile( final_dir + "/a.lines", linesx );
    DumpLineFiles( linesx, hb, inv, paths, final_dir );
    {
        vec<vec<covcount>> covsx;
        {
            ComputeCoverage( hb, inv, paths, linesx, subsam_starts, covsx );
            BinaryWriter::writeFile( final_dir + "/a.covs", covsx );
            vec<int> npairsx;
            {
                vec<int> llensx;
                GetLineLengths( hb, linesx, llensx );
                GetLineNpairs( hb, inv, paths, linesx, npairsx );
                BinaryWriter::writeFile(final_dir + "/a.lines.npairs", npairsx);
                WriteLineStats( final_dir + "/a",
                                linesx, llensx, npairsx, covsx );
            }
            // TestLineSymmetry( linesx, inv );

            // Generate scaffolded assembly fasta file and other output files.

            if (SAVE_FASTA) {
                MakeFinalFasta( hb, inv, linesx, npairsx, covsx, hitsx, final_dir,
                                work_dir, ALIGN_TO_GENOME );
            }
        }

        {
            vecbasevector afinal;
            for ( int e = 0; e < hb.E( ); e++ )
                afinal.push_back( hb.EdgeObject(e) );
            afinal.WriteAll( final_dir + "/a.fastb" );
            paths.WriteAll( final_dir + "/a.paths" );
            VecULongVec invPaths;
            invert( paths, invPaths, hb.E( ) );
            invPaths.WriteAll( final_dir + "/a.paths.inv" );
        }

        // Build genome map.

        if ( IsRegularFile( work_dir + "/genome.names" ) && species == "human" ) {
            fast_ifstream in( work_dir + "/genome.names" );
            vecbitvector genome_amb;
            if ( IsRegularFile( work_dir + "/genome.fastamb" ) )
                genome_amb.ReadAll( work_dir + "/genome.fastamb" );
            else {
                vecbasevector genome( work_dir + "/genome.fastb" );
                Mimic( genome, genome_amb );
            }
            vec<String> genome_names;
            String line;
            while(1) {
                getline( in, line );
                if ( in.fail( ) ) break;
                genome_names.push_back(line);
            }
            BuildGenomeMap( hb, inv, linesx, covsx, hitsx,
                            genome_amb, genome_names, final_dir );
        }
    }

    // Report edge stats.

    String edges_kb_report;
    int64_t genome_size = 0, genome_size0 = 0;
    std::cout << "\n============================================================="
              << "=======================" << std::endl;
    if ( IsRegularFile( work_dir + "/genome.fastb" ) ) {
        vecbasevector genome( work_dir + "/genome.fastb" );

        // Get genome size with and without ambiguous bases.

        genome_size = genome.SizeSum( );
        if ( IsRegularFile( work_dir + "/genome.fastamb" ) ) {
            vecbitvector amb( work_dir + "/genome.fastamb" );
            for ( int g = 0; g < (int) amb.size( ); g++ )
                for ( int64_t p = 0; p < (int64_t) amb[g].size( ); p++ )
                    if ( !amb[g][p] ) genome_size0++;
        } else genome_size0 = genome_size;

        double edges_kb = double( hb.E( ) ) * 1000.0 / double(genome_size);
        std::ostringstream eout;
        eout << edges_kb;
        edges_kb_report = ", edges/kb = " + eout.str( );
        PerfStatLogger::log(
            "edges_per_kb",ToString(edges_kb,2), "edges per kb" );
    }
    {
        vec<int> lens;
        for ( int e = 0; e < hb.EdgeObjectCount( ); e++ ) {
            if ( hb.EdgeLengthBases(e) > 0 )
                lens.push_back( hb.EdgeLengthKmers(e) );
            else lens.push_back(0);
        }
        Sort(lens);
        std::ostringstream out1;
        out1 << Mean(lens);
        std::cout << "\nassembly has " << lens.size( ) << " edges of mean length "
                  << out1.str( ) << edges_kb_report << std::endl;
        if ( res[X].Defined( ) ) {
            std::ostringstream out2;
            out2 << res[X].meanlen;
            if ( lens.isize( ) == res[X].nedges && out1.str( ) == out2.str( ) )
                std::cout << "as EXPECTED." << std::endl;
            else {
                std::cout << "Something changed!  We expected " << res[X].nedges
                          << " edges of mean length " << out2.str( )
                          << "." << std::endl;
            }
        }
    }

    // Remove the gaps from the scaffolded assembly and compute line N50.

    HyperBasevector hby(hb);
    {
        vec<int> dels;
        for ( int e = 0; e < hby.EdgeObjectCount( ); e++ )
            if ( hby.EdgeLengthBases(e) == 0 ) dels.push_back(e);
        hby.DeleteEdges(dels);
    }
    {
        vec<vec<vec<vec<int>>>> linesy;
        FindLines( hby, inv, linesy, MAX_CELL_PATHS, MAX_DEPTH );
        const int min_len = 1000;
        int64_t contigN50 = LineN50( hby, linesy, min_len );
        int64_t scaffoldN50 = LineN50( hb, linesx, min_len );
        vec<int> lens;
        GetLineLengths( hb, linesx, lens );
        int64_t total1 = 0, total10 = 0;
        for ( int i = 0; i < lens.isize( ); i++ ) {
            if ( lens[i] >= 1000 ) total1 += lens[i];
            if ( lens[i] >= 10000 ) total10 += lens[i];
        }
        total1 /= 2;
        total10 /= 2;

        // Estimate chimera rate.

        int64_t chim = 0, chim_total = 0;
        const int min_dist = 2000;
        for ( int64_t i = 0; i < (int64_t) paths.size( ); i += 2 ) {
            const ReadPath &p1 = paths[i], &p2 = paths[i+1];
            if ( p1.size( ) == 0 || p2.size( ) == 0 ) continue;
            int e1 = p1[0], e2 = p2[0];
            int start1 = p1.getOffset( ), start2 = p2.getOffset( );
            if ( start1 < 0 || start2 < 0 ) continue;
            if ( hb.Bases(e1) - start1 >= min_dist ) {
                chim_total++;
                if ( e2 != inv[e1] ) chim++;
            } else if ( hb.Bases(e2) - start2 >= min_dist ) {
                chim_total++;
                if ( e1 != inv[e2] ) chim++;
            }
        }

        int64_t nreads = paths.size( );

        disco_stats stats;
        if ( IsRegularFile( work_dir + "/disco_stats" ) )
            BinaryReader::readFile( work_dir + "/disco_stats", &stats );

        Ofstream( sout, final_dir + "/stats" );
        {
            sout << "# assembly statistics" << std::endl;
            sout << "# please see also frags.dist.png in parent directory"
                 << std::endl;
            sout << "contig line N50: " << ToStringAddCommas(contigN50) << std::endl;
            sout << "scaffold line N50: " << ToStringAddCommas(scaffoldN50)
                 << std::endl;
            sout << "total bases in 1 kb+ scaffolds: "
                 << ToStringAddCommas(total1);
            if ( genome_size > 0 ) {
                sout << " (" << PERCENT_RATIO( 3, total1, genome_size0 )
                     << " of genome)";
            }
            sout << std::endl;
            sout << "total bases in 10 kb+ scaffolds: "
                 << ToStringAddCommas(total10);
            if ( genome_size > 0 ) {
                sout << " (" << PERCENT_RATIO( 3, total10, genome_size0 )
                     << " of genome)";
            }
            sout << std::endl;

            if ( IsRegularFile( work_dir + "/disco_stats" ) ) {
                sout << "There are " << ToStringAddCommas(nreads)
                     << " reads of mean length " << std::setiosflags(std::ios::fixed)
                     << std::setprecision(1) << stats.mean_read
                     << std::resetiosflags(std::ios::fixed) << " and mean base quality "
                     << std::setiosflags(std::ios::fixed) << std::setprecision(1)
                     << stats.mean_qual << "." << std::endl;

                // Compute MPL1.

                vecbasevector some_reads( work_dir + "/sample_reads.fastb" );
                vec<int64_t> sample;
                BinaryReader::readFile(
                    work_dir + "/sample_reads.ids", &sample );
                vec<int> prl;
                for ( int i = 0; i < sample.isize( ); i++ ) {
                    int64_t id = sample[i];
                    if ( id % 2 != 0 ) continue;
                    const basevector& b = some_reads[i];
                    const ReadPath& p = paths[id];
                    if ( p.size( ) == 0 ) {
                        prl.push_back(0);
                        continue;
                    }
                    if ( p.size( ) != 1 ) continue;
                    int offset = p.getOffset( );
                    if ( offset < 0 ) continue;
                    int e = p[0];
                    if ( offset + b.isize( ) > hb.Bases(e) ) continue;
                    int j;
                    for ( j = 0; j < b.isize( ); j++ ) {
                        if ( b[j] != hb.EdgeObject(e)[j+offset] ) break;
                    }
                    prl.push_back(j);
                }
                int MPL1 = int( round( Mean(prl) ) );
                if ( prl.nonempty( ) ) {
                    sout << "MPL1 = mean length of first read in pair up to "
                         << "first error = " << MPL1 << std::endl;
                    sout << "(normal range is 175-225 for 250 base reads)"
                         << std::endl;
                }
            }

            // Compute chimera rate.

            sout << "Estimated chimera rate in read pairs "
                 << "(including mismapping) = "
                 << PERCENT_RATIO( 2, chim, chim_total ) << "." << std::endl;
            if ( double(chim)/double(chim_total) >= 0.05 ) {
                sout << "\nWARNING: based on the observed chimera rate, it "
                     << "looks like something may be\nvery wrong "
                     << "with your input data.\n" << std::endl;
            }
            if ( IsRegularFile( work_dir + "/disco_stats" ) ) {
                sout << "genomic read coverage, using 1 kb+ scaffolds for "
                     << "genome size estimate: " << std::setiosflags(std::ios::fixed)
                     << std::setprecision(1)
                     << double(stats.total_bases)/double(total1)
                     << std::resetiosflags(std::ios::fixed) << std::endl;
            }
            if ( G.size( ) > 0 ) {
                Bool concatenate = ( SAMPLE == "NA12878" || SAMPLE == "HCC1954"
                                     || SAMPLE == "HCC1954BL" );
                int n50 = N50PerfectStretch( hby, G, concatenate );
                sout << "N50 perfect stretch = " << n50 << std::endl;
                PerfStatLogger::log(
                    "N50_perfect",n50,"error-free N50");
            }
        }

        fast_ifstream in( final_dir + "/stats" );
        String line;
        while(1) {
            getline( in, line );
            if ( in.fail( ) ) break;
            if ( !line.Contains( "#", 0 ) ) std::cout << line << std::endl;
        }

        PerfStatLogger::log("N50_contig",contigN50,"Contig Line N50.");
        PerfStatLogger::log("N50_scaffold",scaffoldN50,"Scaffold Line N50.");
    }

    // Assay misassemblies.

    if ( ALIGN_TO_GENOME && IsRegularFile( work_dir + "/genome.fastb" ) ) {
        // std::cout << Date( ) << ": assaying misassemblies" << std::endl;
        AssayMisassemblies( hb, inv, hitsx, linesx, final_dir );
    }

    // Compute edge coverage.

    // std::cout << Date( ) << ": computing edge coverage" << std::endl;
    {
        Ofstream( xout, final_dir + "/a.counts" );
        int ns = subsam_names.size( );
        xout << "#";
        for ( int j = 0; j < ns; j++ )
            xout << " " << subsam_names[j];
        xout << "\n";
        vec<vec<int>> count( ns, vec<int>( hb.E( ), 0 ) );
        for ( int64_t id = 0; id < (int64_t) paths.size( ); id++ ) {
            int ss;
            for ( ss = 0; ss < ns; ss++ )
                if ( ss == ns - 1 || id < subsam_starts[ss+1] ) break;
            const ReadPath& p = paths[id];
            for ( int j = 0; j < (int) p.size( ); j++ ) {
                count[ss][ p[j] ]++;
                if ( inv[ p[j] ] != p[j] )
                    count[ss][ inv[ p[j] ] ]++;
            }
        }
        for ( int e = 0; e < hb.E( ); e++ ) {
            xout << e;
            for ( int j = 0; j < ns; j++ )
                xout << " " << count[j][e];
            xout << "\n";
        }
        BinaryWriter::writeFile( final_dir + "/a.countsb", count );
    }

    // Evaluate assembly.

    if ( EVALUATE == "True" ) {
        double clock = WallClockTime( );
        std::cout << "\n";
        GapToyEvaluate( SAMPLE, species, hby, G, fosmids, work_dir,
                        final_dir, 501, res[X], EVALUATE_VERBOSE );
        // std::cout << "\n" << TimeSince(clock) << " used in assembly evaluation"
        //      << std::endl;
    }
}
