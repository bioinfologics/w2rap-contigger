///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

const char *DOC =
    "CrossOut finds and removes cross-contamination from DISCOVAR de novo assemblies "
    "in cases where samples have been prepared in parallel and in close proximity in "
    "the lab. It cannot be used with a single assembly to remove general contamination.";

// ISSUES
// 1. Only tested on one case.
// 2. Ignores reverse complements: support for a kmer is measured by the number of
//    times the kmer occurs in the reads.
// 3. Slow and uses a lot of memory.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "MainTools.h"
#include "ParallelVecUtilities.h"
#include "kmers/MakeLookup.h"
#include "math/Functions.h"
#include "paths/HyperBasevector.h"
#include "paths/RemodelGapTools.h"
#include "paths/long/ReadPath.h"
#include "paths/long/large/FinalFiles.h"
#include "paths/long/large/GapToyTools.h"
#include "paths/long/large/Lines.h"
#include "bam/ReadBAM.h"

int main(int argc, char *argv[]) {
    RunTime( );

    BeginCommandArguments;
    CommandDoc(DOC);
    CommandArgument_String_Doc(DIR, "parent directory for assemblies");
    EndCommandArguments;

    // Find assembly directories.

    vec<String> all_files = AllFiles(DIR), dirs;
    for ( int i = 0; i < all_files.isize( ); i++ ) {
        if ( IsDirectory( DIR + "/" + all_files[i] ) )
            dirs.push_back( all_files[i] );
    }
    int nsam = dirs.size( );
    cout << Date( ) << ": loading" << std::endl;
    const int verbosity = 1; // 0 or 1 or 2

    // Heuristics.

    const int K = 60;
    const double min_ratio = 20;
    const double min_tag = 0.1;

    // Globals, not sure where these should go.

    const Bool ALIGN_TO_GENOME = False;
    const String EVALUATE = "False";
    const Bool EVALUATE_VERBOSE = False;
    const String X;
    map<String,GapToyResults> res;
    const String SAMPLE;
    const String species;
    const vec<int> fosmids;
    const vecbasevector G;

    // Load data.

    vecbasevector all;
    vec<int64_t> starts( nsam + 1, 0 ), starts2( nsam + 1, 0 );
    for ( int j = 0; j < nsam; j++ ) {
        String work_dir = DIR + "/" + dirs[j];
        all.ReadAll( work_dir + "/data/frag_reads_orig.fastb", True );
        starts[j+1] = all.size( );
    }

    starts2[0] = all.size( );
    int L = -1;
    for ( int j = 0; j < nsam; j++ ) {
        String work_dir = DIR + "/" + dirs[j];
        const String final_dir = work_dir + "/a.final";
        HyperBasevectorX hb;
        BinaryReader::readFile( final_dir + "/a.hbx", &hb );
        L = hb.K( );
        for ( int e = 0; e < hb.E( ); e++ )
            all.push_back( hb.EdgeObject(e) );
        starts2[j+1] = all.size( );
    }

    // Compute total number of kmers in each dataset.

    vec<int64_t> total( nsam, 0 );
    for ( int j = 0; j < nsam; j++ )
        for ( int64_t i = starts[j]; i < starts[j+1]; i++ )
            if ( all[i].isize( ) >= K ) total[j] += all[i].isize( ) - K + 1;

    // Find tags.

    vec< triple<int,int,int> > tags;
    vec<String> prefixes = { "A", "C", "G", "T" };
    for ( int pass = 0; pass < 4; pass++ ) {
        cout << Date( ) << ": making kmers, pass " << pass+1 << std::endl;
        vec< triple<kmer<K>,int,int> > kmers_plus;
        MakeKmerLookup0Pre( all, prefixes[pass], kmers_plus );
        const int64_t batches = 100;
        vec<int64_t> bstart(batches+1);
        for ( int64_t i = 0; i <= batches; i++ )
            bstart[i] = ( (int64_t) kmers_plus.size( ) * i ) / batches;
        for ( int64_t i = 1; i < batches; i++ ) {
            int64_t& s = bstart[i];
            while( s > bstart[i-1]
                    && kmers_plus[s].first == kmers_plus[s-1].first ) {
                s--;
            }
        }
        cout << Date( ) << ": traversing kmers" << std::endl;
        #pragma omp parallel for
        for ( int64_t bi = 0; bi < batches; bi++ ) {
            vec<int> count(nsam);
            vec<Bool> assembled(nsam), tagged(nsam);
            vec<double> molarity(nsam);
            for ( int64_t i = bstart[bi]; i < bstart[bi+1]; i++ ) {
                int64_t j;
                for ( j = i + 1; j < bstart[bi+1]; j++ )
                    if ( kmers_plus[j].first != kmers_plus[i].first ) break;
                if ( j - i == 1 ) continue;
                for ( int l = 0; l < nsam; l++ ) {
                    count[l] = 0;
                    assembled[l] = False;
                    tagged[l] = False;
                    molarity[l] = 0;
                }
                {
                    for ( int64_t k = i; k < j; k++ ) {
                        int id = kmers_plus[k].second, p;
                        for ( p = 0; p < nsam; p++ )
                            if ( id < starts[p+1] ) break;
                        if ( p < nsam ) count[p]++;
                        if ( id >= starts2[0] ) {
                            for ( p = 0; p < nsam; p++ )
                                if ( id < starts2[p+1] ) break;
                            assembled[p] = True;
                        }
                    }
                }
                for ( int l = 0; l < nsam; l++ ) {
                    if ( count[l] > 0 )
                        molarity[l] = double(count[l]) / double(total[l]);
                }
                double high = Max(molarity);
                for ( int l = 0; l < nsam; l++ ) {
                    if ( count[l] > 0 && high/molarity[l] >= min_ratio
                            && assembled[l] ) {
                        tagged[l] = True;
                    }
                }
                if ( Sum(tagged) > 0 ) {
                    #pragma omp critical
                    {
                        if ( verbosity >= 2 ) {
                            cout << kmers_plus[i].first.ToString( );
                            for ( int l = 0; l < count.isize( ); l++ ) {
                                if ( count[l] > 0 ) {
                                    cout << " " << l+1 << "["
                                    << count[l] << "]";
                                }
                            }
                            cout << std::endl;
                        }
                        for ( int64_t k = i; k < j; k++ ) {
                            int id = kmers_plus[k].second, p;
                            if ( id >= starts2[0] ) {
                                for ( p = 0; p < nsam; p++ )
                                    if ( id < starts2[p+1] ) break;
                                if ( tagged[p] ) {
                                    tags.push( p, id - starts2[p],
                                               kmers_plus[k].third );
                                    if ( verbosity >= 2 ) {
                                        cout << "tagging " << p+1 << "."
                                             << id - starts2[p] << "."
                                             << kmers_plus[k].third
                                             << std::endl;
                                    }
                                }
                                assembled[p] = True;
                            }
                        }
                    }
                }
                i = j - 1;
            }
        }
    }

    // Collate tags.

    ParallelSort(tags);
    vec<vec<int>> dels(nsam);
    cout << "\ntag summary\n";
    for ( int64_t i = 0; i < tags.jsize( ); i++ ) {
        int64_t j;
        int id = tags[i].first, e = tags[i].second;
        for ( j = i + 1; j < tags.jsize( ); j++ )
            if ( tags[j].first != id || tags[j].second != e ) break;
        int len = all[ starts2[id] + e ].size( );
        int total = len - K + 1 - 2 * (L-K), hits = 0;
        for ( int64_t k = i; k < j; k++ )
            if ( tags[k].third >= L-K && tags[k].third <= len - L ) hits++;
        if ( total > 0 && double(hits)/double(total) >= min_tag ) {
            dels[id].push_back(e);
            cout << id + 1 << "." << e << ".";
            for ( int64_t k = i; k < j; k++ ) {
                int64_t l;
                for ( l = k + 1; l < j; l++ )
                    if ( tags[l].third != tags[l-1].third + 1 ) break;
                if ( k > i ) std::cout << ",";
                cout << tags[k].third;
                if ( l-1 > k ) std::cout << "-" << tags[l-1].third;
                k = l - 1;
            }
            cout << " (" << PERCENT_RATIO( 3, hits, total ) << ")\n";
        }
        i = j - 1;
    }
    cout << "\n" << Date( ) << ": tagging complete" << std::endl;

    // Delete edges and write final files.

    cout << Date( ) << ": deleting edges and writing files" << std::endl;
    for ( int j = 0; j < nsam; j++ ) {
        if ( dels[j].empty( ) ) continue;
        String work_dir = DIR + "/" + dirs[j];
        const String final_dir = work_dir + "/a.final";
        const String clean_dir = work_dir + "/a.clean";
        Mkdir777(clean_dir);

        // Delete edges.

        HyperBasevector hb;
        BinaryReader::readFile( final_dir + "/a.hbv", &hb );
        vec<int> inv;
        BinaryReader::readFile( final_dir + "/a.inv", &inv );
        ReadPathVec paths( final_dir + "/a.paths" );
        int nd = dels[j].size( );
        for ( int i = 0; i < nd; i++ )
            dels[j].push_back( inv[ dels[j][i] ] );
        hb.DeleteEdges( dels[j] );
        RemoveSmallComponents3( hb, True );
        Cleanup( hb, inv, paths );

        // Write final files.

        const int MAX_CELL_PATHS = 50; // note copied from GapToyCore.cc
        const int MAX_DEPTH = 10;      // note copied from GapToyCore.cc
        vec<String> subsam_names;
        vec<int64_t> subsam_starts;
        BinaryReader::readFile( work_dir + "/subsam.names", &subsam_names );
        BinaryReader::readFile( work_dir + "/subsam.starts", &subsam_starts );

        Bool SAVE_FASTA = True;
        FinalFiles( hb, inv, paths, subsam_names, subsam_starts, work_dir,
                    clean_dir, MAX_CELL_PATHS, MAX_DEPTH, ALIGN_TO_GENOME, EVALUATE,
                    EVALUATE_VERBOSE, X, res, SAMPLE, species, fosmids,
                    G, SAVE_FASTA );
    }

    // Done.

    cout << "\n" << Date( ) << ": DONE" << std::endl;
    Scram(0);
}
