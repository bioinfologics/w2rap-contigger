///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <omp.h>

#include <algorithm>

#include "CoreTools.h"
#include "PrintAlignment.h"
#include "pairwise_aligners/SmithWatBandedA.h"
#include "pairwise_aligners/SmithWatFree.h"
#include "paths/long/CreateGenome.h"
//#include "paths/long/EvalCorrected.h"
#include "paths/long/Friends.h"
#include "paths/long/KmerAlign.h"
#include "paths/long/LongProtoTools.h"
#include "paths/long/MakeKmerStuff.h"
#include "paths/long/OverlapReads.h"
#include "paths/long/ultra/Prefab.h"
#include "random/Random.h"

//#define PREFAB_EXT_PARA            // allow external parameters to be fed in

void SelectInitialReads( const vec<int>& rid, const IAndOsVec& F, vec<int>& rid1,
     const long_heuristics& heur, const long_logging_control& log_control )
{    srandomx(1234567);
     for ( int i = 0; i < rid.isize( ); i++ )
          if ( randomx( ) % 5 == 0 ) rid1.push_back( rid[i] );    }


namespace {  // anonymous namespace for local functions

    // Local data structure to store the mapping of the founder reads to candidate assemblies
    struct AlignRange {
        align a;
        int pos1, Pos1; // positin on the read
        int pos2, Pos2; // position on the assembly
        int n_errors;   // total number of errors
        float error_rate;   // error rate
    };
    bool operator< (const AlignRange& lhs, const AlignRange& rhs) { return lhs.error_rate < rhs.error_rate; }
    std::ostream& operator<< ( std::ostream& out, const AlignRange& r) {
        return out << "[" << r.pos1 << ", " << r.Pos1 << ") to [" << r.pos2 << ", " << r.Pos2 << ") n_error= " 
                   << r.n_errors << "(" << r.error_rate << ")";   }


    // Select the seemingly best assembly based on the alignment results
    // ,which is the first qualifying assembly it finds after sort all 
    // assemblies by error rates.
    int SelectAssembly( const vec<AlignRange>& results, double Max_Err_Rate,
           int Min_Assembly_Lengh, const vec<basevector>& assemblies,
           vec<int>* p_alt_sels = 0 ) 
    {
        if ( p_alt_sels != 0 ) p_alt_sels->clear();
        vec<AlignRange> rr = results;
        vec<int> index( results.size(), vec<int>::IDENTITY );
        vec<int> lens( assemblies.size() );
        for ( size_t i = 0; i < assemblies.size(); ++i ) lens[i] = assemblies[i].size();
        SortSync(rr, index, lens);
        int sel = -1;
        for ( size_t i = 0; i < rr.size(); ++i ) {
            //if ( rr[i].pos2 == 0 || rr[i].Pos2 == lens[i] ) continue; // require full range mapping
            if ( rr[i].error_rate > Max_Err_Rate ) continue;
            if ( (rr[i].Pos1 - rr[i].pos1) < Min_Assembly_Lengh ) continue;
            if ( sel == -1 ) sel = index[i];       
            if ( p_alt_sels != 0 ) p_alt_sels->push_back( index[i] );
        }
        return sel;
    }

    // location of a read on the genome
    struct ReadGLoc {
        size_t start, stop, gid;
    };

    ReadGLoc SimRange( const long_logging_control& log_control, int id ) {
        ReadGLoc ans = { (size_t) (*log_control.readlocs)[id].start, 
               (size_t) (*log_control.readlocs)[id].stop, 
               (size_t) (*log_control.readlocs)[id].id };
        return ans;
    }
    basevector ReadTruth( const long_logging_control& log_control, int id) {
        ReadGLoc loc = SimRange( log_control, id );
        return basevector( (*log_control.G)[loc.gid], loc.start, loc.stop - loc.start );
    }
    std::ostream& operator<< ( std::ostream& out, const ReadGLoc& gloc ) {
        return out << "[" << gloc.start << "," << gloc.stop << ")@" << gloc.gid;
    }
    void PrintPathGLocs( std::ostream& out, const vec<ReadOverlap>& path, 
            const vec<int>& creads_filt_ids, const long_logging_control& log_control ) {
        for ( size_t i = 0; i < path.size(); ++i ) {
            if ( i==0 ) 
                out <<  SimRange(log_control,  creads_filt_ids[ path[i].rid1 ] ) << " ";
            out <<  SimRange(log_control, creads_filt_ids[ path[i].rid2 ] ) << " ";
        }
    }

     bool IsTrueFriend( const long_logging_control& log_control, int id1, int id2 )
     {    const vec<ref_loc>& locs = *log_control.readlocs;
          return ( locs[id1].id == locs[id2].id
               && IntervalOverlap( locs[id1].start, locs[id1].stop,
                    locs[id2].start, locs[id2].stop ) > 0 );    }

     // Quickly check if there are errors in the read

    int HasError( const long_logging_control& log_control, int id, const vec<basevector>& reads ) {
        basevector read_truth = ReadTruth( log_control, id );
        for ( size_t i = 0; i < reads.size(); ++i ) 
            if ( search( read_truth.begin(), read_truth.end(), reads[i].begin(), reads[i].end() ) != read_truth.end() ) 
                return 0;
        return 1;
    }

    // Banded Smith-Waterman alignment 
    AlignRange SWAlign( const basevector& s, const basevector& t) 
    {
        align a;
        int error = 0;
        int extra = s.size() / 10;
        int off_ub = 0 + extra;
        int off_lb = - ( (int)t.size() - (int)s.size() ) - extra;
        int offset = ( off_ub + off_lb ) / 2;
        int bandwidth = ( off_ub - off_lb ) / 2;
        SmithWatBandedA( s, t, offset, bandwidth, a, error, 0, 1, 1 );
        float error_rate = error * 1.0 / ( a.Pos1() - a.pos1() );
        AlignRange result = { a, a.pos1(), a.Pos1(), a.pos2(), a.Pos2(), error, error_rate };
        return result;
    }

    void GetKmersAndSort( const int K,  const basevector& u, 
                          vec<basevector>& kmers, vec<int> & locs )
    {   kmers.resize( u.size() - K + 1 );
        locs.resize( u.size() - K + 1 );
        for ( int j = 0; j <= u.isize( ) - K; j++ )
        {   kmers[j].SetToSubOf( u, j, K );
            locs[j] = j;  } 
        SortSync(kmers, locs);  }

    // Banded Smith-Waterman alignment assisted by kmer align ( for offset and bandwidth calculation )
    AlignRange KmerAlign( const int K, const basevector& s, const basevector t, 
                          const vec<basevector>& kmers, const vec<int> & locs )
    {
        vec< std::pair<int,int> > offsets;
        for ( size_t i = 0; i < t.size()-K+1; ++i ) {
            basevector tk( t, i, K);
            vec<basevector>::const_iterator it 
                = lower_bound( kmers.begin(), kmers.end(), tk );
            while ( it != kmers.end() && *it == tk ) {
                offsets.push( locs[it - kmers.begin()] - i, 
                             locs[it - kmers.begin()] );
                ++it;
            }
        }
        AlignRange result0 = { align(), 0, 0, 0, 0, 100, 100.0 };
        if ( offsets.empty() ) return result0;
        vec< std::pair<int,int> > offsets2;
        LargestOffsetCluster( offsets, offsets2, 200, 0.2 );
        {  // kmeralign-assisted SmithWaterman
            Sort( offsets2 );
            int off_u = offsets2.back().first;
            int off_l = offsets2.front().first;
            int off = ( off_u + off_l ) / 2;
            int bandwidth = ( off_u - off_l ) / 2;
            align a;
            int error;
            SmithWatBandedA( s, t, off, bandwidth, a, error, 0, 1, 1 );
            float error_rate = error * 1.0 / ( a.Pos1() - a.pos1() );
            AlignRange result = { a, a.pos1(), a.Pos1(), a.pos2(), a.Pos2(), error, error_rate };
            return result;
        }
    }

    // If we have several candidate corrected reads, which one is supported more by all
    // friend reads
    void VoteByReads( const vec<basevector>& candidates, const vec<basevector>& reads, 
         vec<int>& votes ) 
    {
        const int K2=14;
        votes.assign( candidates.size(), 0 );
        // kmer-lookup for each candidates
        vec< vec<basevector> > kmers( candidates.size() ); 
        vec< vec<int> >locs( candidates.size() );
        for ( size_t i = 0; i < candidates.size(); ++i ) 
            GetKmersAndSort( K2, candidates[i], kmers[i], locs[i] );
        for ( size_t i = 0; i < reads.size(); ++i ) {
            vec< AlignRange > results;
            vec<int> index;
            for ( size_t j = 0; j < candidates.size(); ++j ) {
                AlignRange ar = KmerAlign( K2, candidates[j], reads[i],
                        kmers[j], locs[j] );
                results.push_back( ar );
                index.push_back(j);
            }
            SortSync(results, index);
            if ( results[0].error_rate > 0.1 ) continue;
            votes[index[0]]++;
        }
    }

    // Given the gang of raw friend reads (threads), find all suggested
    // edits to the founder read (t)
    int FindEdits( const vec<basevector>& threads, const basevector& t, 
         vec< vec<edit0> >& edits, std::ostream& out ) 
    {   
        edits.clear_and_resize( t.size( ) + 1 );
        const int K2=14;
        vec<basevector> kmers; vec<int> locs;
        GetKmersAndSort( K2, t, kmers, locs );
        out << "Align threads to read " << std::endl;
        int n_good_reads = 0;
        for ( int i = 0; i < threads.isize( ); i++ ) {    
            AlignRange ar = KmerAlign( K2, t, threads[i], kmers, locs );
            out << i << " " << ar << std::endl;
            if ( ar.error_rate > 0.1 ) continue;
            n_good_reads++;
            const align& a = ar.a;
            int p1 = a.pos1( ), p2 = a.pos2( );
            for ( int j = 0; j < a.Nblocks( ); j++ ) {    
                if ( a.Gaps(j) > 0 ) {    
                    basevector b( threads[i], p2, a.Gaps(j) );
                    edits[p1].push( INSERTION, b.ToString( ) );
                    p2 += a.Gaps(j);    
                }
                if ( a.Gaps(j) < 0 ) {
                    edits[p1].push( DELETION, -a.Gaps(j) ); 
                    p1 -= a.Gaps(j);    
                }
                for ( int x = 0; x < a.Lengths(j); x++ ) {
                    if ( threads[i][p2] != t[p1] )
                        edits[p1].push( SUBSTITUTION, (char) as_base( threads[i][p2] ) );
                        p1++, p2++;    
                   }    
              }    
         }
         for ( int p = 0; p <= t.isize( ); p++ )
             Sort( edits[p] );
         return n_good_reads;    
    }
    
    // Given all the edits generated by FindEdits, identify signals suggesting
    // errors in the founder reads, and make the correction.
    int  PFMakeEdits( const vec<basevector>& threads, basevector& t, 
         const vec< vec<edit0> >& edits )
    {    const int Threshold_Vote = threads.isize( ) * 0.5;
         vec< std::pair<int,edit0> > keepers;
         for ( int p = 0; p <= t.isize( ); p++ )
         {    for ( int i = 0; i < edits[p].isize( ); i++ )
              {    int j = edits[p].NextDiff(i);
                   if ( (j-i) > Threshold_Vote ) keepers.push( p, edits[p][i] );
                   i = j - 1;    }    }
         int nedits = 0;
         for ( int j = keepers.isize( ) - 1; j >= 0; j-- )
         {    if ( j < keepers.isize( ) - 1 && keepers[j].first == keepers[j+1].first )
                   continue;
              nedits++;
              int p = keepers[j].first;
              const edit0& e = keepers[j].second;
              if ( e.etype == INSERTION )
              {    basevector b1( t, 0, p );
                   basevector b2( e.seq );
                   basevector b3( t, p, t.size( ) - p );
                   t = Cat( b1, b2, b3 );    }
              if ( e.etype == DELETION )
              {    basevector b1( t, 0, p );
                   basevector b2( t, p + e.n, t.isize( ) - ( p + e.n ) );
                   t = Cat( b1, b2 );    }
              if ( e.etype == SUBSTITUTION ) t.Set( p, as_char( e.seq[0] ) );    
         } 
         return nedits;
    }



    // Sometimes efasta contains {,}, creating redundant edges when expanding. 
    // We want to remove those empty brackets.
    efasta CleanEfasta( const efasta& input) {
        efasta output;
        output.reserve( input.size() );
        efasta::const_iterator it = input.begin();
        efasta::const_iterator it2;
        while( it!= input.end() ) {
            if ( *it == '{' ) {
                // point it2 to '}'
                it2 = it+1;
                bool is_good = false;
                while( *it2 != '}' ) {
                    if ( *it2 != ',' ) is_good = true;
                    it2++;
                }
                if ( is_good ) output.append( it, it2+1 );
                it = it2;
            } 
            else output.push_back(*it);
            ++it;
        }
        return output;
    }
    enum CSMR_MESSAGE { CSMR_OK, CSMR_NO_FRIENDS, CSMR_NO_CR_FRIENDS, CSMR_NO_ASSEMBLY, CSMR_MANY_ASSEMBLY, CSMR_SHORT_ASSEMBLY } ;
}


std::pair<efasta,CSMR_MESSAGE> CorrectOneReadFast( int ec, const vecbasevector& reads, 
        const IAndOsVec& F,  const vec<int>& known_friends,
        const VecEFasta& corrected, const vec<int>& cid,
        const vec<bool>& is_corrected, const vec<int>& mapping,
        const long_heuristics& heur, const long_logging_control& log_control,
        int VERBOSITY,
        std::ostream& out ) 
{
#ifdef PREFAB_EXT_PARA
    extern int Csmr_Edits;
    extern int Known_Friends;
    extern int Cheat_Reads;          // Cheating: replace the read-to-be-corrected with true sequence
    extern int Cheat_Friend_Reads;   // Cheating: replace the pre-corrected reads with true sequence
    extern int Cheat_Friend_List;    // Cheating: the friend list is based their true genome locations
#else
    const int Csmr_Edits = 0;         // Make further edits using piled reads
    const int Known_Friends = 1;      // Use only friends known during pre-correction
    const int Cheat_Reads = 0;        // Cheating: replace the read-to-be-corrected with true sequence
    const int Cheat_Friend_Reads = 0; // Cheating: replace the pre-corrected reads with true sequence
    const int Cheat_Friend_List = 0;  // Cheating: the friend list is based their true genome locations
#endif

    // where is the real sequence of the read

    const vec<ref_loc>& rlocs = *log_control.readlocs;
    int sim_gid = rlocs[ec].id; 
    int sim_start = rlocs[ec].start, sim_stop = rlocs[ec].stop;
    Bool sim_rc = rlocs[ec].rc2;
    basevector read_truth( 
         (*log_control.G)[sim_gid], sim_start, sim_stop - sim_start );
    if ( sim_rc ) read_truth.ReverseComplement();
    // The read to be corrected
    basevector read_input( reads[ec] );
    const basevector& read0 = ( Cheat_Reads ? read_truth : read_input );
    // friends of the read
    vec<int> friends;
    vec<bool> friend_rcs;
    if ( Cheat_Friend_List ) {
        for ( size_t fid = 0; fid < reads.size(); ++fid ){
            if ( (int)fid == ec ) continue;
            if( IsTrueFriend( log_control, ec, fid ) ) {
                friends.push_back(fid);
                Bool sim_rc2 = rlocs[fid].rc2;
                friend_rcs.push_back(sim_rc ^ sim_rc2);
            }
        }
    }
    else {
        for ( size_t i = 0; i < known_friends.size(); ++i ) {
            bool rc = ( known_friends[i] < 0 );
            int fid = ( rc ? ~known_friends[i] : known_friends[i] );
            friends.push_back(fid);
            friend_rcs.push_back(rc);
        }
    }
    if ( VERBOSITY >=2 ) {
        out << "Make corrections using " << friends.size() << " friend reads " << std::endl;
    }
    if ( friends.size() < 2 ) {
        if ( VERBOSITY >=1 ) 
            out << "Cannot find at least two friends" << std::endl;
        return std::make_pair( efasta(), CSMR_NO_FRIENDS );
    }
    // Generate the assembly of the only the corrected friend reads
    vec< vec<basevector> > creads_filtered;           // expanded efasta 
    vec< vec<Ambiguity> > creads_ambs;                // ambiguities from the efasta
    vec< int > creads_filt_ids;
    for ( size_t i = 0; i < friends.size(); ++i ) {
        int fid = friends[i];
        if ( ! is_corrected[ fid ] ) continue;
        creads_filt_ids.push_back( fid );
        efasta cread;
        if ( Cheat_Friend_Reads ) {
            int sim_start2 = rlocs[fid].start, sim_stop2 = rlocs[fid].stop;
            int sim_gid2 = rlocs[fid].id;
            Bool sim_rc2 = rlocs[fid].rc2;
            ForceAssertEq( sim_gid2, sim_gid );
            basevector read_truth2( 
                 (*log_control.G)[sim_gid2],  sim_start2, sim_stop2 - sim_start2 );
            if ( sim_rc2 ) read_truth2.ReverseComplement();
            cread = efasta( read_truth2.ToString() );
        }
        else
            cread = corrected[ mapping[fid] ] ;
        if (  friend_rcs[i] ) cread.ReverseComplement();
        creads_filtered.push_back( vec<basevector>() );
        cread.ExpandTo( creads_filtered.back() );
        fastavector fa;
        creads_ambs.resize( creads_ambs.size() + 1 );
        cread.FlattenTo( fa, creads_ambs.back(), False ); // extract the ambiguities
        //if ( creads_ambs.back().size() > 0 ) {
        //    out << "read " << creads_ambs.size()-1 << std::endl;
        //    out << cread << std::endl;
        //    out << corrected[ mapping[fid] ] << std::endl;
        //}
    }
    // Validation: are the friends correct?
    if ( VERBOSITY >= 1 ) {
        unsigned int n_true_friends = 0;
        for ( size_t i = 0; i < creads_filt_ids.size(); ++i ) 
            if( IsTrueFriend( log_control, ec, creads_filt_ids[i] ) ) n_true_friends++;
        if ( n_true_friends < creads_filt_ids.size() && VERBOSITY >=1 ) {
            out << "Wrong friending found at read " << ec << " length= " << read0.size() << std::endl;
            out << "n_friends= " << creads_filt_ids.size() 
                << ", n_true_friends= " << n_true_friends 
                << std::endl;
        }
        vec<int> n_errors( creads_filt_ids.size(), 0 ), n_ambs( creads_filt_ids.size(), 0 );
        for ( size_t i = 0; i < creads_filt_ids.size(); ++i ) {
            n_errors[i] = HasError( log_control, creads_filt_ids[i], creads_filtered[i] );
            n_ambs[i] = corrected[ mapping[ creads_filt_ids[i] ] ].AmbEventCount();
        }
        int n_error_reads = n_errors.size() - n_errors.CountValue(0);
        if ( n_error_reads > 0) {
            out << n_error_reads << " out of " << n_errors.size() << " reads has errors: ";
            for ( size_t i = 0; i < n_errors.size(); ++i ) 
                if ( n_errors[i] > 0 ) out << i << " ";
            out << std::endl;
        }
        int n_amb_reads = n_ambs.size() - n_ambs.CountValue(0);
        if ( n_amb_reads > 0) {
            out << n_amb_reads << " out of " << n_ambs.size() << " reads has ambiguities: ";
            for ( size_t i = 0; i < n_ambs.size(); ++i ) 
                if ( n_ambs[i] > 0 ) out << i << "(" << n_ambs[i] << ") ";
            out << std::endl;
        }
    }
    // More about the reads
    if ( VERBOSITY >= 3 ) {
        out << "Read0 start " << 0  << " to " << sim_stop - sim_start << "@" << sim_gid << std::endl;
        out << "Find " << creads_filt_ids.size() << " corrected friend reads" << std::endl;
        vec<int> starts;
        vec<int> gids;
        vec<int> index;
        for ( size_t i = 0; i < creads_filt_ids.size(); ++i ) {
            int fid = creads_filt_ids[i];
            int sim_start2 = rlocs[fid].start, sim_stop2 = rlocs[fid].stop;
            int sim_gid2 = rlocs[fid].id;
            Bool sim_rc2 = rlocs[fid].rc2;
            int start = sim_start2 - sim_start;
            starts.push( start );
            index.push( i );
            gids.push( sim_gid2 );
        }
        SortSync(starts, gids, index);
        for ( size_t i = 0; i < starts.size(); ++i ) 
            out << "    read " << index[i] << "_0 " << starts[i]
                << " to " << starts[i] + creads_filtered[ index[i] ][0].size() 
                << "@" << gids[i] << std::endl;
    }
    if ( creads_filtered.size() < 2 ) {
        if ( VERBOSITY >=1 ) 
            out << "Cannot find at least two corrected friends" << std::endl;
        return std::make_pair( efasta(), CSMR_NO_CR_FRIENDS );
    }
    // assembly of the reads
    const int Max_Num_Assembly = 5;
    SomeReadsAssembler assembler( creads_filtered );
    assembler.SetMinOverlap(800);
    assembler.SetMaxAssemblyNum( Max_Num_Assembly );
    assembler.SetCombineAmbReads( creads_ambs );
    assembler.InitGraph();
    assembler.Assembly();
    const vec<basevector>& assemblies = assembler.GetAssemblies();
    if ( assemblies.empty() ) {
        if ( VERBOSITY >=1 ) 
            out << "Cannot generate any assembly" << std::endl;
        return std::make_pair( efasta(), CSMR_NO_ASSEMBLY );
    }
    if ( assemblies.isize() >= Max_Num_Assembly ) {
        if ( VERBOSITY >=1 ) 
            out << "Generated too many assemblies" << std::endl;
        return std::make_pair( efasta(), CSMR_MANY_ASSEMBLY );
    }
    // kmer-lookup of read0
    const int K2=14;
    vec<basevector> kmers; 
    vec<int> locs;
    GetKmersAndSort( K2, read0, kmers, locs );
    vec<basevector> kmers0; vec<int> locs0;
    if ( VERBOSITY >= 1 )
        GetKmersAndSort( K2, read_truth, kmers0, locs0 );
    vec< AlignRange > results;
    vec< AlignRange > results0; // validation
    for ( size_t i = 0; i < assemblies.size(); ++i ) {
        results.push_back( KmerAlign( K2, read0, assemblies[i], kmers, locs ) );
        //results.push_back( SWFAlign( read0, assemblies[i] ) );
        if ( VERBOSITY >= 1 )
            results0.push_back( KmerAlign( K2, read_truth, assemblies[i], kmers0, locs0 ) );
            //results0.push_back( SWFAlign( read_truth, assemblies[i] ) );
    }
    if ( VERBOSITY >=2 ) {
        if ( VERBOSITY >=3  ) 
            assembler.graph_.PrintGraph(out);
        for ( size_t i = 0; i < assemblies.size(); ++i ) {
            out << "assembly " << i << " (" << assemblies[i].size() << ") "
                << assembler.GetAssemblyPath(i) << std::endl;
            out << " / " << results[i]
                << " / " << results0[i] << std::endl; 
        }
    }
    int sel = SelectAssembly( results, 0.1, read0.size() * 0.8, assemblies );
    if ( sel == -1 ) {
        if ( VERBOSITY >= 1 ) 
            out << "Warning: read mapping on assembly is poor" << std::endl;
        return std::make_pair( efasta(), CSMR_SHORT_ASSEMBLY );
    }
    // Validate the selected assembly path: 
    // Do they cover the whole read and do they come from the same chromosome?
    if ( VERBOSITY >=1 ) {
        int start0 = sim_start;
        int end0 = sim_stop;
        vec<ReadOverlap> path = assembler.GetAssemblyPathRegion(sel, results[sel].pos2, results[sel].Pos2 );
        int rid1 = creads_filt_ids[ path.front().rid1 ];
        int rid2 = creads_filt_ids[ path.back().rid2 ];
        int start = (sim_rc ? SimRange(log_control, rid2 ).start
                            : SimRange(log_control, rid1 ).start );
        int end = (sim_rc ? SimRange(log_control, rid1 ).stop
                          : SimRange(log_control, rid2 ).stop );
        bool path_incomplete = ( start > start0 || end < end0 );
        bool path_invalid  = false;
        for ( size_t i = 0; i < path.size(); ++i ) {
            int rid1 = creads_filt_ids[ path[i].rid1 ];
            int rid2 = creads_filt_ids[ path[i].rid2 ];
            if ( i==0 && (int) SimRange( log_control, rid1 ).gid != sim_gid
                      || (int) SimRange( log_control, rid2 ).gid != sim_gid ) {
                path_invalid = true;
                break;
            }
        }
        if ( path_invalid || path_incomplete ) {
            out << "Read0 " << SimRange(log_control, ec) << std::endl;
            String ee = "";
            if ( path_incomplete ) ee += "incomplete";
            if ( path_invalid ) {
                if ( ! ee.empty() ) ee += "|";
                ee += "invalid";
            }
            out << "Full assembly " << sel << " size " << assemblies[sel].size() 
                << " " << ee 
                << " mapped range " << start - start0 << " " << end - start0 << std::endl;
            out << "Map  range: " ;
            PrintPathGLocs( out, path, creads_filt_ids, log_control );
            out << std::endl;
            for ( size_t i = 0; i < assemblies.size(); ++i ) {
                out << "Assembly " << i << " size " << assemblies[i].size() 
                    << " Path: " << assembler.GetAssemblyPath(i) << std::endl;
                out << " Locs: ";
                PrintPathGLocs( out, assembler.GetAssemblyPath(i), creads_filt_ids, log_control );
                out << std::endl;
            }
        }
    }
    const basevector& assembly = assemblies[sel];
    basevector provisional( assembly, results[sel].pos2, 
                           results[sel].Pos2 - results[sel].pos2 );
    if ( VERBOSITY >= 1 ) {
        AlignRange ar = KmerAlign( K2, read_truth, provisional, kmers0, locs0 );
        if ( ar.n_errors > 0 ) {
            out << "Error found in provisional "
                << ar << " from assembly " << sel << std::endl ;
            for ( size_t i = 0; i < assemblies.size(); ++i ) {
                out << "assembly " << i << " (" << assemblies[i].size()
                    << ") /" << results[i]
                    << " /" << results0[i] << std::endl;
            }
            if ( VERBOSITY >=2 ) {
                out << "read_truth: " << read_truth.ToString() << std::endl;
                out << "read_provi: " << provisional.ToString() << std::endl;
                PrintVisualAlignment( true, out, read_truth, provisional, ar.a  );
                out << "Check accuracy of corrected reads" << std::endl;
                vec< std::pair<int,int> > rid_alts;
                vec<ReadOverlap> path = assembler.GetAssemblyPathRegion(sel, results[sel].pos2, results[sel].Pos2 );
                for ( size_t i = 0; i < path.size(); ++i ) {
                    if ( i== 0 )
                        rid_alts.push( path[i].rid1, path[i].alt1 );
                    rid_alts.push( path[i].rid2, path[i].alt2 );
                }
                for ( size_t i = 0; i < rid_alts.size(); ++i ) {
                    basevector query2 = creads_filtered[ rid_alts[i].first ][ rid_alts[i].second ];
                    basevector read_truth2 = ReadTruth( log_control, creads_filt_ids[ rid_alts[i].first ] );
                    AlignRange ar2 = SWAlign(query2,read_truth2);
                    if ( ar2.n_errors > 0 ) {
                        out << "read " << rid_alts[i].first << "_" << rid_alts[i].second << std::endl;
                        PrintVisualAlignment( true, out, query2, read_truth2, ar2.a  );
                    }
                }
            }
        }
    }
    efasta provisional_efasta = assembler.GetEfastaAssembly(sel, results[sel].pos2, results[sel].Pos2);

    if ( Csmr_Edits == 0 )
        return std::make_pair( provisional_efasta, CSMR_OK );

    // Additional editing seems bring no added value to the error correction, given that
    // those errors have survived multiple rounds of error corrections.

    // collecting all the uncorrected friend reads
    vec<basevector> other_reads;
    other_reads.push_back( read0 );
    for ( size_t i = 0; i < friends.size(); ++i ) {
        int fid = friends[i];
        other_reads.push_back( reads[fid] );
        if ( friend_rcs[i] ) 
            other_reads.back().ReverseComplement();
    }
    if ( VERBOSITY >=2 ) 
        out << "Use " << other_reads.size() << " original reads for further correction " << std::endl;
    if ( other_reads.size() < 10 ) {
        if (VERBOSITY>=1 ) out << "Skipped extra editing due to insufficient reads " << other_reads.size() << std::endl;
        return std::make_pair( provisional_efasta, CSMR_OK );
    }

    // Editing method 1
    //    If there are alternative assemblies ( eg. from other chromosome ), use the original
    //    reads to vote and make the selection.
    //
    if ( Csmr_Edits == 1 ) { 
        // find alternative assemblies
        int match_size = results[sel].Pos1 - results[sel].pos1;
        vec<int> alt_sels; // other possible assemblies
        SelectAssembly( results, 0.1, match_size * 0.9, assemblies, &alt_sels );
        vec<basevector> candidates;
        vec<int> index;
        for ( size_t i = 0; i < alt_sels.size(); ++i ) {
            int aid = alt_sels[i];
            candidates.push_back( basevector(assemblies[aid], results[aid].pos2, 
                           results[aid].Pos2 - results[aid].pos2) );
            index.push_back(aid);
        }
        UniqueSortSync(candidates, index);
        if ( candidates.solo() ) return std::make_pair( provisional_efasta, CSMR_OK );

        if ( VERBOSITY >= 2 )
            out << "Check alternative = " << alt_sels.size() << std::endl;
        // voting
        vec<int> votes( candidates.size(), 0 );
        VoteByReads( candidates, other_reads, votes );
        ReverseSortSync(votes, candidates, index);
        if ( VERBOSITY >=2 ) {
            out << "Votes " << std::endl;
            for ( size_t i = 0; i < votes.size(); ++i ) 
                out << "candidate " << i << " " << votes[i] 
                    << ( candidates[i] == provisional ? "*": "" ) 
                    << std::endl;
        }
        if ( VERBOSITY >= 1 ) {
            AlignRange ar = KmerAlign( K2, read_truth, candidates[0], kmers0, locs0 );
            if ( ar.n_errors > 0 ) {
                out << "Error found in corrected " << ar << std::endl;
                if ( VERBOSITY >=2 ) 
                    PrintVisualAlignment( true, out, read_truth, candidates[0], ar.a  );
            }
        }
        int sel2 = index[0];
        efasta edited_efasta = assembler.GetEfastaAssembly(sel2, results[sel2].pos2, results[sel2].Pos2);
        return std::make_pair( edited_efasta, CSMR_OK ); 
    }

    // Editing method 2:
    //    Pile up all the friend reads on the provisional one, and make edits if more than half reads suggest so.
    //
    if ( Csmr_Edits == 2) {
        basevector to_correct = provisional;
        vec< vec<edit0> > edits;
        int n_good_reads = FindEdits( other_reads, to_correct, edits, out );
        if ( VERBOSITY >= 2 ) {
            out << "#good reads" << n_good_reads << std::endl;
            for ( int p = 0; p < edits.isize( ); p++ )
                if (  edits[p].size() >=1 )
                    out << "site " << p << " has " << edits[p].size() << std::endl;
        }
        int nedits = PFMakeEdits( other_reads, to_correct, edits );
        if ( nedits == 0 ) 
            return std::make_pair( provisional_efasta, CSMR_OK );
        if ( VERBOSITY >= 1 && nedits >0 ) {
            out << "Made " << nedits << "edits" << std::endl;
            AlignRange ar = KmerAlign( K2, read_truth, to_correct, kmers0, locs0 );
            if ( ar.n_errors > 0 ) {
                out << "Error found in corrected " << ar << std::endl;
                if ( VERBOSITY >=2 ) 
                    PrintVisualAlignment( true, out, read_truth, to_correct, ar.a  );
            }
        }
        return std::make_pair( efasta(to_correct), CSMR_OK ); 
    }

    return std::make_pair( provisional_efasta, CSMR_OK );
}

