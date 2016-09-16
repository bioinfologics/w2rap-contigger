///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "Basevector.h"
#include "CoreTools.h"
#include "PairsManager.h"
#include "ParseSet.h"
#include "Qualvector.h"
#include "lookup/LookAlign.h"
#include "paths/FindErrorsCore.h"
#include "paths/HyperKmerPath.h"
//#include "kmers/KmerSpectrumCore.h"
#include "kmers/naif_kmer/KernelKmerStorer.h"
#include "lookup/LibInfo.h"
#include "lookup/SAM2CRD.h"
#include "paths/MergeReadSetsCore.h"
//#include "paths/long/Correct1.h"
#include "paths/long/Correct1Pre.h"
#include "paths/long/CorrectPairs1.h"
#include "paths/long/DataSpec.h"
#include "paths/long/FillPairs.h"
#include "paths/long/Heuristics.h"
#include "paths/long/LoadCorrectCore.h"
#include "paths/long/Logging.h"
#include "paths/long/LongProtoTools.h"
#include "paths/long/LongReadsToPaths.h"
#include "paths/long/PreCorrectAlt1.h"
#include "paths/long/PreCorrectOldNew.h"
#include "paths/long/DiscovarTools.h"
#include <numeric>
#include <type_traits>

void SamIAm( const int i, const String& getsam, const String& TMP, bool keepLocs, 
     String const& dexterLibs, Bool PF_ONLY, const Bool KEEP_NAMES )
{
    std::ofstream ofs((TMP+"/SAMProcessing"+ToString(i)+".log").c_str());
    ofs << "Running: " << getsam << std::endl;
    Logger logger(ofs);
    procbuf buf(getsam.c_str(),std::ios::in);
    std::istream is(&buf);
    SAM::SAMFile samfile(is,logger);
    vecbasevector seqs;
    vecqualvector quals;
    vec<look_align_x> alns;
    vec<pairinfo> pairsVec;
    vecString namesv;
    vecString libNames;
    bool use_OQ = true;
    vec<Bool> first_in_pair;
    SAM2CRD(samfile, seqs, quals, alns, pairsVec, namesv, first_in_pair, libNames,
                false, true, use_OQ, PF_ONLY, false, KEEP_NAMES );
    if (KEEP_NAMES)
    {    for ( int64_t i = 0; i < (int64_t) namesv.size( ); i++ )
         {    if ( first_in_pair[i] ) namesv[i] += ".1";
              else namesv[i] += ".2";    }    }
    int samtoolsResult = buf.close();
    if ( samtoolsResult )
    {
        ofs << "Error: samtools returned with non-zero status="
                << samtoolsResult << '\n';
        ofs.close();
        DiscovarTools::ExitSamtoolsFailed();
    }

    String outHead = TMP+'/'+ToString(i);
    if ( keepLocs )
    {
        Ofstream( aout, outHead + ".qltout");
        Ofstream( mapout, outHead + ".mapq");
        for ( auto itr=alns.begin(),end=alns.end(); itr != end; ++itr )
        {
            itr->PrintParseable(aout);
            mapout << (int)itr->mapQ() << "\n";
        }
    }

    int const SEP = -15;
    int const DEV = 12;
    int const NOMINAL_READ_LEN = 251;
    PairsManager pairs(seqs.size());
    if ( dexterLibs.empty() )
        for ( auto itr(libNames.begin()),end(libNames.end()); itr!=end; ++itr )
            pairs.addLibrary(SEP, DEV, *itr);
    else
    {
        LibInfoDB libInfoDB(dexterLibs);
        for ( auto itr(libNames.begin()),end(libNames.end()); itr!=end; ++itr )
        {
            LibInfo const* pInfo = libInfoDB.getInfo(*itr);
            if ( pInfo )
            {
                int sep = pInfo->mMean-2*NOMINAL_READ_LEN;
                pairs.addLibrary(sep,pInfo->mStdDev,*itr);
            }
            else
            {
                ofs << "Warning: Cannot find entry for library: "
                        << (*itr) << " in library information file: "
                        << dexterLibs << '\n';
                ofs << "Using default SEP and DEV.\n";
                pairs.addLibrary(SEP, DEV, *itr);
            }
        }
    }

    for ( auto itr(pairsVec.begin()),end(pairsVec.end()); itr != end; ++itr )
        pairs.addPairToLib(itr->readID1,itr->readID2,itr->libraryID,false);

    pairs.Write( outHead + ".pairs" );
    seqs.WriteAll(outHead + ".fastb");
    quals.WriteAll(outHead + ".qualb");
    if (KEEP_NAMES) namesv.WriteAll( outHead + ".names" );
    ofs.close();
}

void MergeReadSets( const vec<String>& heads, const String& TMP,
     const long_logging& logc, const Bool KEEP_NAMES )
{    if ( heads.empty( ) )
     {    std::cout << "\nInternal error, heads is empty." << std::endl;
          std::cout << "Please send us a report.\n" << std::endl;
          _exit(1);    }
     else if ( heads.solo( ) )
     {    System( "cp " + TMP + "/" + heads[0] + ".fastb "
               + TMP + "/frag_reads_orig.fastb" );
          System( "cp " + TMP + "/" + heads[0] + ".qualb "
               + TMP + "/frag_reads_orig.qualb" );
          System( "cp " + TMP + "/" + heads[0] + ".pairs "
               + TMP + "/frag_reads_orig.pairs" );
          if (KEEP_NAMES)
          {    System( "cp " + TMP + "/" + heads[0] + ".names "
                    + TMP + "/frag_reads_orig.names" );    }
          const String& qltout = TMP + "/" + heads[0] + ".qltout";
          if ( logc.KEEP_LOCS ) {
              if ( IsRegularFile( qltout ) )
                  System( "cp " + qltout + " " +
                        TMP + "/frag_reads_orig.qltout" );
              else
                  std::cout << "warning: missing aligns file: " << qltout << std::endl;
          }
     }
     else
     {    String head_out = TMP + "/frag_reads_orig";
          vec<String> headsp( heads.size( ) );
          for ( int i = 0; i < heads.isize( ); i++ )
               headsp[i] = TMP + "/" + heads[i];
          if (CheckFileSetExists( headsp, ".fastb", true ) == false)
               FatalErr("Fastb file missing - see above for details");
          //[GONZA] changed from size_t to uint64_t to match the function definition, it's not used anywhere else
          vec<uint64_t> sizes;
          for ( size_t i = 0; i < heads.size( ); i++ )
               sizes.push_back(MastervecFileObjectCount( headsp[i] + ".fastb") );
          //[GONZA] commented, never used
          //size_t size_sum = BigSum(sizes);
          MergeFeudal( head_out, headsp, ".fastb", sizes );
          MergeFeudal( head_out, headsp, ".qualb", sizes );
          if (KEEP_NAMES) MergeFeudal( head_out, headsp, ".names", sizes );
          PairsManager pairs;
          for ( size_t i = 0; i < heads.size( ); ++i )
          {    PairsManager pairs_loc;
               pairs_loc.Read( headsp[i] + ".pairs" );
               ForceAssertEq( pairs_loc.nReads( ), sizes[i] );
               pairs.Append(pairs_loc);    }
          pairs.Write( head_out + ".pairs" );
          if ( logc.KEEP_LOCS && CheckFileSetExists( headsp, ".qltout", true ) )
               MergeQltout( head_out, headsp, sizes );
      }  }



void PopulateSpecials( const vecbasevector& creads, const PairsManager& pairs,
     const vecbasevector& creads_done, const vec<Bool>& done,
     const VecEFasta& corrected, const int NUM_THREADS, vec<Bool>& special//,
     /*const long_logging& logc*/ )
{
     //if (logc.STATUS_LOGGING) std::cout << Date( ) << ": computing kmers_plus" << std::endl;
     const int M = 40;
     const int min_strong = 5;

     // only implmemented for K<=60 right now; need a wider Kmer below for higher K.
     ForceAssertLe(M, 60);

     typedef Kmer60 Kmer_t;
     typedef KmerKmerFreq<Kmer_t> KmerRec_t;            // Kmer + frequency
     vec<KmerRec_t> kmer_vec;

     // calculate the kmer frequency db, thresholding at min_freq
     Validator valid( min_strong, 0 );
     KernelKmerStorer<KmerRec_t> storer( creads, M, &kmer_vec, &valid );
     naif_kmerize( &storer, NUM_THREADS, false );

     // this is a kludge, but for now we're just going to convert the
     // kmers to the style expected by the code below.
     vec< kmer<M> > kmers;
     for ( size_t i = 0; i < kmer_vec.size(); ++i ) {
         basevector bv(M);
         for ( int j = 0; j < M; ++j )
              bv.Set(j, kmer_vec[i][j] );
         kmer<M> x(bv);
         kmers.push_back(x);
         x.ReverseComplement();
         kmers.push_back(x);
     }
     vec<KmerRec_t>().swap( kmer_vec );                 // return memory for kmer_vec.
     std::sort( kmers.begin(), kmers.end() );           // sort for later lookup, serial to avoid Parallel sort memory penalty

     //if (logc.STATUS_LOGGING)
     //     std::cout << Date( ) << ": computing right extensions" << std::endl;
     const int min_ext = 200;
     vec<Bool> right_ext( kmers.size( ), False );
     //#pragma omp parallel for
     for ( size_t id = 0; id < corrected.size( ); id++ )           // calculating right_ext[i] for kmers[i]
     {    vec<basevector> v;
          corrected[id].ExpandTo(v);
          if ( done[id] ) v.push_back( creads_done[id] );
          for ( int j = 0; j < v.isize( ); j++ )
          {    kmer<M> x;
               for ( int s = 0; s <= v[j].isize( ) - M; s++ )
               {    int ext = v[j].isize( ) - s;
                    x.SetToSubOf( v[j], s );
                    int64_t p;
                    if ( ext >= min_ext )
                    {    p = BinPosition( kmers, x );
                         if ( p >= 0 && !right_ext[p] )
                              right_ext[p] = True;    }
                    ext = s + M;
                    if ( ext >= min_ext )
                    {    x.ReverseComplement( );
                         p = BinPosition( kmers, x );
                         if ( p >= 0 && !right_ext[p] )
                              right_ext[p] = True;    }    }    }    }

     //if (logc.STATUS_LOGGING) std::cout << Date( ) << ": finding specials" << std::endl;
     special.resize( creads.size(), False );
     vec< kmer<M> > fails;
     for ( int64_t i = 0; i < kmers.jsize( ); i++ )
          if ( !right_ext[i] ) fails.push_back( kmers[i] );           // fails is list of kmers for which !right_ext[]
     //#pragma omp parallel for
     for ( int64_t id = 0; id < (int64_t) creads.size( ); id++ )
     {    const int64_t idp = pairs.getPartnerID(id);
          kmer<M> x;
          for ( int s = 0; s <= creads[id].isize( ) - M; s++ )
          {    x.SetToSubOf( creads[id], s );                         // kmerize read and populate "special"
               if ( BinMember( fails, x ) )
               {
                    special[id] = True;
                    special[idp] = True;    }
               if ( s + M >= min_ext )
               {    x.ReverseComplement( );
                    if ( BinMember( fails, x ) )
                    {
                         special[id] = True;
                         special[idp] = True;    }    }    }    }
}

// to switch between VirtualMasterVec<bvec> and vecbasevector. The former cannot be const&
template<typename T>
struct ZeroCorrectedQuals_impl{
    static void do_it(T oreads, vecbvec const& creads, vecqvec* pQuals){
        vecqvec& cquals = *pQuals;
        ForceAssertEq(oreads.size(),creads.size());
        ForceAssertEq(oreads.size(),cquals.size());
        auto iOBV = oreads.begin();
        auto iCBV = creads.begin();
        auto iEnd = cquals.end();
        for ( auto iQV = cquals.begin(); iQV != iEnd; ++iOBV,++iCBV,++iQV )
        {
            bvec const& bvOrig = *iOBV;
            bvec const& bvCorr = *iCBV;
            qvec& qv = *iQV;
            ForceAssertEq(bvOrig.size(), bvCorr.size());
            ForceAssertEq(bvOrig.size(), qv.size());
            auto iO = bvOrig.cbegin();
            auto iC = bvCorr.cbegin();
            for ( auto iQ = qv.begin(), iE = qv.end(); iQ != iE; ++iQ,++iO,++iC )
                if ( *iO != *iC )
                    *iQ = 0;
        }
    }
};
// zero all quality scores associated with corrections (for reads in a file)
void ZeroCorrectedQuals( String const& readsFile, vecbvec const& creads,
                            vecqvec* pQuals )
{
    VirtualMasterVec<bvec> oreads(readsFile);
    ZeroCorrectedQuals_impl<typename std::add_lvalue_reference<decltype(oreads)>::type>::do_it(oreads,creads,pQuals);
}
// zero all quality scores associated with corrections (for reads already in memory)
void ZeroCorrectedQuals( vecbasevector const& oreads, vecbvec const& creads,
                            vecqvec* pQuals )
{
    ZeroCorrectedQuals_impl<decltype(oreads)>::do_it(oreads,creads,pQuals);
}

void CapQualityScores( vecqualvector& cquals, const vec<Bool>& done )
{    const int cap_radius = 4;
     for ( int64_t id = 0; id < (int64_t) cquals.size( ); id++ )
     {    if ( done[id] ) continue;
          vec<int> q( cquals[id].size( ), 1000000 );
          for ( int j = 0; j < (int) cquals[id].size( ); j++ )
          {    int start = Max( 0, j - cap_radius );
               int stop = Min( (int) cquals[id].size( ) - 1, j + cap_radius );
               for ( int l = start; l <= stop; l++ )
                    q[j] = Min( q[j], (int) cquals[id][l] );    }
          for ( int j = 0; j < (int) cquals[id].size( ); j++ )
               cquals[id][j] = q[j];    }    }


void CorrectionSuite(vecbasevector &gbases, vecqualvector &gquals, PairsManager &gpairs,
                     const long_heuristics &heur,
                     vecbasevector &creads,
                     VecEFasta &corrected, vec<int> &cid, vec<pairing_info> &cpartner,
                     const uint NUM_THREADS, const String &EXIT, const double clock,
                     bool useOldLRPMethod,  LongProtoTmpDirManager &tmp_mgr) {
    // Run Correct1.

    vec<int> trace_ids, precorrect_seq;
    //ParseIntSet( logc.TRACE_IDS, trace_ids );

    vec<int> trace_pids;
    //ParseIntSet( logc.TRACE_PIDS, trace_pids );
    for (int i = 0; i < trace_pids.isize(); i++) {
        int64_t pid = trace_pids[i];
        trace_ids.push_back(2 * pid, 2 * pid + 1);
    }

    ParseIntSet("{" + heur.PRECORRECT_SEQ + "}", precorrect_seq, false);
    const int max_freq = heur.FF_MAX_FREQ;

    const String sFragReadsOrig = "frag_reads_orig";
    const String sFragReadsMod0 = "frag_reads_mod0";
    const bool bOrgReadsInMem = tmp_mgr[sFragReadsOrig].inMem();


    double bclock = WallClockTime();
    vecqualvector cquals;
    creads = tmp_mgr[sFragReadsOrig].reads();
    cquals = tmp_mgr[sFragReadsOrig].quals();
    size_t nReads = creads.size();
    ForceAssertEq(nReads, cquals.size());
    size_t nBases = 0, qualSum = 0;
    for (qvec const &qv : cquals) {
        nBases += qv.size();
        qualSum = std::accumulate(qv.begin(), qv.end(), qualSum);
    }
    ForceAssertEq(nBases, creads.SizeSum());

    if (heur.PRECORRECT_ALT1)
        precorrectAlt1(&creads);
    else if (heur.PRECORRECT_OLD_NEW)
        PreCorrectOldNew(&creads, cquals, trace_ids);
    else {
        PC_Params pcp;
        const int K_PC = 25;
        KmerSpectrum kspec(K_PC);
        pre_correct_parallel(pcp, K_PC, &creads, &cquals, &kspec,
                             -1, NUM_THREADS);
    }

    ZeroCorrectedQuals(tmp_mgr[sFragReadsOrig].reads(), creads, &cquals);

    PairsManager const &pairs = tmp_mgr[sFragReadsOrig].pairs();
    pairs.makeCache();

    // Carry out initial pair filling.

    vecbasevector creads_done;
    vec<Bool> to_edit(nReads, True);
    vec<Bool> done(nReads, False);

    double fclock = WallClockTime();
    creads_done = creads;
    vecbasevector filled;
    //if (logc.STATUS_LOGGING)
    //     std::cout << Date( ) << ": start initial pair filling" << std::endl;
    const int MIN_FREQ = 5;
    FillPairs(creads, pairs, MIN_FREQ, filled, heur.FILL_PAIRS_ALT);
    //REPORT_TIME( fclock, "used in FillPairs" );
    double f2clock = WallClockTime();
    int64_t fill_count = 0;
    for (int64_t id = 0; id < (int64_t) filled.size(); id++) {
        if (filled[id].size() == 0) continue;
        fill_count++;
        int n = creads[id].size();
        creads_done[id] = filled[id];
        cquals[id].resize(0);
        cquals[id].resize(filled[id].size(), 40);
        creads[id] = creads_done[id];
        if (n < creads[id].isize()) {
            cquals[id].resize(n);
            if (pairs.getPartnerID(id) >= id) creads[id].resize(n);
            else {
                creads[id].SetToSubOf(creads[id],
                                      creads[id].isize() - n, n);
            }
        }

        done[id] = True;
        if (pairs.getPartnerID(id) < id) { creads_done[id].resize(0); }
        to_edit[id] = False;
    }
    //if (logc.STATUS_LOGGING)
    //{    std::cout << Date( ) << ": "
    //          << PERCENT_RATIO( 3, fill_count, (int64_t) filled.size( ) )
    //          << " of pairs filled" << std::endl;    }
    //REPORT_TIME( f2clock, "used in filling tail" );


    // New precorrection.

    vec<int> trim_to;

    double mclock = WallClockTime();
    //if (logc.STATUS_LOGGING)
    //     std::cout << Date( ) << ": begin new precorrection" << std::endl;

    // Cap quality scores.

    CapQualityScores(cquals, done);
    //REPORT_TIME( mclock, "used capping" );

    // Do precorrection.

    for (int j = 0; j < precorrect_seq.isize(); j++) {
        Correct1Pre(tmp_mgr.dir(), precorrect_seq[j], max_freq, creads, cquals,
                    pairs, to_edit, trim_to, trace_ids, /*logc,*/ heur);
    }
    /*
    Correct1( 40, max_freq, creads, cquals, pairs, to_edit, trim_to,
         trace_ids, log_control, logc );
    */

    double nclock = WallClockTime();

    tmp_mgr[sFragReadsMod0].reads(true) = creads;
    tmp_mgr[sFragReadsMod0].quals(true) = cquals;
    tmp_mgr[sFragReadsMod0].pairs(true);



    // Path the precorrected reads.

    //if (logc.STATUS_LOGGING)
    //     std::cout << Date( ) << ": pathing precorrected reads" << std::endl;
    unsigned const COVERAGE = 50u;
    const int K2 = 80; // SHOULD NOT BE HARDCODED!
    vecbasevector correctedv(creads);
    for (int64_t id = 0; id < (int64_t) creads.size(); id++)
        correctedv[id].resize(trim_to[id]);
    HyperBasevector hb;
    HyperKmerPath h;
    vecKmerPath paths, paths_rc;
    LongReadsToPaths(correctedv, K2, COVERAGE, &hb, &h, &paths, &paths_rc);
    vecKmerPath hpaths;
    vec<tagged_rpint> hpathsdb;
    for (int e = 0; e < h.EdgeObjectCount(); e++)
        hpaths.push_back_reserve(h.EdgeObject(e));
    CreateDatabase(hpaths, hpathsdb);
    //if (logc.STATUS_LOGGING) std::cout << Date( ) << ": done" << std::endl;

    // Close pairs that we're done with.  Code copied with minor
    // changes from LongHyper.cc.  Should be completely rewritten.

    //if (logc.STATUS_LOGGING)
    //     std::cout << Date( ) << ": initially closing pairs" << std::endl;
    //#pragma omp parallel for
    for (int64_t id1 = 0; id1 < (int64_t) creads.size(); id1++) {
        if (done[id1]) continue;
        const int id2 = pairs.getPartnerID(id1);
        if (id2 < id1) continue;
        vec<vec<int> > u(2);
        vec<int> left(2);

        for (int pass = 0; pass < 2; pass++) {
            const KmerPath &p
                    = (pass == 0 ? paths[id1] : paths_rc[id2]);
            vec<triple<ho_interval, int, ho_interval> > M, M2;
            int rpos = 0;
            for (int j = 0; j < p.NSegments(); j++) {
                const KmerPathInterval &I = p.Segment(j);
                vec<longlong> locs;
                Contains(hpathsdb, I, locs);
                for (int l = 0; l < locs.isize(); l++) {
                    const tagged_rpint &t = hpathsdb[locs[l]];
                    int hid = t.PathId();
                    if (hid < 0) continue;
                    longlong hpos = I.Start() - t.Start();
                    longlong start = Max(I.Start(), t.Start());
                    longlong stop = Min(I.Stop(), t.Stop());
                    longlong hstart = start - t.Start();
                    for (int r = 0; r < t.PathPos(); r++)
                        hstart += hpaths[hid].Segment(r).Length();
                    longlong hstop = hstart + stop - start;
                    longlong rstart = rpos + start - I.Start();
                    longlong rstop = rstart + stop - start;
                    M.push(ho_interval(rstart, rstop), hid,
                           ho_interval(hstart, hstop));
                }
                rpos += I.Length();
            }
            Bool bad = False;
            for (int i = 0; i < M.isize(); i++) {
                int j;
                for (j = i + 1; j < M.isize(); j++) {
                    if (M[j].first.Start()
                        != M[j - 1].first.Stop() + 1) { break; }
                    if (M[j].second != M[j - 1].second) break;
                    if (M[j].third.Start()
                        != M[j - 1].third.Stop() + 1) { break; }
                }
                u[pass].push_back(M[i].second);
                Bool incomplete = False;
                if (i > 0 && M[i].third.Start() > 0)
                    incomplete = True;
                if (j < M.isize() && M[j - 1].third.Stop()
                                     != hpaths[M[i].second].KmerCount() - 1) {
                    incomplete = True;
                    bad = True;
                }
                if (i == 0 && j == M.isize() && !incomplete) {
                    i = j - 1;
                    continue;
                }
                int last = (i == 0 ? -1 : M2.back().first.Stop());
                if (M[i].first.Start() > last + 1) bad = True;
                M2.push(ho_interval(M[i].first.Start(),
                                    M[j - 1].first.Stop()), M[i].second,
                        ho_interval(M[i].third.Start(),
                                    M[j - 1].third.Stop()));
                if (j == M.isize() && M[j - 1].first.Stop()
                                      < p.KmerCount() - 1) { bad = True; }
                i = j - 1;
            }
            if (bad) u[pass].clear();
            if (u[pass].nonempty()) { left[pass] = M.front().third.Start(); }
        }
        if (u[0].solo() && u[1].solo() && u[0][0] == u[1][0]) {
            int b1siz = correctedv[id1].isize();
            int b2siz = correctedv[id2].isize();
            int offset = left[1] - left[0];
            if (b1siz == creads[id1].isize()
                && b2siz == creads[id2].isize() && offset >= 0) {
                auto beg = hb.EdgeObject(u[0][0]).begin() + left[0];
                auto end = beg + (left[1] - left[0] + b2siz);
                creads_done[id1].assign(beg, end);
                creads_done[id2] = creads_done[id1];
                creads_done[id2].ReverseComplement();

                creads[id1] = creads_done[id1];
                creads[id1].resize(b1siz);
                creads[id2] = creads_done[id2];

                creads[id2].SetToSubOf(creads[id2],
                                       creads[id2].isize() - b2siz, b2siz);
                cquals[id1].resize(0);
                cquals[id1].resize(creads[id1].size(), 40);
                cquals[id2].resize(0);
                cquals[id2].resize(creads[id2].size(), 40);

                done[id1] = done[id2] = True;
                creads_done[id2].resize(0);
                to_edit[id1] = False;
                to_edit[id2] = False;
            }
        }
    }



    corrected.clear().resize(creads.size());
    CorrectPairs1(tmp_mgr.dir(), 40, max_freq, creads, cquals, pairs, to_edit,
                  trace_ids, heur, /*log_control, logc,*/ corrected);
    for (size_t id = 0; id < corrected.size(); id++) {
        if (corrected[id].size() > 0) {
            to_edit[id] = False;
            const int64_t idp = pairs.getPartnerID(id);
            to_edit[idp] = False;
        }
    }

    if (heur.CP2) {
        double cp2_clock = WallClockTime();
        vec<Bool> special;
        PopulateSpecials(creads, pairs, creads_done, done, corrected,
                         NUM_THREADS, special/*, logc*/ );

        for (size_t id = 0; id < corrected.size(); id++)
            if (!special[id]) to_edit[id] = False;

        long_heuristics heur2(heur);
        // heur2.CP_MIN_GLUE = 5;
        heur2.CP_MIN_GLUE = 15;
        heur2.CP_MINQ_FLOOR = 0;
        heur2.CP_RAISE_ZERO = True;
        heur2.CP_MAX_QDIFF = 25.0;

        CorrectPairs1(tmp_mgr.dir(), 40, max_freq, creads, cquals, pairs, to_edit,
                      trace_ids, heur2, /*log_control, logc,*/ corrected);
    } // end of heur.CP2

    double pclock = WallClockTime();
    for (int64_t id = 0; id < done.jsize(); id++) { if (done[id]) corrected[id] = creads_done[id]; }



}

// Define pairing info.  Note that for now we set all the library ids to 0.

void DefinePairingInfo( const LongProtoTmpDirManager& tmp_mgr, const vecbasevector& creads,
     const vec<Bool>& to_delete, vec<int>& cid, VecEFasta& corrected,
     vec<pairing_info>& cpartner/*, const long_logging& logc*/ )
{    double clock = WallClockTime( );
     PairsManager const& pairs = tmp_mgr.get("frag_reads_orig").pairs();
//     pairs.Read( TMP + "/frag_reads_orig.pairs" );
     for ( int64_t id = 0; id < (int64_t) creads.size( ); id++ )
          if ( !to_delete[id] ) cid.push_back(id);
     corrected.EraseIf(to_delete);
     cpartner.resize( creads.size( ) );
     for ( int64_t xid1 = 0; xid1 < (int64_t) corrected.size( ); xid1++ )
     {    int id1 = cid[xid1];
          int64_t pid = pairs.getPairID(id1);
          if ( pairs.isUnpaired(id1) ) cpartner[xid1] = pairing_info(0,-1,-1);
          else
          {    int id2 = pairs.getPartnerID(id1);
               int xid2 = BinPosition( cid, id2 );
               if ( xid2 < 0 ) cpartner[xid1] = pairing_info(0,-1,-1);
               else
               {    if ( pairs.ID1(pid) == id1 )
                         cpartner[xid1] = pairing_info(1,xid2,0);
                    else cpartner[xid1] = pairing_info(2,xid2,0);    }    }    }
     /*REPORT_TIME( clock, "used in load tail" );*/
}

