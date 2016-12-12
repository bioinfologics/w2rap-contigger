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
#include "FetchReads.h"
#include "Intvector.h"
#include "IteratorRange.h"
#include "PairsManager.h"
#include "ParallelVecUtilities.h"
#include "Qualvector.h"
#include "TokenizeString.h"
#include "efasta/EfastaTools.h"
#include "feudal/HashSet.h"
#include "graph/DigraphTemplate.h"
#include "kmers/KMerHasher.h"
#include "paths/HyperBasevector.h"
#include "paths/ReadsToPathsCoreX.h"
#include "paths/RemodelGapTools.h"
#include "paths/Unipath.h"
//#include "paths/long/EvalAssembly.h"
#include "paths/long/Heuristics.h"
#include "paths/long/LoadCorrectCore.h"
#include "paths/long/Logging.h"
#include "paths/long/LongHyper.h"
#include "paths/long/LongProtoTools.h"
#include "paths/long/MakeKmerStuff.h"
#include "paths/long/ReadPath.h"
#include "paths/long/RefTrace.h"
#include "paths/long/SupportedHyperBasevector.h"
#include "paths/long/large/GapToyTools.h"
#include "random/Random.h"
#include <fstream>
#include <ctime>

PerfStatLogger PerfStatLogger::gInst;

namespace
{

const int MaxIter = 10000;

bool EndsWith( vec<int> const& whole, vec<int> const& e )
{
     if ( whole.size() < e.size() ) return false;
     return std::equal( e.begin(), e.end(), whole.end() - e.size() );
}

bool BeginsWith( vec<int> const& whole, vec<int> const& b )
{
     if ( whole.size() < b.size() ) return false;
     return std::equal( b.begin(), b.end(), whole.begin() );
}

void TailHeadOverlap( vec<int> const& x1, vec<int> const& x2, vec<int>& overlap )
{
     overlap.clear();
     if ( x1.size() == 0 || x2.size() == 0 ) return;

     vec<int> s1(x1), s2(x2);
     UniqueSort(s1);
     UniqueSort(s2);
     if (!Meet(s1,s2)) return;

     int i = x1.size()-1;
     int j = 0;

     // find head of x2 in x1
     for ( ; i >=0; --i )
          if ( x1[i] == x2[j] ) break;

     // now walk forward in x1 and x2 pushing back matches
     for ( ; i >= 0 && i < x1.isize() && j < x2.isize(); ++i, ++j )
          if ( x1[i] == x2[j] ) overlap.push_back( x1[i] );
          else break;

     if ( i < x1.isize() )
          overlap.clear();      // we didn't make it back to the end

     if ( overlap.size() > 1 ) {
          overlap.clear();
     }
}

typedef enum { AC_NOCYCLE, AC_CYCLE, AC_MAXITER } AcReason;

bool AcyclicSubgraphWithLimits( HyperBasevector const& hbv,
          vec<int> const& edge_path,
          vec<int> const& to_left, vec<int> const& to_right,
          const int max_bases, int max_iter = MaxIter,
          const bool verbose = false, AcReason *reasonp = 0 )
{
     ForceAssertGt(edge_path.size(), 0u);
     int K = hbv.K();

     if ( verbose ) std::cout << "edge_path = " << printSeq(edge_path) << std::endl;

     vec<bool> seen(hbv.N(), false);
     vec<bool> special(hbv.N(), false);



     // first validate overlap path
     for ( size_t i = 1; i < edge_path.size(); ++i )
          if ( to_left[edge_path[i]] != to_right[edge_path[i-1]] )
               FatalErr("edge_path is not a path");

     // now mark overlap path as seen
     int len = hbv.K()-1;
     for ( auto e : edge_path ) {
          len += hbv.EdgeLengthKmers(e);
          seen[to_left[e]] = special[to_left[e]] = true;
     }

     // now explore outward from the end
     vec<triple<int,int,int>> stack;
     stack.push( to_right[edge_path.back()], len, edge_path.size() );
     while ( stack.size() ) {
          // we return false if it takes more than max_iter iterations to
          // figure this out
          if ( max_iter-- <= 0 ) {
              if ( reasonp ) *reasonp = AC_MAXITER;
              return false;
          }

          auto rec = stack.back();
          stack.pop_back();
          int rec_v = rec.first;
          int rec_len = rec.second;
          int rec_depth = rec.third;

          if ( verbose ) PRINT3_TO(std::cout, rec_v, rec_len, rec_depth );

          // we *ignore* parts of the graph that are more than max_bases away
          if ( rec_len > max_bases ) continue;

/*
          // We return false if we are able to drift more than rec_depth from
          // our starting point *without* violating the previous rule about length.
          // This is a threshold set for computational reasons -- we return false
          // because we simply can't determine if it's true without exceeding
          // our traversal limit
          if ( rec_depth > max_depth ) return false;
*/

          if ( verbose ) PRINT2_TO(std::cout, seen[rec_v], special[rec_v] );

          // now we now check for a cycle
          if ( seen[ rec_v ] ) {
               // if the cycle hits a vertex in the path, that's bad
               if ( special[ rec_v ] ) {
                   if ( reasonp ) *reasonp = AC_CYCLE;
                   return false;
               }
               // if not, who cares
               else continue;
          }

          seen[rec_v] = true;

          for ( size_t i = 0; i < hbv.From(rec_v).size(); ++i ) {
               int v = hbv.From(rec_v)[i];
               int e = hbv.FromEdgeObj(rec_v)[i];
               stack.push( v, rec_len+hbv.EdgeLengthKmers(e), rec_depth+1 );
          }
     }

     if ( reasonp ) *reasonp = AC_NOCYCLE;
     return true;
}

template <unsigned K>
class IndirectKmer
{
public:
    struct Hasher
    {
        typedef IndirectKmer const& argument_type;
        size_t operator()( IndirectKmer const& kmer ) const
        { return kmer.getHash(); }
    };
    typedef HashSet<IndirectKmer,IndirectKmer::Hasher> Dictionary;

    IndirectKmer() = default;
    IndirectKmer( bvec const* pBV, unsigned offset, size_t hash )
    : mpBV(pBV), mOffset(offset), mHash(hash) {}

    bvec::const_iterator getBases() const { return mpBV->begin(mOffset); }
    size_t getId( bvec const* pBV ) const { return mpBV-pBV; }
    size_t getHash() const { return mHash; }

    static void kmerize( bvec const& gvec, Dictionary* pDict )
    { if ( gvec.size() < K ) return;
      auto end = gvec.end()-K+1;
      auto itr = gvec.begin();
      BuzHasher<K> h;
      uint64_t hVal = h.hash(itr);
      // note: using insertUniqueValue may result in two equal keys being
      // present, but we're actually exploiting that as a feature.
      // (see BizarreComparator, below)
      pDict->insertUniqueValue(IndirectKmer(&gvec,itr.pos(),hVal));
      while ( ++itr != end )
      { hVal = h.step(itr,hVal);
        pDict->insertUniqueValue(IndirectKmer(&gvec,itr.pos(),hVal)); } }

    static bool matches( Dictionary const& dict, bvec const& edge )
    { if ( edge.size() < K ) return false;
      auto end = edge.end()-K+1;
      auto itr = edge.begin();
      BuzHasher<K> h;
      uint64_t hVal = h.hash(itr);
      if ( dict.lookup(IndirectKmer(&edge,itr.pos(),hVal)) ) return true;
      while ( ++itr != end )
      { hVal = h.step(itr,hVal);
        if ( dict.lookup(IndirectKmer(&edge,itr.pos(),hVal)) ) return true; }
      return false; }

    friend bool operator==( IndirectKmer const& k1, IndirectKmer const& k2 )
    { if ( k1.getHash() != k2.getHash() ) return false;
      auto itr = k1.getBases();
      return std::equal(itr,itr+K,k2.getBases()); }

    // Always fails, but accumulates a set of ids associated with the
    // IndirectKmers for some hash value, that are, in fact, equivalent.
    // In other words, it's useful only for an odd side-effect.
    class BizarreComparator
    {
    public:
        BizarreComparator( bvec const* pBV0 ) : mpBV0(pBV0) {}

        bool operator()( IndirectKmer const& k1, IndirectKmer const& k2 ) const
        { if ( k1.getHash() == k2.getHash() )
          { auto itr = k1.getBases();
            if ( std::equal(itr,itr+K,k2.getBases()) )
              mIds.insert(k2.getId(mpBV0)); }
          return false; }

        std::set<size_t> const& getIds() const { return mIds; }

    private:
        bvec const* mpBV0;
        unsigned mK;
        mutable std::set<size_t> mIds;
    };

    static void findMatches( Dictionary const& dict, bvec const* pG0,
                                bvec const& edge, ULongVec& invKeys )
    { if ( edge.size() < K ) return;
      BizarreComparator comp(pG0);
      auto end = edge.end()-K+1;
      auto itr = edge.begin();
      BuzHasher<K> h;
      uint64_t hVal = h.hash(itr);
      dict.lookup(hVal,IndirectKmer(&edge,itr.pos(),hVal),comp);
      while ( ++itr != end )
      { hVal = h.step(itr,hVal);
        dict.lookup(hVal,IndirectKmer(&edge,itr.pos(),hVal),comp); }
      std::set<size_t> const& keys = comp.getIds();
      invKeys.assign(keys.begin(),keys.end()); }

private:
    bvec const* mpBV;
    unsigned mOffset;
    size_t mHash;
};

class Dotter
{
public:
    Dotter( HyperBasevector const& hbv, unsigned adjacencyDepth )
    : mHBV(hbv), mAdjacencyDepth(adjacencyDepth)
    { hbv.ToLeft(mFromVtx); hbv.ToRight(mToVtx); }

    void add( int edgeId )
    { mEdgeMap[edgeId];
      mIncompleteVertexSet.insert(mFromVtx[edgeId]);
      mIncompleteVertexSet.insert(mToVtx[edgeId]); }


private:

    struct EdgeStatus
    {
        EdgeStatus()
        : mUsed(false), mIsAdjacent(false) {}

        void setAdjacent( bool isAdjacent = true ) { mIsAdjacent = isAdjacent; }
        bool isAdjacent() const { return mIsAdjacent; }

        void setUsed( bool used = true ) { mUsed = used; }
        bool isUsed() const { return mUsed; }

        bool mUsed;
        bool mIsAdjacent;
    };


    typedef std::pair<int const,EdgeStatus> MapEntry;



    HyperBasevector const& mHBV;
    vec<int> mFromVtx;
    vec<int> mToVtx;
    unsigned mAdjacencyDepth;
    std::map<int,EdgeStatus> mEdgeMap;
    std::set<int> mIncompleteVertexSet;
};


}; // end of anonymous namespace


void FixPaths(HyperBasevector const& hbv, ReadPathVec& paths)
{
     vec<int> to_right, to_left;
     hbv.ToRight(to_right);
     hbv.ToLeft(to_left);
     #pragma omp parallel for
     for ( int64_t m = 0; m < (int64_t) paths.size( ); m++ )
     {    ReadPath& p = paths[m];
          for ( int i = 0; i < ( (int) p.size( ) ) - 1; i++ )
          {    int e1 = p[i], e2 = p[i+1];
               if ( to_right[e1] != to_left[e2] )
               {    p.resize(i+1);
                    break;    }    }    }
}


void Dot( const int nblobs, int& nprocessed, int& dots_printed, 
     const Bool ANNOUNCE, const int bl )
{
     #pragma omp critical
     {    nprocessed++;
          double done_percent = 100.0 * double(nprocessed) / double(nblobs);
          while ( done_percent >= dots_printed+1)
          {    if ( dots_printed % 10 == 0 && dots_printed > 0 
                    && dots_printed != 50 )
               {    std::cout << " ";    }
                    std::cout << ".";
               dots_printed++;
               if ( dots_printed % 50 == 0 ) std::cout << std::endl;    
               flush(std::cout);    }    }
     if (ANNOUNCE)
     {
          #pragma omp critical
          {    std::cout << "\n" << Date( ) << ": STOP " << bl << std::endl;    }    }    }


void BasesToGraph( vecbasevector& bpathsx, const int K, HyperBasevector& hb )
{    HyperKmerPath h;
     vecKmerPath paths, paths_rc, unipaths;
     vec<big_tagged_rpint> pathsdb, unipathsdb;
     ReadsToPathsCoreY( bpathsx, K, paths );
     CreateDatabase( paths, paths_rc, pathsdb );
     Unipath( paths, paths_rc, pathsdb, unipaths, unipathsdb );
     digraph A; 
     BuildUnipathAdjacencyGraph( paths, paths_rc, pathsdb, unipaths,
          unipathsdb, A );
     BuildUnipathAdjacencyHyperKmerPath( K, A, unipaths, h );
     KmerBaseBrokerBig kbb( K, paths, paths_rc, pathsdb, bpathsx );
     hb = HyperBasevector( h, kbb );    }

void MakeLocalAssembly2(VecEFasta &corrected,
                        const vec<int> &lefts, const vec<int> &rights,
                        SupportedHyperBasevector &shb, const int K2_FLOOR,
                        vecbasevector &creads, vec<int> &cid,
                        vec<pairing_info> &cpartner) {
    long_logging logc("", "");
    logc.STATUS_LOGGING = False;
    logc.MIN_LOGGING = False;
    ref_data ref;
    vec<ref_loc> readlocs;
    long_logging_control log_control(ref, &readlocs, "", "");
    long_heuristics heur("");
    heur.K2_FLOOR = K2_FLOOR;
    int count = 0;
    for (int l = 0; l < (int) corrected.size(); l++)
        if (corrected[l].size() > 0) count++;
    if (count == 0) {
        //mout << "No reads were corrected." << std::endl;
    } else {
        if (!LongHyper(corrected, cpartner, shb, heur, log_control, logc, False)) {
            //mout << "No paths were found." << std::endl;
            SupportedHyperBasevector shb0;
            shb = shb0;
        } else {
            // heur.LC_CAREFUL = True;
            shb.DeleteLowCoverage(heur, log_control, logc);
            if (shb.NPaths() == 0) {
                //mout << "No paths were found." << std::endl;
                SupportedHyperBasevector shb0;
                shb = shb0;
            } //TODO bj, check this is really not needed!
            // else shb.TestValid(logc);
        }
    }
    /*mout << "using K2 = " << shb.K( ) << "\n";
    mout << "local assembly has " << shb.EdgeObjectCount( ) << " edges" << "\n";
    mout << "assembly time 2 = " << TimeSince(clock) << std::endl;*/
}

void LogTime( const double clock, const String& what, const String& work_dir )
{    static String dir;
     if ( work_dir != "" ) dir = work_dir;
     else if ( dir != "" )
     {    Echo( TimeSince(clock) + " used " + what, dir + "/clock.log" );    }    }

void CleanupCore( HyperBasevector& hb, vec<int>& inv, ReadPathVec& paths )
{
     vec<Bool> used;
     hb.Used(used);
     vec<int> to_new_id( used.size( ), -1 );
     {    int count = 0;
          for ( int i = 0; i < used.isize( ); i++ )
               if ( used[i] ) to_new_id[i] = count++;
     }

     vec<int> inv2;
     for ( int i = 0; i < hb.EdgeObjectCount( ); i++ )
     {    if ( !used[i] ) continue;
          if ( inv[i] < 0 ) inv2.push_back(-1);
          else inv2.push_back( to_new_id[ inv[i] ] );
     }
     inv = inv2;

     vec<Bool> to_delete( paths.size( ), False );

     #pragma omp parallel for
     for ( int64_t i = 0; i < (int64_t) paths.size( ); i++ )
     {
          std::vector<int>& p = paths[i];
          for ( int j = 0; j < (int) p.size( ); j++ )
          {    int n = to_new_id[ p[j] ];
               if ( n < 0 ) to_delete[i] = True;
               else p[j] = n;
          }
     }

     hb.RemoveDeadEdgeObjects( );

     hb.RemoveEdgelessVertices( );

     
     }

void Cleanup( HyperBasevector& hb, vec<int>& inv, ReadPathVec& paths )
{    
     //XXX: truncates paths, which should be already done by the graph-modifying functions
    {
        vec<Bool> used;
        hb.Used(used);
        for (int64_t i = 0; i < (int64_t) paths.size(); i++) {
            for (int64_t j = 0; j < (int64_t) paths[i].size(); j++) {
                if (paths[i][j] < 0 || paths[i][j] >= hb.EdgeObjectCount() || !used[paths[i][j]]) {
                    paths[i].resize(j);
                    break;
                }
            }
        }
    }
    RemoveUnneededVertices2( hb, inv, paths );
    CleanupCore( hb, inv, paths );
}

void CleanupLoops( HyperBasevector& hb, vec<int>& inv, ReadPathVec& paths )
{    RemoveUnneededVerticesLoopsOnly( hb, inv, paths );
     CleanupCore( hb, inv, paths );    }
