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
#include "MapReduceEngine.h"
#include "ParallelVecUtilities.h"
#include "Qualvector.h"
#include "Set.h"
#include "feudal/Algorithms.h"
#include "feudal/HashSet.h"
#include "feudal/Mempool.h"
#include "kmers/KMer.h"
#include "kmers/KmerRecord.h"
#include "kmers/MakeLookup.h"
#include "math/Functions.h"
#include "paths/HyperBasevector.h"
#include "paths/RemodelGapTools.h"
#include "paths/long/ReadPath.h"
#include "paths/long/large/GapToyTools.h"
//#include "system/RunTime.h"
#include "system/SpinLockedData.h"
#include "system/WorklistN.h"
#include <algorithm>
#include <ctime>
#include <util/OutputLog.h>

// AlignToGenome.  Find alignments of assembly edges to the genome.  Currently this
// maps edges to pairs (g,p) consisting of a genome contig g and an inferred start
// position p of the edge on g.
//
// First we find 60-mer matches between the assembly and the genome.  We require
// that the 60-mer occur exactly once in the assembly and once in the genome.  Each
// such match is required to exactly perfectly to a total length of 500 bases.



void ReroutePaths( const HyperBasevector& hb, const vec<int>& inv,
     ReadPathVec& paths, const vecbasevector& bases, const VecPQVec& quals )
{
     // Create indices.

     double clock = WallClockTime( );
     vec<int> to_left, to_right;
     hb.ToLeft(to_left), hb.ToRight(to_right);

     // Begin reroute.

     const int max_depth = 3;
     const int max_paths = 200;
     const int max_qsum = 100;

     int improveds = 0;
     #pragma omp parallel for schedule(dynamic, 1000)
     for ( int64_t id = 0; id < (int64_t) paths.size( ); id++ )
     {    ReadPath& p = paths[id];
          // Only consider full placements.

          if ( p.size( ) == 0 ) continue;
          if ( p.getOffset( ) < 0 ) continue;
          vec<int> s( p.size( ) );
          s[0] = p.getOffset( );
          for ( int j = 1; j < (int) p.size( ); j++ )
               s[j] = s[j-1] - hb.EdgeLengthKmers( p[j-1] );
          int n = bases[id].size( );
          if ( s.back( ) + n > hb.EdgeLengthBases( p.back( ) ) ) continue;

          // Find possible starts for the read.

          vec< std::pair<int,int> > starts = { std::make_pair( p[0], p.getOffset( ) ) };
          std::set< std::pair<int,int> > startsx;
          startsx.insert( std::make_pair( p[0], p.getOffset( ) ) );
          vec<int> depth = {0};
          for ( int i = 0; i < starts.isize( ); i++ )
          {    if ( depth[i] == max_depth ) continue;
               int e = starts[i].first, start = starts[i].second;
               int v = to_left[e], w = to_right[e];
               for ( int j = 0; j < hb.To(v).isize( ); j++ )
               {    int ex = hb.EdgeObjectIndexByIndexTo( v, j );
                    int startx = start + hb.EdgeLengthKmers(ex);
                    if ( !Member( startsx, std::make_pair( ex, startx ) ) )
                    {    starts.push( ex, startx );
                         startsx.insert( std::make_pair( ex, startx ) );
                         depth.push_back( depth[i] + 1 );    }    }
               for ( int j = 0; j < hb.From(w).isize( ); j++ )
               {    int ex = hb.EdgeObjectIndexByIndexFrom( w, j );
                    int startx = start - hb.EdgeLengthKmers(e);
                    if ( !Member( startsx, std::make_pair( ex, startx ) ) )
                    {    starts.push( ex, startx );
                         startsx.insert( std::make_pair( ex, startx ) );
                         depth.push_back( depth[i] + 1 );    }    }    }

          // Create initial paths.

          vec<ReadPath> ps;
          for ( int i = 0; i < starts.isize( ); i++ )
          {    if ( starts[i].second < 0 
                    || starts[i].second >= hb.EdgeLengthBases( starts[i].first ) )
               {    continue;    }
               ReadPath q;
               q.push_back( starts[i].first );
               q.setOffset( starts[i].second );
               ps.push_back(q);    }

          // Extend the paths.

          vec<Bool> to_delete( ps.size( ), False );
          for ( int i = 0; i < ps.isize( ); i++ )
          {    if ( i >= max_paths ) break;
               vec<int> s( ps[i].size( ) );
               s[0] = ps[i].getOffset( );
               for ( int j = 1; j < (int) ps[i].size( ); j++ )
                    s[j] = s[j-1] - hb.EdgeLengthKmers( ps[i][j-1] );
               int n = bases[id].size( );
               if ( s.back( ) + n <= hb.EdgeLengthBases( ps[i].back( ) ) ) continue;
               to_delete[i] = True;
               int v = to_right[ ps[i].back( ) ];
               for ( int j = 0; j < hb.From(v).isize( ); j++ )
               {    ReadPath r(ps[i]);
                    r.push_back( hb.EdgeObjectIndexByIndexFrom( v, j ) );
                    ps.push_back(r);    
                    to_delete.push_back(False);    }    }
          if ( ps.isize( ) > max_paths ) continue;
          EraseIf( ps, to_delete );

          // Score the paths.

          const int K = hb.K( );
          vec< std::pair<int,int> > qsum( ps.size( ), std::make_pair(0,0) );
          const basevector& r = bases[id];
          QualVec qv;
          quals[id].unpack(&qv);
          for ( int i = 0; i < ps.isize( ); i++ )
          {    const ReadPath& q = ps[i];
               qsum[i].second = -q.size( );
               int start = q.getOffset( );
               basevector b = hb.EdgeObject( q[0] );
               
               // XXX: Similar opt as in Cleanup (?) 
               for ( int l = 1; l < (int) q.size( ); l++ )
               {    b.resize( b.isize( ) - (K-1) );
                    b = Cat( b, hb.EdgeObject( q[l] ) );    }
               
               /* 
               // Measure total size
               int tot_b = 0;
               for ( int l = 1; l < (int) q.size( ); l++ )
               {  tot_b += b.isize( ) - (K-1); }
               // Allocate total size
               b.resize( tot_b );
               // Put all b togheter in b (using cat??)
               for ( int l = 1; l < (int) q.size( ); l++ )
               {  b.append(hb.EdgeObject( q[l] )); }
               */                       
               for ( int m = 0; m < r.isize( ); m++ )
                    if ( r[m] != b[start+m] ) qsum[i].first += qv[m];    }
          int qorig = qsum[0].first;
          SortSync( qsum, ps );
          for ( int i = 0; i < ps.isize( ); i++ )
               qsum[i].second = -qsum[i].second;
          Bool ok = False;
          for ( int j = 0; j < ps.isize( ); j++ )
          {    if ( ps[j] == p && qsum[j].first == qsum[0].first )
               {    ok = True;    }    }
          if (ok) continue;

          if ( qsum[0].first > max_qsum ) continue;

          #pragma omp critical
          {    improveds++;    }
     
          int ooo = qsum[0].first;
          while( ps.size( ) >= 2 && qsum[0] == qsum[1] )
          {    ps.SetToSubOf( ps, 2, ps.size( ) - 2 );
               qsum.SetToSubOf( qsum, 2, qsum.size( ) - 2 );    }

          vec<Bool> del( ps.size( ), False );
          for ( int j = 1; j < ps.isize( ); j++ )
          {    if ( qsum[j].first > qsum[0].first ) break;
               if ( qsum[j].second < qsum[0].second ) del[j] = True;    }
          EraseIf( ps, del );
          EraseIf( qsum, del );

          if ( ooo < qsum[0].first ) continue;

          Bool verbose = False;
          if (verbose)
          {    std::cout << "\n[" << id << "] ";
               for ( int j = 0; j < ps.isize( ); j++ )
               {    std::cout << ps[j].getOffset( ) << ":" << printSeq(ps[j])
                         << " --> " << qsum[j].first;
                    std::cout << std::endl;
                    break;    }    }
          p = ps[0];
     }

     //std::cout << Date() << ": " << improveds << " / " << paths.size( ) << " paths improved by rerouting" << std::endl;
     LogTime( clock, "rerouting paths" );
     // std::cout << "\n" << Date( ) << ": done" << std::endl;
          }

//    Tamp down hanging ends.
//
//    Under appropriate conditions, replace
//
//    v-----------------e1--------------->w
//    v------------e2--------->x (where no edges emanate from x)
//
//    by
//
//    v--------e1a------------>x----e1b-->w
//    v--------e2a------------>x
//
//    This has the effect of eliminating a hanging end by allowing it to
//    'continue on'.
//
//    Conditions:
//    (a) picture complete except for edges entering v, exiting w
//    (b) e2 matches e1 at end for >= 40 bases
//    (c) e2 mismatches e1 at <= 4 bases
//    (d) |e1| - |e2| + match > K - 1.

void Tamp( HyperBasevector& hb, vec<int>& inv, ReadPathVec& paths,
     const int max_shift )
{
     // Set up.  We build the paths index, but the updates to it are commented
     // out.  Since the bulk of the run time in this module is spent building the
     // paths index, it might be more efficient to pass the paths index to this
     // module, and have it update the index.
     time_t now = time(0);
     //std::cout << "[GapToyTools5.cc] Begining Tamp: " << ctime(&now) << std::endl;
     double clock = WallClockTime( );
     VecULongVec paths_index;
     invert( paths, paths_index, hb.EdgeObjectCount( ) );
     int K = hb.K( );
     vec<int> to_left, to_right;
     hb.ToLeft(to_left), hb.ToRight(to_right);
     auto edges_before=hb.EdgeObjectCount( );
     // Heuristics.

     const int min_match = 40;
     const int max_mismatches = 4;

     // Identify loci to be edited.

     int count = 0;
     vec< triple<int,int,int> > vj;
     vec<int> shift_vj;
     vec<Bool> touched( hb.EdgeObjectCount( ), False );
     for ( int v = 0; v < hb.N( ); v++ )
     {    if ( hb.From(v).size( ) != 2 ) continue;
          for ( int j = 0; j < 2; j++ )
          {    
               // Find candidate e1 and e2.

               int e1 = hb.EdgeObjectIndexByIndexFrom( v, j );
               int e2 = hb.EdgeObjectIndexByIndexFrom( v, 1-j );

               // Test them.

               int n1 = hb.EdgeLengthBases(e1), n2 = hb.EdgeLengthBases(e2);
               if ( n1 <= n2 ) continue;
               int x = to_right[e2], w = to_right[e1];
               if ( !hb.From(x).empty( ) || !hb.To(x).solo( ) ) continue;
               if ( !hb.To(w).solo( ) ) continue;
               const basevector &x1 = hb.EdgeObject(e1), &x2 = hb.EdgeObject(e2);
               if ( !IsUnique( v, x, w ) ) continue;

               int mis = 0, match = 0;
               for ( int l = n2 - 1; l >= 0; l-- )
               {    if ( x1[l] != x2[l] ) 
                    {    mis++;
                         if ( mis > max_mismatches ) break;    }
                    else if ( mis == 0 ) match++;    }
               int shift = 0;
               if ( max_shift == 0 && K - 1 - match < 0 ) continue;
               if ( max_shift == 0 && ( mis > max_mismatches || match < min_match ) )
                    continue;
               if ( max_shift > 0 ) // note does not use max_mismatches!!!!!!!!!!!!!
               {    vec<int> shifts;
                    for ( int s = -max_shift; s <= max_shift; s++ )
                    {    Bool mis = False;
                         for ( int l = n2-1; l > n2-1-min_match; l-- )
                         {    if ( l+s >= x1.isize( ) || x1[l+s] != x2[l] ) 
                              {    mis = True;
                                   break;    }    }
                         if ( !mis ) shifts.push_back(s);    }
                    if ( !shifts.solo( ) ) continue;
                    shift = shifts[0];    
                    match = min_match;    }

               if ( n1 - n2 - shift + match <= K - 1 ) continue;
               int re1 = inv[e1], re2 = inv[e2];
               if ( !IsUnique( e1, e2, re1, re2 ) ) continue; // not sure possible
               if ( touched[e1] || touched[e2] || touched[re1] || touched[re2] )
                    continue;
               touched[e1] = touched[e2] = touched[re1] = touched[re2] = True;
               vj.push( v, j, match );    
               shift_vj.push_back(shift);
               count++;    }    }

     // Edit loci.

     inv.resize( hb.EdgeObjectCount( ) + 4*count );
     // paths_index.reserve( hb.EdgeObjectCount( ) + 4*count );
     for ( int u = 0; u < vj.isize( ); u++ )
     {    int v = vj[u].first, j = vj[u].second, match = vj[u].third;
          int e1 = hb.EdgeObjectIndexByIndexFrom( v, j );
          int e2 = hb.EdgeObjectIndexByIndexFrom( v, 1-j );
          int n1 = hb.EdgeLengthBases(e1), n2 = hb.EdgeLengthBases(e2);
          int x = to_right[e2], w = to_right[e1];
          const basevector &x1 = hb.EdgeObject(e1), &x2 = hb.EdgeObject(e2);
          int re1 = inv[e1], re2 = inv[e2];
          int shift = shift_vj[u];

          // Fix the graph.  Note that we're not checking for interaction between
          // the edit and the rc edit, and this could lead to problems.

          basevector x2a = Cat( x2, basevector( x1, n2 + shift, K - 1 - match ) );
          basevector x1a = basevector( x1, 0, x2a.size( ) + shift );
          basevector x1b = basevector( x1, x2a.isize( ) - (K-1) + shift, 
               x1.isize( ) - ( x2a.isize( ) - (K-1) + shift ) );

          hb.DeleteEdgeFrom( v, j );
          hb.EdgeObjectMutable(e2) = x2a;
          int e1a = hb.AddEdge( v, x, x1a ), e1b = hb.AddEdge( x, w, x1b );

          // Fix the reverse complement instance on the graph.

          x2a.ReverseComplement( );
          x1a.ReverseComplement( ), x1b.ReverseComplement( );
          int rv = to_right[re1], rw = to_left[re1], rx = to_left[re2], rj;
          for ( rj = 0; rj < 2; rj++ )
               if ( hb.EdgeObjectIndexByIndexTo( rv, rj ) == re1 ) break;
          hb.DeleteEdgeTo( rv, rj );
          hb.EdgeObjectMutable(re2) = x2a;
          int re1a = hb.AddEdge( rx, rv, x1a ), re1b = hb.AddEdge( rw, rx, x1b );

          // Update inversion.

          inv[e1] = inv[re1] = -1;
          inv[e1a] = re1a, inv[re1a] = e1a, inv[e1b] = re1b, inv[re1b] = e1b;

          // Update paths and paths_index.

          for ( int64_t l = 0; l < (int64_t) paths_index[e1].size( ); l++ )
          {    ReadPath& p = paths[ paths_index[e1][l] ];
               for ( int m = 0; m < (int) p.size( ); m++ )
               {    if ( p[m] == e1 )
                    {    if ( m > 0 || p.getOffset( ) < hb.EdgeLengthBases(e1a) )
                         {    p[m] = e1a;
                              int p1a = p.getOffset( );
                              for ( int j = 0; j <= m; j++ )
                                   p1a -= hb.EdgeLengthKmers( p[j] );
                              if ( m < (int) p.size( ) - 1 || p1a >= 0 )
                              {    p.insert( p.begin( ) + m + 1 , e1b );
                                   m++;    }    }
                         else
                         {    p[m] = e1b;
                              p.setOffset( p.getOffset( ) 
                                   - hb.EdgeLengthKmers(e1a) );    }    }    }
               // int pi = paths_index[e1][l];
               // paths_index[e1a].push_back(pi);
               // paths_index[e1b].push_back(pi);    
                    }
          for ( int64_t l = 0; l < (int64_t) paths_index[re1].size( ); l++ )
          {    ReadPath& p = paths[ paths_index[re1][l] ];
               for ( int m = 0; m < (int) p.size( ); m++ )
               {    if ( p[m] == re1 )
                    {    if ( m > 0 || p.getOffset( ) < hb.EdgeLengthBases(re1b) )
                         {    p[m] = re1b;
                              int p1b = p.getOffset( );
                              for ( int j = 0; j <= m; j++ )
                                   p1b -= hb.EdgeLengthKmers( p[j] );
                              if ( m < (int) p.size( ) - 1 || p1b >= 0 )
                              {    p.insert( p.begin( ) + m + 1 , re1a );
                                   m++;    }    }
                         else
                         {    p[m] = re1a;
                              p.setOffset( p.getOffset( )
                                   - hb.EdgeLengthKmers(re1b) );    }    }    }
               // int pi = paths_index[re1][l];
               // paths_index[re1a].push_back(pi);
               // paths_index[re1b].push_back(pi);    
                    }
          // paths_index[e1].resize(0), paths_index[re1].resize(0);
          Bool verbose = False;
          if (verbose) PRINT2( e1, e2 );    }

     // Clean up.

     LogTime( clock, "tamping" );
     Cleanup( hb, inv, paths );
     OutputLog(2) << count << " / " << edges_before << " edges tamped, " << hb.EdgeObjectCount() << " edges a after tamping" << std::endl;
     //std::cout << "[GapToyTools5.cc] Finished Tamping: " << ctime(&now) << std::endl;
}





namespace
{
unsigned const KLEN = 28;

size_t findInterestingReadIds( HyperBasevector const& hbv,
                                ReadPathVec const& paths,
                                vecbvec const& reads,
                                vec<size_t>* pReadIds )
{
    vec<int> endEdges;
    if ( true )
    {
        int const MAX_DIST = 10000000; // dangerous!
        int const GOOD_DIST = 500;
        vec<int> D;
        int hbk = hbv.K();
        DistancesToEndFunc( hbv, [hbk](bvec const& bv){return bv.size()-hbk+1;},
                                MAX_DIST, True, D );
        hbv.ToRight(endEdges);
        for ( int& vvv : endEdges )
            vvv = (D[vvv] <= GOOD_DIST);
    }

    size_t nReads = reads.size();
    pReadIds->clear();
    pReadIds->reserve(nReads/200);
    size_t nKmers = 0;
    for ( size_t readId=0; readId != nReads; ++readId )
    {
        if ( paths[readId].empty() )
        {
            ReadPath const& matePath = paths[readId^1];
            if ( !matePath.empty() && endEdges[matePath.back()] )
            {
                bvec const& read = reads[readId];
                if ( read.size() >= KLEN )
                {
                    pReadIds->push_back(readId);
                    nKmers += read.size()-KLEN+1;
                }
            }
        }
    }
    return nKmers;
}

class Loc
{
public:
    Loc()=default;
    Loc( unsigned readId, int offset )
    : mReadId(readId), mOffset(offset) {}

    unsigned getReadId() const { return mReadId; }
    int getOffset() const { return mOffset; }

    void nextOffset() { mOffset += 1; }

    friend bool operator<( Loc const& loc1, Loc const& loc2 )
    { if ( loc1.mReadId != loc2.mReadId ) return loc1.mReadId < loc2.mReadId;
      return loc1.mOffset < loc2.mOffset; }
    friend bool operator==( Loc const& loc1, Loc const& loc2 )
    { return loc1.mReadId == loc2.mReadId && loc1.mOffset == loc2.mOffset; }
    friend std::ostream& operator<<( std::ostream& os, Loc const& loc )
    { return os << loc.getReadId() << '[' << loc.getOffset() << ']'; }

private:
    unsigned mReadId;
    int mOffset;
};

typedef KMer<KLEN> GT_Kmer;
class KmerLoc : public GT_Kmer
{
public:
    KmerLoc()=default;
    template <class Itr>
    KmerLoc( Itr itr, unsigned readId )
    : GT_Kmer(itr), mLoc(readId,0u) {}

    void nextLoc( unsigned char nextBase )
    { toSuccessor(nextBase); mLoc.nextOffset(); }

    Loc const& getLoc() const { return mLoc; }

private:
    Loc mLoc;
};

class KmerLocs : public GT_Kmer
{
public:
    KmerLocs()=default;

    KmerLocs( GT_Kmer const& kmer, Loc* pLocs, Loc* pEnd )
    : GT_Kmer(kmer), mpLocs(pLocs), mNLocs(pEnd-pLocs), mELocs(0)
    { if ( mNLocs != 1 ) std::sort(pLocs,pEnd);
      else
      { ForceAssertGe(sizeof(Loc*),sizeof(Loc));
        *reinterpret_cast<Loc*>(&mpLocs) = *pLocs; } }

    KmerLocs( KmerLoc const& kloc )
    : GT_Kmer(kloc), mNLocs(1u), mELocs(0)
    { ForceAssertGe(sizeof(Loc*),sizeof(Loc));
      *reinterpret_cast<Loc*>(&mpLocs) = kloc.getLoc(); }

    size_t getNLocs() const { return mNLocs; }

    Loc const* locsBegin() const
    { return mNLocs!=1 ? mpLocs : reinterpret_cast<Loc const*>(&mpLocs); }

    Loc const* locsEnd() const
    { return locsBegin()+mNLocs; }

    unsigned getELocs() const { return mELocs; }
    void setELocs( unsigned eLocs ) { mELocs = eLocs; }

    size_t getTotalLocs() const { return mNLocs+mELocs; }

private:
    Loc* mpLocs;
    unsigned mNLocs;
    unsigned mELocs;
};

typedef HashSet<KmerLocs,GT_Kmer::Hasher,std::equal_to<GT_Kmer>> GT_Dict;

class MREReadProc
{
public:
    MREReadProc( vec<size_t> const& ids, vecbvec const& reads, Mempool& alloc,
                    size_t maxMultiplicity, GT_Dict* pDict )
    : mIds(ids), mReads(reads), mAlloc(alloc),
      mMaxMultiplicity(maxMultiplicity), mDict(*pDict),
      mpNext(nullptr), mpEnd(nullptr)
    {}

    template <class OItr>
    void map( size_t idxId, OItr oItr )
    { size_t readId = mIds[idxId];
      bvec const& read = mReads[readId];
      if ( read.size() < KLEN ) return;
      KmerLoc loc(read.begin(),readId);
      *oItr = loc; ++oItr;
      for ( auto itr=read.begin(KLEN),end=read.end(); itr != end; ++itr )
      { loc.nextLoc(*itr); *oItr = loc; ++oItr; } }

    void reduce( KmerLoc const* beg, KmerLoc const* end )
    { size_t nLocs = end-beg;
      if ( nLocs > mMaxMultiplicity )
        return;
      if ( nLocs == 1 )
        mDict.insertUniqueValue(KmerLocs(*beg));
      else
      { Loc* pLocs = alloc(nLocs);
        KmerLocs kLocs(*beg,pLocs,pLocs+nLocs);
        while ( beg != end ) *pLocs++ = beg++->getLoc();
        mDict.insertUniqueValue(kLocs); } }

    KmerLoc* overflow( KmerLoc* beg, KmerLoc* end )
    { return beg+std::min(size_t(end-beg),mMaxMultiplicity+1); }

private:
    Loc* alloc( size_t nnn )
    { if ( mpNext+nnn > mpEnd )
      { size_t nAlloc = std::max(10000ul,nnn);
        mpNext = static_cast<Loc*>(mAlloc.allocate(nAlloc*sizeof(Loc),4));
        mpEnd = mpNext+nAlloc; }
      Loc* result = mpNext; mpNext += nnn; return result; }

    vec<size_t> const& mIds;
    vecbvec const& mReads;
    Mempool& mAlloc;
    size_t mMaxMultiplicity;
    GT_Dict& mDict;
    Loc* mpNext;
    Loc* mpEnd;
};

typedef MapReduceEngine<MREReadProc,KmerLoc,GT_Kmer::Hasher,std::less<GT_Kmer>> RMRE;

class MREEdgeProc
{
public:
    MREEdgeProc( HyperBasevector const& hbv, GT_Dict* pDict )
    : mHBV(hbv), mDict(*pDict)
    {}

    template <class OItr>
    void map( size_t edgeId, OItr oItr )
    { bvec const& edge = mHBV.EdgeObject(edgeId);
      if ( edge.size() < KLEN ) return;
      GT_Kmer kmer(edge.begin());
      *oItr = kmer; ++oItr;
      for ( auto itr=edge.begin(KLEN),end=edge.end(); itr != end; ++itr )
      { kmer.toSuccessor(*itr); *oItr = kmer; ++oItr; } }

    void reduce( GT_Kmer const* beg, GT_Kmer const* end )
    { bumpCount(*beg,end-beg); }

    GT_Kmer* overflow( GT_Kmer* beg, GT_Kmer* end )
    { bumpCount(*beg,end-beg); return beg; }

private:
    void bumpCount( GT_Kmer const& kmer, size_t count )
    { KmerLocs const* pLocs = mDict.lookup(kmer);
      if ( pLocs )
          const_cast<KmerLocs*>(pLocs)->setELocs(pLocs->getELocs()+count); }

    HyperBasevector const& mHBV;
    GT_Dict& mDict;
};

typedef MapReduceEngine<MREEdgeProc,GT_Kmer,GT_Kmer::Hasher> EMRE;

class EdgeProc
{
    static int const WINDOW = 60;
    static int const MAX_MISMATCHES = 4;
    static unsigned char const TRUSTED_QUAL = 30;
public:
    static int const NOT_AN_EDGE = -1;

    EdgeProc( HyperBasevector const& hbv, GT_Dict const& dict,
                vecbvec const& reads, VecPQVec const& quals,
                ReadPathVec& paths )
    : mHBV(hbv), mDict(dict), mReads(reads), mQuals(quals), mPaths(paths) {}

    void operator()( size_t edgeId )
    { bvec const& edge = mHBV.EdgeObject(edgeId);
      if ( edge.size() < KLEN ) return;
      GT_Kmer kmer(edge.begin());
      int eOffset = 0;
      mLocs.clear();
      KmerLocs const* pLocs;
      if ( (pLocs = mDict.lookup(kmer)) )
        addLocs(*pLocs,eOffset);
      for ( auto itr=edge.begin(KLEN),end=edge.end(); itr != end; ++itr )
      { kmer.toSuccessor(*itr); eOffset += 1;
        if ( (pLocs = mDict.lookup(kmer)) )
          addLocs(*pLocs,eOffset); }
      std::sort(mLocs.begin(),mLocs.end());
      mLocs.erase(std::unique(mLocs.begin(),mLocs.end()),mLocs.end());
      for ( auto const& loc : mLocs )
      { ReadPath& path = mPaths[loc.getReadId()];
        if ( !path.empty() && path[0] == NOT_AN_EDGE )
          continue; // read placement already known not to be unique
        if ( isGood(edgeId,loc) )
        { static SpinLockedData gLock;
          SpinLocker locker(gLock);
          if ( !path.empty() ) path[0] = NOT_AN_EDGE; // mark as not unique
          else
          { path.push_back(edgeId);
            path.setOffset(-loc.getOffset()); } } }
      mLocs.clear(); }

    void cleanAmbiguousPlacements( vec<size_t> const& readIds )
    { for ( size_t readId : readIds )
      { ReadPath& path = mPaths[readId];
        if ( !path.empty() && path[0] == EdgeProc::NOT_AN_EDGE )
        { path.clear(); path.setOffset(0); } } }

private:
    void addLocs( KmerLocs const& locs, int eOffset )
    { auto lItr=locs.locsBegin(), lEnd=locs.locsEnd();
      auto beg=mLocs.begin(),end=mLocs.end();
      while ( beg != end && lItr != lEnd )
      { --lEnd; Loc loc(lEnd->getReadId(),lEnd->getOffset()-eOffset);
        if ( !(loc == *--end) ) { ++lEnd; break; } }
      for ( ; lItr != lEnd; ++lItr )
        mLocs.push_back(Loc(lItr->getReadId(),lItr->getOffset()-eOffset)); }

    // a loc is good if there are no high-quality mismatches, and there's at
    // least one window with few enough mismatches
    bool isGood( size_t edgeId, Loc const& loc )
    { size_t readId = loc.getReadId();
      bvec const& read = mReads[readId];
      mQuals[readId].unpack(&mQVec);
      bvec const& edge = mHBV.EdgeObject(edgeId);
      int offset = -loc.getOffset();
      auto rBeg = read.begin(), rEnd = read.end();
      auto eBeg = edge.begin(), eEnd = edge.end();
      auto qItr = mQVec.begin();
      if ( offset >= 0 ) eBeg += offset;
      else { rBeg -= offset; qItr -= offset; }
      if ( eEnd-eBeg < WINDOW || rEnd-rBeg < WINDOW ) return false;
      auto rItr = rBeg, eItr = eBeg;
      int misMatches = 0;
      for ( auto rWin=rBeg+WINDOW; rItr != rWin; ++rItr,++eItr,++qItr )
        if ( *rItr != *eItr )
        { if ( *qItr >= TRUSTED_QUAL ) return false;
          misMatches += 1; }
      bool good = misMatches <= MAX_MISMATCHES;
      for ( ; rItr != rEnd && eItr != eEnd; ++rItr,++eItr,++qItr,++rBeg,++eBeg )
      { if ( *rItr != *eItr )
        { if ( *qItr >= TRUSTED_QUAL ) return false;
          misMatches += 1; }
        if ( *rBeg != *eBeg ) misMatches -= 1;
        if ( misMatches <= MAX_MISMATCHES ) good = true; }
      return good; }

    HyperBasevector const& mHBV;
    GT_Dict const& mDict;
    vecbvec const& mReads;
    VecPQVec const& mQuals;
    ReadPathVec& mPaths;
    vec<Loc> mLocs;
    QualVec mQVec;
};

}

void PartnersToEnds(const HyperBasevector &hbv, ReadPathVec &paths,
                    const vecbasevector &reads, const VecPQVec &quals) {
     double clock = WallClockTime();

     // find unplaced partners of reads near sinks

     vec<size_t> readIds;
     OutputLog(2) << "finding interesting reads" << std::endl;

     size_t nKmers = findInterestingReadIds(hbv, paths, reads, &readIds);
     size_t nReads = readIds.size();
     if (nReads == 0) return;

     // kmerize those reads, and reduce them into a dictionary of KmerLocs

     OutputLog(2) << "building dictionary" << std::endl;
     GT_Dict *pDict = new GT_Dict(nKmers);

     Mempool locsAlloc;
     size_t const MAX_MULTIPLICITY = 80;
     OutputLog(2) << "reducing" << std::endl;
     {
          RMRE rmre(MREReadProc(readIds, reads, locsAlloc, MAX_MULTIPLICITY, pDict));
          rmre.run(nKmers, 0ul, nReads, RMRE::VERBOSITY::QUIET);
     }

     // kmerize edges, setting the multiplicity for existing dictionary entries

     OutputLog(2) << "kmerizing" << std::endl;
     {
          EMRE emre(MREEdgeProc(hbv, pDict));
          auto const &edges = hbv.Edges();
          size_t edgeKmers = kmerCount(edges.begin(), edges.end(), KLEN);
          emre.run(edgeKmers, 0ul, size_t(hbv.E()), EMRE::VERBOSITY::QUIET);
     }

     // remove dictionary entries having too great a kmer multiplicity

     OutputLog(2) << "cleaning" << std::endl;
     pDict->remove_if([](KmerLocs const &kLocs) { return kLocs.getTotalLocs() > MAX_MULTIPLICITY; });

     // find a uniquely aligning edge and path the read on that edge

     OutputLog(2) << "finding uniquely aligning edges" << std::endl;
     EdgeProc proc(hbv, *pDict, reads, quals, paths);
     #pragma omp parallel
     {
        auto p=proc;
        #pragma omp for
        for (auto i = 0; i < hbv.E(); ++i) {
            p(i);
        }
     }
     proc.cleanAmbiguousPlacements(readIds);
     delete pDict;
}


double CNIntegerFraction(const HyperBasevector& hb, const vec<vec<covcount>>& covs,
			 const double frac, const int min_edge_size) {
    
    size_t edge_count = 0;  // all edges larger than min_edge_size
    size_t good_count = 0;  // edges within frac of an int

    for ( size_t ss = 0; ss < covs.size(); ss++ ) // for each sample
	for ( int e = 0; e < hb.E(); e++ ) {
	    int len = hb.EdgeLength(e);
	    if (len >= min_edge_size && covs[ss][e].Def( ) ) {
		edge_count++;
		double frac_cov = covs[ss][e].Cov();
		int int_cov = frac_cov + 0.5;
		if ( abs(int_cov - frac_cov) <= frac)
		    good_count++;
	    }
	}
    return static_cast<double>(good_count) / edge_count;
}
