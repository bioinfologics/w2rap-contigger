/*
 * SupportedHyperBasevector6.cc
 *
 *  Created on: Apr 3, 2013
 *      Author: tsharpe
 */

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "Basevector.h"
#include "CoreTools.h"
#include "IteratorRange.h"
#include "ParallelVecUtilities.h"
#include "Qualvector.h"
#include "paths/long/CreateGenome.h"
//#include "paths/long/EvalByReads.h"
#include "paths/long/KmerCount.h"
#include "paths/long/SupportedHyperBasevector.h"
#include "random/Bernoulli.h"
#include "system/WorklistN.h"
#include <algorithm>
#include <set>

namespace { // open anonymous namespace

class BubbleProc
{
public:
    BubbleProc( SupportedHyperBasevector const& shbv,
                SupportedHyperBasevector::BubbleParams const& parms,
                SupportedHyperBasevector::BubbleAux const& aux,
                long_logging_control const& log_control,
                long_logging const& logc,
                vec< std::pair<int,vec<int> > >* pSubs,
                const long_heuristics& heur )
    : mSHBV(shbv), mParms(parms), mAux(aux), mLogControl(log_control),
      mLogC(logc), mpSubs(pSubs), mheur(&heur) {}

    void operator()( size_t edgeId )
    { vec<int> substPath;
      if ( mSHBV.IsPoppable(edgeId,mParms,mAux,mLogControl,mLogC,&substPath,*mheur) )
      { SpinLocker lock(gLockedData); mpSubs->push(edgeId,substPath); } }

private:
    SupportedHyperBasevector const& mSHBV;
    SupportedHyperBasevector::BubbleParams const& mParms;
    SupportedHyperBasevector::BubbleAux const& mAux;
    long_logging_control const& mLogControl;
    long_logging const& mLogC;
    vec< std::pair<int,vec<int> > >* mpSubs;
    static SpinLockedData gLockedData;
    const long_heuristics* mheur;
};

SpinLockedData BubbleProc::gLockedData;

} // close anonymous namespace

void TruncateMe( vec<int>& p, Bool& del, int& start, int& stop, 
     const vec<Bool>& used, const SupportedHyperBasevector& shb )
{    
     del = False;
     vec< vec<int> > subs, subs_origin;
     vec<int> s, s_origin;
     for ( int j = 0; j <= p.isize( ); j++ )
     {    if ( j == p.isize( ) || ( p[j] >= 0 && !used[ p[j] ] ) )
          {    if ( s.nonempty( ) ) 
               {    if ( s.back( ) < 0 ) 
                    {    s.pop_back( );
                         s_origin.pop_back( );    }
                    subs.push_back(s);
                    subs_origin.push_back(s_origin);    }
               s.clear( );
               s_origin.clear( );    }
          else 
          {    s.push_back( p[j] );
               s_origin.push_back(j);    }    }
     vec<int> nkmers( subs.size( ), 0 ), ids( subs.size( ), vec<int>::IDENTITY );
     for ( int j = 0; j < subs.isize( ); j++ )
     {    for ( int l = 0; l < subs[j].isize( ); l++ )
               nkmers[j] += shb.EdgeLengthKmers( subs[j][l] );    }
     ReverseSortSync( nkmers, ids );
     if ( nkmers.empty( ) || ( nkmers.size( ) >= 2 && nkmers[0] == nkmers[1] ) )
          del = True;
     else
     {    p = subs[ ids[0] ];
          start = subs_origin[ ids[0] ].front( );
          stop = subs_origin[ ids[0] ].back( );    }    }

void SupportedHyperBasevector::TruncatePaths( const long_logging& logc )
{    double clock = WallClockTime( );
     vec<Bool> used, to_delete( NPaths( ), False );
     Used(used);
     for ( int i = 0; i < NPaths( ); i++ )
     {    Bool del;
          int start, stop;
          TruncateMe( PathMutable(i), del, start, stop, used, *this );
          if (del) to_delete[i] = True;    }
     EraseIf( PathsMutable( ), to_delete );
     EraseIf( WeightsFwMutable( ), to_delete );
     EraseIf( WeightsRcMutable( ), to_delete );
     // Temporary if.
     if ( WeightsFwOrigin( ).size( ) == WeightsFw( ).size( ) )
     {    EraseIf( WeightsFwOriginMutable( ), to_delete );
          EraseIf( WeightsRcOriginMutable( ), to_delete );    }
     to_delete.resize_and_set( NPairs( ), False );
     for ( int i = 0; i < NPairs( ); i++ )
     {    Bool del1, del2;
          int start1, stop1, start2, stop2;
          vec<int> left = PairLeftMutable(i), right = PairRightMutable(i);
          TruncateMe( PairLeftMutable(i), del1, start1, stop1, used, *this );
          TruncateMe( PairRightMutable(i), del2, start2, stop2, used, *this );
          if ( del1 || del2 ) to_delete[i] = True;
          else
          {    int sub = 0;
               for ( int k = stop1; k < left.isize( ); k++ )
                    sub += EdgeLengthKmers( left[k] );
               for ( int k = 0; k < start2; k++ )
                    sub += EdgeLengthKmers( right[k] );
               AddTrim( i, -sub );    }    }
     EraseIf( PairsMutable( ), to_delete );
     EraseIf( PairDataMutable( ), to_delete );
     UniqueOrderPaths( );
     REPORT_TIME( clock, "used in TruncatePaths" );    }


bool SupportedHyperBasevector::IsPoppable( int e, BubbleParams const& parms,
                                        BubbleAux const& aux,
                                        long_logging_control const& log_control,
                                        long_logging const& logc,
                                        vec<int>* pSubstPath,
     const long_heuristics& heur ) const
{
     int v = aux.to_left[e];
     int w = aux.to_right[e];

     // don't consider certain edges with strong evidence
     if ( aux.mult_fw[e] > parms.max_pop_del &&
             aux.mult_rc[e] > parms.max_pop_del ) return false;
     if ( Max( aux.mult_fw[e], aux.mult_rc[e] ) > parms.max_pop_del2 ) return false;

     // establish a min and max path length for alternate paths
     int ne = EdgeObject(e).isize( ) - K( ) + 1;
     int low = ne - int( floor( parms.min_pop_ratio * parms.delta_kmers ) );
     int high = ne + int( floor( parms.min_pop_ratio * parms.delta_kmers ) );

     vec< vec<int> > paths, paths_final;
     vec<int> pathlens, pathlens_final;
     vec<fix64_6> mult_fw_final, mult_rc_final;

     // beginning of depth-first recursion from the starting vertex
     // to look for alternate paths
     // NOTE: it would probably be a little faster to go breadth-first
     // instead.  you might more quickly discover multiple alternatives
     // without fully exploring tiny cycles.
     for ( int l = 0; l < From(v).isize( ); l++ )
     {    int f = EdgeObjectIndexByIndexFrom( v, l );
          if ( f != e )
          { paths.push_back(vec<int>(1,f));
            pathlens.push_back( EdgeObject(f).isize( ) - K( ) + 1 );   } }

     // while there are alternatives to consider
     while( paths.nonempty( ) )
     {    vec<int> p = paths.back( );
          int n = pathlens.back( );
          paths.pop_back( );
          // std::cout << "found path " << printSeq(p) << std::endl; // XXXXXXXXXXXXXXXXX
          pathlens.pop_back( );
          int x = aux.to_right[ p.back( ) ];

          // if the path under consideration has the right length, and
          // ends up at the same vertex as the subject edge
          if ( n >= low && n <= high && x == w )
          {    fix64_6 pmult_fw = 0, pmult_rc = 0;
               for ( int l = 0; l < aux.paths_index[ p[0] ].isize( ); l++ )
               {    int id = aux.paths_index[ p[0] ][l].first;
                    int pos = aux.paths_index[ p[0] ][l].second;
                    if ( Path(id).Contains( p, pos ) )
                    {    pmult_fw += WeightFw(id);
                         pmult_rc += WeightRc(id);    }    }
               // if it has appropriate weight, it's a possibility
               if ( pmult_fw <= parms.max_pop_del || pmult_rc <= parms.max_pop_del )
                   continue; // if weight insufficient, no need to extend further
               paths_final.push_back(p);
               // std::cout << "accepting path " << printSeq(p) << std::endl; // XXX
               pathlens_final.push_back(n);
               mult_fw_final.push_back(pmult_fw);
               mult_rc_final.push_back(pmult_rc);
               // we're only looking for situations where there's a
               // single valid alternative path, so quit early when
               // we've discovered two alternatives
               if ( paths_final.size() > 1 ) break;    }

          // if the path is already too long, don't extend it further
          if ( n > high ) continue;

          // extend path and continue recursive exploration
          for ( int l = 0; l < From(x).isize( ); l++ )
          {    int f = EdgeObjectIndexByIndexFrom( x, l );
               if ( f != e )
               { paths.push_back(p);
                 paths.back().push_back(f);
                 pathlens.push_back(
                    n + EdgeObject(f).isize( ) - K( ) + 1 );    }   }   }

     if ( !paths_final.solo( ) ) return false;
     if ( Abs( pathlens_final[0] - ne ) > parms.delta_kmers ) return false;
     /*
     if ( ( mult_fw_final[0] < parms.min_pop_ratio * aux.mult_fw[e] )
          && ( mult_rc_final[0] < parms.min_pop_ratio * aux.mult_rc[e] ) )
     {    continue;    }
     */

     // Proper statistical test - note missing from other documentation.
     // Compare similar code in ImproveLongHyper.cc.

     const double max_asym_rarity = 0.00001;
     double f1 = mult_fw_final[0].ToDouble(), r1 = mult_rc_final[0].ToDouble();;
     double f2 = aux.mult_fw[e].ToDouble(), r2 = aux.mult_rc[e].ToDouble();
     if ( f2 > r2 || ( f2 == r2 && f1 > r1 ) )
     {    std::swap( f1, r1 );
          std::swap( f2, r2 );    }
     int n = int(floor(f1+r1+f2+r2));
     if ( f1 == 0 || n == 0 || n > 10000 ) return false;

     long double p;
     double q;
     if ( !heur.FIX_ASYMMETRY_BUG )
     {    p = Min( 0.5, f1/(f1+r1) ) / 2;
          q = BinomialSum( n, int(ceil(f2)), p );    }
     else
     {    p = 0.5;
          q = BinomialSum( n, int(ceil( Min(f1+r1,f2+r2) )), p );    }
     if ( q >= max_asym_rarity ) return false;

     // Report result.

     if ( logc.verb[ "POP_BUBBLES" ] >= 1 )
     {    static SpinLockedData gLockedData;
          SpinLocker lock(gLockedData);
          std::cout << "\nPopBubbles: replace e = " << e << " by "
               << printSeq( paths_final[0] ) << std::endl;
          PRINT2( aux.mult_fw[e], aux.mult_rc[e] );
          PRINT2( mult_fw_final[0], mult_rc_final[0] );    }
     *pSubstPath = paths_final[0];
     return true; }

namespace {
    typedef int EdgeId;
    typedef int VertexId;
    typedef unsigned SegId;

// An ordered sequence of EdgeIds that describe, e.g., a traversal of a graph
    typedef vec<EdgeId> SHBVPath;

// A collection of edgeIds describing, e.g., the set of edges departing a vertex
    typedef vec<EdgeId> EdgeCollection;

// A pair of edgeIds
    typedef std::pair<EdgeId, EdgeId> EdgePair;

// An element of some read's path on the graph
    struct ReadSegment {
        ReadSegment() = default;

        ReadSegment(EdgeId readId, SegId segId)
                : mReadId(readId), mSegId(segId) {}

        EdgeId mReadId;
        SegId mSegId;
    };

// an operation we think is a valid simplification of the graph
    typedef triple<EdgeCollection, EdgeCollection, EdgeCollection> SHBV_Join;

// instructions about what pull-aparts to do
    struct Recommendations : private SpinLockedData {
        void addRecommendation(std::set<SHBVPath> const &bridges11,
                               std::set<SHBVPath> const &bridges22,
                               EdgeCollection const &reach);

        vec<SHBV_Join> mJoins;
        vec<vec<SHBVPath>> mPaths1, mPaths2;
    };

// adding a recommendation is thread safe
    void Recommendations::addRecommendation(std::set<SHBVPath> const &bridges11,
                                            std::set<SHBVPath> const &bridges22,
                                            EdgeCollection const &reach) {
        ForceAssert(bridges11.size());
        ForceAssert(bridges22.size());
        EdgeId e1 = bridges11.begin()->front(), e2 = bridges22.begin()->front();
        EdgeId f1 = bridges11.begin()->back(), f2 = bridges22.begin()->back();
        EdgeCollection eee{e1, e2};
        EdgeCollection fff{f1, f2};
        SHBV_Join join(eee, reach, fff);
        SpinLocker lock(*this);
        mJoins.push_back(join);
        mPaths1.push(bridges11.begin(), bridges11.end());
        mPaths2.push(bridges22.begin(), bridges22.end());
    }

// paths and cum weights for all paths from some edge e to edges f1 or f2
    struct PathWeights {
        PathWeights(SupportedHyperBasevector const &hbv,
                    vec<ReadSegment> const &segs,
                    EdgeId e, EdgeId f1, EdgeId f2);

        // get the total length of a portion of a path in kmers
        static int getSpan(SupportedHyperBasevector const &hbv,
                           SHBVPath::const_iterator itr,
                           SHBVPath::const_iterator end) {
            int span = 0;
            while (itr != end) {
                span += hbv.EdgeLengthKmers(*itr);
                ++itr;
            }
            return span;
        }

        fix64_6 weight1, weight2; // cum weights of paths to f1 and f2
        std::set<SHBVPath> bridges1, bridges2; // distinct paths to f1 and f2
        int mSpan; // max distance in kmers on any path from the edge to edge f1 or f2
    };

    PathWeights::PathWeights(SupportedHyperBasevector const &hbv,
                             vec<ReadSegment> const &segs,
                             EdgeId e, EdgeId f1, EdgeId f2)
            : weight1(0), weight2(0), mSpan(0) {
        auto segsEnd = segs.end();
        for (auto segsItr = segs.begin(); segsItr != segsEnd; ++segsItr) {
            SHBVPath const &path = hbv.Path(segsItr->mReadId);
            fix64_6 const &weight = hbv.Weight(segsItr->mReadId);
            auto from = path.begin() + segsItr->mSegId;
            auto beg = from + 1;
            auto end = path.end();
            auto itr = std::find(beg, end, f1);
            if (itr != end) {
                if (std::find(beg, itr, e) != itr) // skip it if there's a later
                    continue;                      // instance of e
                weight1 += weight;
                bridges1.insert(SHBVPath(from, itr + 1));
            } else if ((itr = std::find(beg, end, f2)) != end) {
                if (std::find(beg, itr, e) != itr) // skip if there's a later
                    continue;                      // instance of e
                weight2 += weight;
                bridges2.insert(SHBVPath(from, itr + 1));
            } else
                continue;
            mSpan = std::max(mSpan, getSpan(hbv, beg, itr));
        }
    }

// a thing that examines a vertex to decide whether it can be pulled apart
    class JoinDiscoverer {
    public:
        JoinDiscoverer(SupportedHyperBasevector const &hbv,
                       vec<int> const &toLeft, vec<int> const &toRight,
                       int verbosity, double minWeightSplit)
                : mHBV(hbv), mToLeft(toLeft), mToRight(toRight),
                  mPathsIndex(hbv.EdgeObjectCount()),
                  mReadLenFudge(hbv.MedianCorrectedReadLength()),
                  mVerbosity(verbosity),
                  mMinWeightSplit(minWeightSplit), mMinWeightSplitLow(2) {
            EdgeId nEdges = hbv.Paths().isize();
            for (EdgeId edgeId = 0; edgeId != nEdges; edgeId++) {
                vec<EdgeId> const &path = hbv.Path(edgeId);
                SegId nSegs = path.size();
                for (SegId segId = 0; segId != nSegs; segId++)
                    mPathsIndex[path[segId]].push_back(ReadSegment(edgeId, segId));
            }
        }

        void processVertex(int vertexId, Recommendations *) const;

    private:
        static size_t const MAX_AFFIX_LEN = 5;

        // scan read paths to find unique prefixes of paths that enter a vertex
        vec<SHBVPath> findPrefixes(EdgeCollection const &leftAlts) const {
            vec<SHBVPath> prefixes;
            auto altsEnd = leftAlts.end();
            for (auto altsItr = leftAlts.begin(); altsItr != altsEnd; ++altsItr) {
                vec<ReadSegment> const &segs = mPathsIndex[*altsItr];
                auto segsEnd = segs.end();
                for (auto segsItr = segs.begin(); segsItr != segsEnd; ++segsItr) {
                    SHBVPath const &path = mHBV.Path(segsItr->mReadId);
                    size_t nextEle = segsItr->mSegId + 1;
                    auto pathEnd = path.begin() + nextEle;
                    auto pathBeg = pathEnd - std::min(nextEle, MAX_AFFIX_LEN);
                    bool found = false;
                    auto end = prefixes.end();
                    for (auto itr = prefixes.begin(); itr != end; ++itr) {
                        SHBVPath const &prefix = *itr;
                        size_t len = std::min(nextEle, prefix.size());
                        if (std::equal(pathEnd - len, pathEnd, prefix.end() - len)) {
                            if (len < MAX_AFFIX_LEN && nextEle > len)
                                *itr = SHBVPath(pathBeg, pathEnd);
                            found = true;
                            break;
                        }
                    }
                    if (!found)
                        prefixes.push_back(SHBVPath(pathBeg, pathEnd));
                }
            }
            return prefixes;
        }

        // scan read paths to find unique suffixes of paths that leave a vertex
        vec<SHBVPath> findSuffixes(EdgeCollection const &rightAlts) const {
            vec<SHBVPath> suffixes;
            auto altsEnd = rightAlts.end();
            for (auto altsItr = rightAlts.begin(); altsItr != altsEnd; ++altsItr) {
                vec<ReadSegment> const &segs = mPathsIndex[*altsItr];
                auto segsEnd = segs.end();
                for (auto segsItr = segs.begin(); segsItr != segsEnd; ++segsItr) {
                    SHBVPath const &path = mHBV.Path(segsItr->mReadId);
                    auto pathBeg = path.begin() + segsItr->mSegId;
                    size_t pathLen = std::min(MAX_AFFIX_LEN,
                                              path.size() - segsItr->mSegId);
                    bool found = false;
                    auto end = suffixes.end();
                    for (auto itr = suffixes.begin(); itr != end; ++itr) {
                        SHBVPath const &suffix = *itr;
                        size_t len = std::min(pathLen, suffix.size());
                        if (std::equal(pathBeg, pathBeg + len, suffix.begin())) {
                            if (pathLen > len)
                                *itr = SHBVPath(pathBeg, pathBeg + pathLen);
                            found = true;
                            break;
                        }
                    }
                    if (!found)
                        suffixes.push_back(SHBVPath(pathBeg, pathBeg + pathLen));
                }
            }
            return suffixes;
        }

        // proper edge pairs tile the set of path prefixes or suffixes, i.e.,
        // every path contains one or the other edge, and no path contains both
        static vec<EdgePair> getProperEdgePairs(vec<SHBVPath> const &paths) {
            std::map<EdgeId, BitVec> map;
            size_t nPaths = paths.size();
            for (size_t idx = 0; idx != nPaths; ++idx) {
                SHBVPath const &path = paths[idx];
                for (auto itr = path.begin(), end = path.end(); itr != end; ++itr) {
                    BitVec &bits = map[*itr];
                    bits.resize(nPaths, false);
                    bits.set(idx, true);
                }
            }
            vec<EdgePair> pairs;
            for (auto itr = map.begin(), end = map.end(); itr != end; ++itr) {
                auto itr2 = itr;
                while (++itr2 != end) {
                    if ((itr->second & itr2->second).isUniformlyFalse() &&
                        nor(itr->second, itr2->second).isUniformlyFalse())
                        pairs.push_back(EdgePair(itr->first, itr2->first));
                }
            }
            return pairs;
        }

        // gets a sorted list of EdgeIds that lie between e1 or e2 and f1 or f2.
        void getReach(PathWeights const &pw1, PathWeights const &pw2,
                      EdgeCollection *pReach) const {
            pReach->clear();
            for (SHBVPath const &path : pw1.bridges1) {
                ForceAssertGe(path.size(), 2u);
                pReach->insert(pReach->end(), path.begin() + 1, path.end() - 1);
            }
            for (SHBVPath const &path : pw1.bridges2) {
                ForceAssertGe(path.size(), 2u);
                pReach->insert(pReach->end(), path.begin() + 1, path.end() - 1);
            }
            for (SHBVPath const &path : pw2.bridges1) {
                ForceAssertGe(path.size(), 2u);
                pReach->insert(pReach->end(), path.begin() + 1, path.end() - 1);
            }
            for (SHBVPath const &path : pw2.bridges1) {
                ForceAssertGe(path.size(), 2u);
                pReach->insert(pReach->end(), path.begin() + 1, path.end() - 1);
            }

            UniqueSort(*pReach);
        }

        // returns false if the reach is not a complete sub-graph delimited by
        // e1, e2, f1, and f2 (i.e., if non-reach, non-boundary edges enter or exit
        // from reach edges).
        bool checkReach(EdgeCollection const &reach,
                        EdgeId e1, EdgeId e2, EdgeId f1, EdgeId f2) const {
            vec<int> reachVertices;
            reachVertices.reserve(reach.size() + 2);
            for (EdgeId e : reach) {
                reachVertices.push_back(mToLeft[e]);
            }
            reachVertices.push_back(mToLeft[f1]);
            reachVertices.push_back(mToLeft[f2]);
            UniqueSort(reachVertices);

            auto beg = reach.begin(), end = reach.end();
            for (VertexId v : reachVertices) {
                for (EdgeId e : mHBV.ToEdgeObj(v))
                    if (e != e1 && e != e2 && !std::binary_search(beg, end, e))
                        return false;
                for (EdgeId e : mHBV.FromEdgeObj(v))
                    if (e != f1 && e != f2 && !std::binary_search(beg, end, e))
                        return false;
            }
            return true;
        }

        // do the weights allow a pull apart?
        bool pullApart(fix64_6 const &weight_11, fix64_6 const &weight_12,
                       fix64_6 const &weight_21, fix64_6 const &weight_22,
                       int span) const {
            bool result = false;
            if (weight_11 >= mMinWeightSplitLow && weight_22 >= mMinWeightSplitLow &&
                weight_12 == 0 && weight_21 == 0)
                result = true;
            if (weight_11 >= mMinWeightSplit && weight_22 >= mMinWeightSplit &&
                weight_12 + weight_21 <= 2 && span <= mReadLenFudge)
                result = true;
            if (weight_11 >= mMinWeightSplit / 2 && weight_22 >= mMinWeightSplit / 2 &&
                weight_12 + weight_21 < 2 && span <= mReadLenFudge)
                result = true;
            return result;
        }

        // thread-safe logging of a join
        void log(EdgeId e1, EdgeId e2, EdgeId f1, EdgeId f2,
                 EdgeCollection const &reach,
                 std::set<SHBVPath> const &bridge11,
                 std::set<SHBVPath> const &bridge22,
                 fix64_6 const &weight_11, fix64_6 const &weight_12,
                 fix64_6 const &weight_21, fix64_6 const &weight_22,
                 int span, bool result) const {
            SpinLocker lock(gLockCOUT);
            std::cout << "e1=" << e1
                      << " e2=" << e2
                      << " f1=" << f1
                      << " f2=" << f2
                      << "\nreach: " << rangePrinter(reach.begin(), reach.end(), ",")
                      << "\nbridge11:";
            int idx = 0;
            for (auto itr = bridge11.begin(), end = bridge11.end(); itr != end; ++itr)
                std::cout << "\n[" << ++idx << "] " << rangePrinter(itr->begin(), itr->end(), ",");
            std::cout << "\nbridge22:";
            idx = 0;
            for (auto itr = bridge22.begin(), end = bridge22.end(); itr != end; ++itr)
                std::cout << "\n[" << ++idx << "] " << rangePrinter(itr->begin(), itr->end(), ",");
            std::cout << "\nw11=" << weight_11
                      << " w22=" << weight_22
                      << " w12=" << weight_12
                      << " w21=" << weight_21
                      << " span=" << span << '\n';
            if (!result) std::cout << "rejected\n";
            std::cout << std::endl;
        }

        SupportedHyperBasevector const &mHBV;
        vec<int> const &mToLeft;
        vec<int> const &mToRight;
        vec<vec<ReadSegment>> mPathsIndex;
        int mReadLenFudge;
        int mVerbosity;
        double mMinWeightSplit;
        fix64_6 mMinWeightSplitLow;
        static SpinLockedData gLockCOUT; // logging lock
    };

    SpinLockedData JoinDiscoverer::gLockCOUT;
    size_t const JoinDiscoverer::MAX_AFFIX_LEN;

// can this vertex be pulled apart?
    void JoinDiscoverer::processVertex(int vertexId, Recommendations *pRecs) const {
        EdgeCollection const &leftAlts = mHBV.ToEdgeObj(vertexId);
        size_t nLeftAlts = leftAlts.size();
        // only process nodes with some left-diversity to cut down on re-processing
        if (nLeftAlts <= 1)
            return;

        EdgeCollection const &rightAlts = mHBV.FromEdgeObj(vertexId);
        size_t nRightAlts = rightAlts.size();
        if (!nRightAlts) // if its a dead-end, skip it
            return;

        // get all distinct paths that lead to this vertex
        vec<SHBVPath> prefixes = findPrefixes(leftAlts);
        size_t nPrefixes = prefixes.size();
        if (nPrefixes <= 1)
            return;

        // get all distinct paths that lead from this vertex
        vec<SHBVPath> suffixes = findSuffixes(rightAlts);
        size_t nSuffixes = suffixes.size();
        if (nSuffixes <= 1)
            return;

        // get candidate pairs of "e" edges upstream of vertex
        vec<EdgePair> ePairs = getProperEdgePairs(prefixes);
        if (ePairs.empty())
            return;

        // get candidate pairs of "f" edges downstream of vertex
        vec<EdgePair> fPairs = getProperEdgePairs(suffixes);
        if (fPairs.empty())
            return;

        // for each of the valid pairs of e's and f's, find all the reads that
        // go from one of the e's to one of the f's, and accumulate the read weights
        EdgeCollection reach;
        for (auto eItr = ePairs.begin(), eEnd = ePairs.end(); eItr != eEnd; ++eItr) {
            EdgeId e1 = eItr->first, e2 = eItr->second;
            for (auto fItr = fPairs.begin(), fEnd = fPairs.end(); fItr != fEnd; ++fItr) {
                EdgeId f1 = fItr->first, f2 = fItr->second;
                if (e1 == f1 || e1 == f2 || e2 == f1 || e2 == f2)
                    continue;

                PathWeights pw1(mHBV, mPathsIndex[e1], e1, f1, f2);
                PathWeights pw2(mHBV, mPathsIndex[e2], e2, f1, f2);

                getReach(pw1, pw2, &reach);
                // make sure {e1,e2,f1,f2} bound a complete sub-graph
                if (!checkReach(reach, e1, e2, f1, f2))
                    continue;

                int span = std::max(pw1.mSpan, pw2.mSpan);
                if (pullApart(pw1.weight1, pw1.weight2,
                              pw2.weight1, pw2.weight2, span)) {
                    pRecs->addRecommendation(pw1.bridges1, pw2.bridges2, reach);
                    if (mVerbosity >= 2)
                        log(e1, e2, f1, f2, reach, pw1.bridges1, pw2.bridges2,
                            pw1.weight1, pw1.weight2, pw2.weight1, pw2.weight2,
                            span, true);
                } else if (pullApart(pw1.weight2, pw1.weight1,
                                     pw2.weight2, pw2.weight1, span)) {
                    pRecs->addRecommendation(pw1.bridges2, pw2.bridges1, reach);
                    if (mVerbosity >= 2)
                        log(e1, e2, f2, f1, reach, pw1.bridges2, pw2.bridges1,
                            pw1.weight2, pw1.weight1, pw2.weight2, pw2.weight1,
                            span, true);
                } else if (mVerbosity >= 3 &&
                           !pw1.bridges1.empty() && !pw2.bridges2.empty()) {
                    log(e1, e2, f1, f2, reach, pw1.bridges1, pw2.bridges2,
                        pw1.weight1, pw1.weight2, pw2.weight1, pw2.weight2,
                        span, false);
                } else if (mVerbosity >= 3 &&
                           !pw1.bridges2.empty() && !pw2.bridges1.empty()) {
                    log(e1, e2, f2, f1, reach, pw1.bridges2, pw2.bridges1,
                        pw1.weight2, pw1.weight1, pw2.weight2, pw2.weight1,
                        span, false);
                }
            }
        }
    }


    void SupportedHyperBasevector::DeleteLowCoverage(const long_heuristics &heur,
                                                     const long_logging_control &log_control,
                                                     const long_logging &logc) {
        double mclock = WallClockTime();
        vec<int> to_left, to_right, to_delete;
        ToLeft(to_left), ToRight(to_right);
        const double low_cov = 2.0;
        vec<fix64_6> cov(EdgeObjectCount(), 0.0);
        for (int i = 0; i < NPaths(); i++) {
            for (int j = 0; j < Path(i).isize(); j++)
                cov[Path(i)[j]] += Weight(i);
        }

        const double pceil = 0.99;
        const double minp = 0.0001;
        vec<vec<fix64_6> > weights(EdgeObjectCount());

        vec<triple<int64_t, int, fix64_6> > wid;
        vec<double> p(EdgeObjectCount(), 1.0);
        if (heur.NEW_LC_FILT) {
            for (int i = 0; i < NPaths(); i++) {
                for (int k = 0; k < WeightFwOrigin(i).isize(); k++) {
                    wid.push(WeightFwOrigin(i)[k].second, i,
                             WeightFwOrigin(i)[k].first);
                }
                for (int k = 0; k < WeightRcOrigin(i).isize(); k++) {
                    wid.push(WeightRcOrigin(i)[k].second, i,
                             WeightRcOrigin(i)[k].first);
                }
            }
            ParallelSort(wid);
            for (int64_t i = 0; i < wid.jsize(); i++) {
                int64_t j;
                for (j = i + 1; j < wid.jsize(); j++)
                    if (wid[j].first != wid[i].first) break;
                vec<std::pair<int, fix64_6> > w;
                for (int64_t k = i; k < j; k++) {
                    vec<int> p = Path(wid[k].second);
                    UniqueSort(p);
                    for (int l = 0; l < p.isize(); l++)
                        w.push(p[l], wid[k].third);
                }
                Sort(w);
                for (int r = 0; r < w.isize(); r++) {
                    int s;
                    for (s = r + 1; s < w.isize(); s++)
                        if (w[s].first != w[r].first) break;
                    fix64_6 x = 0;
                    for (int t = r; t < s; t++)
                        x += w[t].second;
                    weights[w[r].first].push_back(x);
                    r = s - 1;
                }
                i = j - 1;
            }
            for (int e = 0; e < EdgeObjectCount(); e++) {
                Sort(weights[e]);
                for (int j = 0; j < weights[e].isize(); j++) {
                    if (p[e] < minp) break;
                    p[e] *= 1.0
                            - Min(pceil, weights[e][j].ToDouble());
                }
            }
        }

        const int min_mult = 5;
        for (int e = 0; e < EdgeObjectCount(); e++) {
            int re = Inv(e), v = to_left[e], w = to_right[e];
            fix64_6 c = cov[e];
            fix64_6 rc = (re < 0 ? 1000000000 : cov[re]);
            fix64_6 alt_c = 0, alt_rc = 0;
            for (int j = 0; j < From(v).isize(); j++)
                alt_c = Max(alt_c, cov[EdgeObjectIndexByIndexFrom(v, j)]);
            for (int j = 0; j < To(w).isize(); j++)
                alt_c = Max(alt_c, cov[EdgeObjectIndexByIndexTo(w, j)]);
            if (re >= 0) {
                int v = to_left[re], w = to_right[re];
                for (int j = 0; j < From(v).isize(); j++) {
                    alt_rc = Max(alt_rc, cov[
                            EdgeObjectIndexByIndexFrom(v, j)]);
                }
                for (int j = 0; j < To(w).isize(); j++) {
                    alt_rc = Max(alt_rc, cov[
                            EdgeObjectIndexByIndexTo(w, j)]);
                }
            }
            if (heur.NEW_LC_FILT) {
                if (((c <= low_cov || p[e] >= minp)
                     && alt_c >= min_mult * c)
                    || ((rc <= low_cov || (re >= 0 && p[re] >= minp))
                        && alt_rc >= min_mult * rc)) {
                    if (logc.verb["LOW_COV"] >= 1) {
                        std::cout << "deleting low-coverage edge " << e
                                  << " (c = " << c << ", alt_c = " << alt_c
                                  << ", p[e] = " << p[e] << ")" << std::endl;
                    }
                    to_delete.push_back(e);
                }
            } else {
                if (heur.LC_CAREFUL && alt_c < low_cov) continue;
                if ((c <= low_cov && alt_c >= min_mult * c)
                    || (rc <= low_cov && alt_rc >= min_mult * rc)) {
                    if (logc.verb["LOW_COV"] >= 1) {
                        vec<int> comp;
                        for (int j = 0; j < From(v).isize(); j++) {
                            int e = EdgeObjectIndexByIndexFrom(v, j);
                            if (cov[e] == alt_c) comp.push_back(e);
                        }
                        for (int j = 0; j < To(w).isize(); j++) {
                            int e = EdgeObjectIndexByIndexTo(w, j);
                            if (cov[e] == alt_c) comp.push_back(e);
                        }
                        UniqueSort(comp);
                        std::cout << "deleting low-coverage edge " << e
                                  << " using competing edge(s) " << printSeq(comp)
                                  << ", c = " << c << ", alt_c = " << alt_c << std::endl;
                    }
                    to_delete.push_back(e);
                }
            }
        }
        DeleteEdges(to_delete);
        vec<Bool> p_to_delete(NPaths(), False);
        for (int i = 0; i < NPaths(); i++) {
            Bool OK = True;
            for (int j = 0; j < Path(i).isize(); j++)
                if (BinMember(to_delete, Path(i)[j])) OK = False;
            if (!OK) p_to_delete[i] = True;
        }
        EraseIf(PathsMutable(), p_to_delete);
        // Temporary if.
        if (WeightsFwOrigin().size() == WeightsFw().size()) {
            EraseIf(WeightsFwOriginMutable(), p_to_delete);
            EraseIf(WeightsRcOriginMutable(), p_to_delete);
        }
        EraseIf(WeightsFwMutable(), p_to_delete);
        EraseIf(WeightsRcMutable(), p_to_delete);
        p_to_delete.resize_and_set(NPairs(), False);
        for (int i = 0; i < NPairs(); i++) {
            Bool OK = True;
            for (int j = 0; j < PairLeft(i).isize(); j++)
                if (BinMember(to_delete, PairLeft(i)[j])) OK = False;
            if (!OK) p_to_delete[i] = True;
            for (int j = 0; j < PairRight(i).isize(); j++)
                if (BinMember(to_delete, PairRight(i)[j])) OK = False;
            if (!OK) p_to_delete[i] = True;
        }
        EraseIf(PairsMutable(), p_to_delete);
        EraseIf(PairDataMutable(), p_to_delete);
        RemoveEdgelessVertices();
        RemoveUnneededVertices();
        REPORT_TIME(mclock, "used deleting low-coverage edges");
        RemoveDeadEdgeObjects();
    };
}