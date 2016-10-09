///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * BuildReadQGraph.cc
 *
 *  Created on: Jan 22, 2014
 *      Author: tsharpe
 */

#include "paths/long/BuildReadQGraph.h"
#include "Basevector.h"
#include "FastaFileset.h"
#include "Intvector.h"
#include "IteratorRange.h"
#include "MapReduceEngine.h"
#include "PairsManager.h"
#include "ParallelVecUtilities.h"
#include "Qualvector.h"
#include "Vec.h"
#include "dna/Bases.h"
#include "feudal/BinaryStream.h"
#include "feudal/VirtualMasterVec.h"
//#include "kmers/BigKPather.h"
#include "kmers/ReadPatherDefs.h"
#include "math/Functions.h"
#include "paths/KmerBaseBroker.h"
#include "paths/UnibaseUtils.h"
#include "paths/long/HBVFromEdges.h"
#include "paths/long/KmerCount.h"
#include "system/SortInPlace.h"
#include "system/SpinLockedData.h"
#include "system/WorklistN.h"
#include <algorithm>
#include <atomic>
#include <fstream>
#include <iostream>
#include <numeric>
#include <utility>
#include <vector>
#include "paths/HyperBasevector.h"
#include "paths/long/ExtendReadPath.h"
#include "paths/long/ShortKmerReadPather.h"

namespace
{
    const unsigned K = 60;
    typedef KMer<K> BRQ_Kmer;
    typedef KMer<K-1> BRQ_SubKmer;
    typedef KmerDictEntry<K> BRQ_Entry;
    typedef KmerDict<K> BRQ_Dict;




    inline void summarizeEntries( BRQ_Entry* e1, BRQ_Entry* e2 )
    {
        KMerContext kc;
        size_t count = 0;
        while ( e2-- != e1 )
        {
            KDef const& kDef = e2->getKDef();
            kc |= kDef.getContext();
            count += std::max(kDef.getCount(),1ul);
        }
        KDef& kDef = e1->getKDef();
        kDef.setContext(kc);
        kDef.setCount(count);
    }

    class EdgeBuilder
    {
    public:
        EdgeBuilder( BRQ_Dict const& dict, vecbvec* pEdges )
                : mDict(dict), mEdges(*pEdges) {}

        void buildEdge( BRQ_Entry const& entry )
        { if ( isPalindrome(entry) )
                make1KmerEdge(entry);
            else if ( upstreamExtensionPossible(entry) )
            { if ( downstreamExtensionPossible(entry) )
                    return;
                extendUpstream(entry); }
            else if ( downstreamExtensionPossible(entry) )
                extendDownstream(entry);
            else
                make1KmerEdge(entry); }

        // not thread-safe
        void simpleCircle( BRQ_Entry const& entry )
        { BRQ_Entry const* pFirstEntry = &entry;
            mEdgeSeq.assign(entry.begin(),entry.end());
            mEdgeEntries.push_back(pFirstEntry);
            KMerContext context = entry.getKDef().getContext();
            BRQ_Kmer kmer(entry);
            while ( true )
            { ForceAssertEq(context.getPredecessorCount(),1u);
                ForceAssertEq(context.getSuccessorCount(),1u);
                unsigned char succCode = context.getSingleSuccessor();
                kmer.toSuccessor(succCode);
                BRQ_Entry const* pEntry = lookup(kmer,&context);
                if ( pEntry == pFirstEntry )
                    break;
                if ( !pEntry->getKDef().isNull() )
                { std::cout << "Failed to close circle.\n";
                    for ( auto beg=mEdgeEntries.begin(),end=mEdgeEntries.end();
                          beg != end; ++beg )
                        std::cout << *beg << ' ' << **beg << '\n';
                    std::cout << pEntry << ' ' << *pEntry << std::endl;
                    CRD::exit(1); }
                mEdgeSeq.push_back(succCode);
                mEdgeEntries.push_back(pEntry); }
            canonicalizeCircle();
            addEdge(); }

    private:
        void canonicalizeCircle()
        { auto itr = std::min_element(mEdgeEntries.begin(),mEdgeEntries.end(),
                                      []( BRQ_Kmer const* pEnt1, BRQ_Kmer const* pEnt2 )
                                      { return *pEnt1 < *pEnt2; });
            BRQ_Kmer const& minKmer = **itr;
            size_t idx = itr - mEdgeEntries.begin();
            if ( CF<K>::getForm(mEdgeSeq.begin(idx))==CanonicalForm::REV )
            { mEdgeSeq.ReverseComplement();
                std::reverse(mEdgeEntries.begin(),mEdgeEntries.end());
                idx = mEdgeSeq.size()-idx-K;
                Assert(std::equal(minKmer.begin(),minKmer.end(),mEdgeSeq.begin(idx))); }
            if ( !idx )
                return;
            bvec bv; bv.reserve(mEdgeSeq.size());
            bv.assign(mEdgeSeq.begin(idx),mEdgeSeq.end());
            bv.append(mEdgeSeq.begin(K-1),mEdgeSeq.begin(K+idx-1));
            Assert(std::equal(mEdgeSeq.begin(),mEdgeSeq.begin(K-1),
                              mEdgeSeq.end()-(K-1)));
            mEdgeSeq = bv;
            mEdgeEntries.reserve(mEdgeEntries.size()+idx);
            std::copy(mEdgeEntries.begin(),mEdgeEntries.begin()+idx,
                      std::back_inserter(mEdgeEntries));
            mEdgeEntries.erase(mEdgeEntries.begin(),mEdgeEntries.begin()+idx); }

        bool isPalindrome( BRQ_Kmer const& kmer )
        { if ( !(K&1) )
                return kmer.isPalindrome();
            BRQ_SubKmer subKmer(kmer);
            if ( subKmer.isPalindrome() )
                return true;
            subKmer.toSuccessor(kmer.back());
            return subKmer.isPalindrome(); }

        bool upstreamExtensionPossible( BRQ_Entry const& entry )
        { KMerContext context = entry.getKDef().getContext();
            if ( context.getPredecessorCount() != 1 )
                return false;
            BRQ_Kmer pred(entry);
            pred.toPredecessor(context.getSinglePredecessor());
            if ( isPalindrome(pred) )
                return false;
            lookup(pred,&context);
            return context.getSuccessorCount() == 1; }

        bool downstreamExtensionPossible( BRQ_Entry const& entry )
        { KMerContext context = entry.getKDef().getContext();
            if ( context.getSuccessorCount() != 1 )
                return false;
            BRQ_Kmer succ(entry);
            succ.toSuccessor(context.getSingleSuccessor());
            if ( isPalindrome(succ) )
                return false;
            lookup(succ,&context);
            return context.getPredecessorCount() == 1; }

        void make1KmerEdge( BRQ_Entry const& entry )
        { mEdgeSeq.assign(entry.begin(),entry.end());
            mEdgeEntries.push_back(&entry);
            addEdge(); }

        void extendUpstream( BRQ_Entry const& entry )
        { mEdgeSeq.assign(entry.rcbegin(),entry.rcend());
            mEdgeEntries.push_back(&entry);
            extend(BRQ_Kmer(entry).rc(),entry.getKDef().getContext().rc()); }

        void extendDownstream( BRQ_Entry const& entry )
        { mEdgeSeq.assign(entry.begin(),entry.end());
            mEdgeEntries.push_back(&entry);
            extend(entry,entry.getKDef().getContext()); }

        void extend( BRQ_Kmer const& kmer, KMerContext context )
        { BRQ_Kmer next(kmer);
            while ( context.getSuccessorCount() == 1 )
            { unsigned char succCode = context.getSingleSuccessor();
                next.toSuccessor(succCode);
                if ( isPalindrome(next) )
                    break;
                BRQ_Entry const* pEntry = lookup(next,&context);
                if ( context.getPredecessorCount() != 1 )
                    break;
                mEdgeSeq.push_back(succCode);
                mEdgeEntries.push_back(pEntry); }
            switch (mEdgeSeq.getCanonicalForm())
            {case CanonicalForm::PALINDROME:
                    ForceAssertEq(mEdgeSeq.size(),K);
                    // allow flow-through to FWD case
                case CanonicalForm::FWD:
                    addEdge(); break;
                case CanonicalForm::REV:
                    mEdgeSeq.clear(); mEdgeEntries.clear(); break; } }

        BRQ_Entry const* lookup( BRQ_Kmer const& kmer, KMerContext* pContext )
        { BRQ_Entry const* result;
            if ( kmer.isRev() )
            { result = mDict.findEntryCanonical(BRQ_Kmer(kmer).rc());
                ForceAssert(result);
                *pContext = result->getKDef().getContext().rc(); }
            else
            { result = mDict.findEntryCanonical(kmer);
                ForceAssert(result);
                *pContext = result->getKDef().getContext(); }
            return result; }

        void addEdge()
        { static SpinLockedData gLock;
            EdgeID edgeID;
            if ( mEdgeSeq.getCanonicalForm() == CanonicalForm::REV )
            {
                mEdgeSeq.ReverseComplement();
                std::reverse(mEdgeEntries.begin(),mEdgeEntries.end());
            }
            if ( true )
            { SpinLocker lock(gLock);
                edgeID.setVal(mEdges.size());
                mEdges.push_back(mEdgeSeq); }
            unsigned offset = 0;
            bool err = false;
            for ( BRQ_Entry const* pEnt : mEdgeEntries )
            { KDef const& kDef = pEnt->getKDef();
                if ( kDef.isNull() )
                { Assert(std::equal(pEnt->begin(),pEnt->end(),mEdgeSeq.begin(offset)) ||
                         std::equal(pEnt->rcbegin(),pEnt->rcend(),mEdgeSeq.begin(offset)));
                    const_cast<KDef&>(kDef).set(edgeID,offset++); }
                else
                { std::cout << edgeID << ':' << offset++ << ' ' << pEnt
                  << " Already occupied as " << kDef.getEdgeID()
                  << ':' << kDef.getEdgeOffset() << std::endl;
                    err = true; } }
            if ( err )
                FatalErr("Having trouble with preoccupied kmers.");
            mEdgeSeq.clear();
            mEdgeEntries.clear(); }

        BRQ_Dict const& mDict;
        vecbvec& mEdges;
        std::vector<BRQ_Entry const*> mEdgeEntries;
        bvec mEdgeSeq;
    };

    void buildEdges( BRQ_Dict const& dict, vecbvec* pEdges )
    {
        EdgeBuilder eb(dict,pEdges);
        dict.parallelForEachHHS(
                [eb]( BRQ_Dict::Set::HHS const& hhs ) mutable
                { for ( BRQ_Entry const& entry : hhs )
                    if ( entry.getKDef().isNull() )
                        eb.buildEdge(entry); });

        size_t nRegularEdges = pEdges->size();
        size_t nTotalLength = pEdges->SizeSum();
        //std::cout << "There are " << nRegularEdges
        //          << " edges of total length " << nTotalLength << '.' << std::endl;

        // if a kmer isn't marked as being on an edge, it must be a part of a smooth
        // circle: add those edges, too.  simpleCircle method isn't thread-safe, so
        // this part is single-threaded.
        //std::cout << Date() << ": finding smooth circles." << std::endl;
        for ( auto const& hhs : dict )
            for ( auto const& entry : hhs )
                if ( entry.getKDef().isNull() )
                    eb.simpleCircle(entry);
        //std::cout << "There are " << pEdges->size()-nRegularEdges
        //          << " circular edges of total length "
        //          << pEdges->SizeSum()-nTotalLength << '.' << std::endl;
    }

    template <class Itr1, class Itr2>
    size_t matchLen( Itr1 itr1, Itr1 end1, Itr2 itr2, Itr2 end2 )
    {
        size_t result = 0;

        while ( itr1 != end1 && itr2 != end2 && *itr1 == *itr2 )
        { ++result; ++itr1; ++itr2; }

        return result;
    }


    template <class Itr1, class Itr2, class ItrW, typename sumType >
    size_t fuzzyMatchLen( Itr1 itr1, Itr1 end1, Itr2 itr2, Itr2 end2,
                          ItrW itrw, ItrW endw, sumType& maxWeight)
    {
        size_t result = 0;

        while ( itr1 != end1 && itr2 != end2 && itrw != endw )
        {
            if ( *itr1 != *itr2 ) maxWeight -= *itrw;

            if ( maxWeight < 0 ) break;             // EARLY EXIT!

            ++result; ++itr1; ++itr2; ++itrw;
        }

        return result;
    }

// fuzzyMatchLenBi - bidirectional matcher
//
// note the semantics:
//
// forward is the usual iterator traversal from its initial position
// up to, but not including, the end.
//
// backward first compares *(--itr) and continues backwards up to and
// including the passed in "end" position (which should be something like the
// beginning of the read or edge)
    template <class Itr1, class Itr2, class ItrW, typename sumType >
    size_t fuzzyMatchLenBi( Itr1 itr1, Itr1 end1, Itr2 itr2, Itr2 end2,
                            ItrW itrw, ItrW endw, sumType& maxWeight, bool backward=false)
    {
        size_t result = 0;
        int offset  = backward ? -1 : 0 ;
        int incr    = backward ? -1 : 1 ;

        while ( itr1 != end1 && itr2 != end2 && itrw != endw )
        {
            if ( itr1[offset] != itr2[offset] ) maxWeight -= *itrw;

            if ( maxWeight < 0 ) break;             // EARLY EXIT!

            ++result;
            itr1 += incr; itr2 += incr; itrw += incr;
        }

        return result;
    }


    class EdgeLoc
    {
    public:
        EdgeLoc() : mEdgeOffset(0) {}
        EdgeLoc( EdgeID const& edgeID, bool rc, unsigned edgeOffset )
                : mEdgeID(edgeID), mEdgeOffset(edgeOffset)
        { if ( rc ) mEdgeOffset = ~edgeOffset; }

        EdgeID const& getEdgeID() const { return mEdgeID; }
        bool isRC() const { return mEdgeOffset < 0; }
        unsigned getEdgeOffset() const
        { return mEdgeOffset < 0 ? ~mEdgeOffset : mEdgeOffset; }

        size_t hash() const { return (mEdgeID.val() << 32) | mEdgeOffset; }

        friend bool operator<( EdgeLoc const& el1, EdgeLoc const& el2 )
        { if ( el1.mEdgeID < el2.mEdgeID ) return true;
            if ( el2.mEdgeID < el1.mEdgeID ) return false;
            return el1.mEdgeOffset < el2.mEdgeOffset; }

        friend std::ostream& operator<<( std::ostream& os, EdgeLoc const& el )
        { if ( el.isRC() ) os << '~';
            return os << el.getEdgeID() << '[' << el.getEdgeOffset() << ']'; }

    private:
        EdgeID mEdgeID;
        int mEdgeOffset;
    };

    class PathPart
    {
    public:
        PathPart() : mLen(0u), mEdgeLen(0) {}
        PathPart( unsigned gapLen ) : mLen(gapLen), mEdgeLen(0) {}
        PathPart( EdgeID const& edgeID, bool rc, unsigned edgeOffset,
                  unsigned len, unsigned edgeLen )
                : mEdgeLoc(edgeID,rc,edgeOffset), mLen(len), mEdgeLen(edgeLen)
        { AssertLt(getEdgeOffset(),mEdgeLen); }

        EdgeID const& getEdgeID() const { return mEdgeLoc.getEdgeID(); }
        bool isRC() const { return mEdgeLoc.isRC(); }
        unsigned getEdgeOffset() const { return mEdgeLoc.getEdgeOffset(); }

        // if true, then only the length is meaningful
        bool isGap() const { return !mEdgeLen; }
        unsigned getLength() const { return mLen; }
        void incrLength(unsigned len) { mLen += len; }

        unsigned getEndOffset() const { return getEdgeOffset()+mLen; }

        // note: edge length is measured in kmers, not bases
        unsigned getEdgeLen() const { return mEdgeLen; }

        EdgeLoc firstLoc() const
        { return mEdgeLoc; }

        EdgeLoc lastLoc() const
        { EdgeLoc el(getEdgeID(),isRC(),getEdgeOffset()+mLen-1);
            AssertLt(el.getEdgeOffset(),mEdgeLen);
            return el; }

        bool isSameEdge( PathPart const& that ) const
        { return getEdgeID() == that.getEdgeID() && isRC() == that.isRC(); }

        bool isConformingCapturedGap( unsigned maxJitter ) const
        { PathPart const& prev = this[-1];
            PathPart const& next = this[1];
            unsigned graphDist = next.getEdgeOffset()-prev.getEndOffset();
            if ( !prev.isSameEdge(next) )
                graphDist += prev.getEdgeLen();
            // if the gap size is consistent with the graph, return true
            return unsigned(std::abs(int(mLen-graphDist))) <= maxJitter; }

        PathPart rc() const
        { return isGap() ? PathPart(mLen) :
                 PathPart(getEdgeID(),!isRC(),mEdgeLen-getEndOffset(),mLen,mEdgeLen); }

        friend std::ostream& operator<<( std::ostream& os, PathPart const& ppt )
        { if ( ppt.isGap() ) return os << ppt.getLength() << 'X';
            if ( ppt.isRC() ) os << '~';
            int offset = ppt.getEdgeOffset();
            return os << ppt.getEdgeID()
                   << '[' << offset << '-' << offset+ppt.getLength()
                   << ':' << ppt.getEdgeLen() << ']'; }

    private:
        EdgeLoc mEdgeLoc;
        unsigned mLen;
        unsigned mEdgeLen;
    };

    class BRQ_Pather
    {
    public:
        BRQ_Pather( BRQ_Dict const& dict, vecbvec const& edges )
                : mDict(dict), mEdges(edges) {}

        std::vector<PathPart> const& path( bvec const& read )
        { mPathParts.clear();

            if ( read.size() < K )
            { mPathParts.emplace_back(read.size());
                return mPathParts; } // EARLY RETURN!

            auto itr = read.begin();
            auto end = read.end()-K+1;
            while ( itr != end )
            { BRQ_Kmer kmer(itr);
                BRQ_Entry const* pEnt = mDict.findEntry(kmer);
                if ( !pEnt )
                {
                    unsigned gapLen = 1u;
                    auto itr2 = itr+K; ++itr;
                    auto end2 = read.end();
                    while ( itr2 != end2 )
                    { kmer.toSuccessor(*itr2); ++itr2;
                        if ( (pEnt = mDict.findEntry(kmer)) )
                            break;
                        ++gapLen; ++itr; }
                    mPathParts.emplace_back(gapLen); }
                if ( pEnt )
                { KDef const& kDef = pEnt->getKDef();
                    bvec const& edge = mEdges[kDef.getEdgeID().val()];
                    int offset = kDef.getEdgeOffset();
                    auto eBeg(edge.begin(offset));
                    size_t len = 1u;
                    bool rc = CF<K>::isRC(itr,eBeg);
                    if ( !rc )
                        len += matchLen(itr+K,read.end(),eBeg+K,edge.end());
                    else
                    { offset = edge.size() - offset;
                        auto eBegRC(edge.rcbegin(offset));
                        len += matchLen(itr+K,read.end(),eBegRC,edge.rcend());
                        offset = offset - K; }
                    unsigned edgeKmers = edge.size()-K+1;

                    mPathParts.emplace_back(kDef.getEdgeID(),rc,offset,len,edgeKmers);
                    itr += len;
                } }
            return mPathParts; }

        bool isJoinable( PathPart const& pp1, PathPart const& pp2 ) const
        { if ( pp1.getEdgeID() == pp2.getEdgeID() ) return true;
            bvec const& e1 = mEdges[pp1.getEdgeID().val()];
            bvec const& e2 = mEdges[pp2.getEdgeID().val()];
            BRQ_SubKmer k1 = pp1.isRC() ? BRQ_SubKmer(e1.rcend()-K+1) : BRQ_SubKmer(e1.end()-K+1);
            BRQ_SubKmer k2 = pp2.isRC() ? BRQ_SubKmer(e2.rcend()-K+1) : BRQ_SubKmer(e2.end()-K+1);
            return k1==k2; }

    private:
        BRQ_Dict const& mDict;
        vecbvec const& mEdges;
        std::vector<PathPart> mPathParts;
    };

    class GapFiller
    {
    public:
        GapFiller( vecbvec const& reads, vecbvec const& edges, unsigned maxGapSize,
                   unsigned minFreq, BRQ_Dict* pDict )
                : mReads(reads), mDict(*pDict), mPather(*pDict,edges),
                  mMaxGapSize(maxGapSize), mMinFreq(minFreq) {}

        template <class OItr>
        void map( size_t readId, OItr oItr )
        { bvec const& read = mReads[readId];
            std::vector<PathPart> const& parts = mPather.path(read);
            if ( parts.size() < 3 ) return; // EARLY RETURN!
            auto rItr = read.begin(parts[0].getLength());
            auto end = parts.end()-1;
            // for each path part sandwiched between two others
            for ( auto pPart=parts.begin()+1; pPart != end; ++pPart )
            { // if it's not a gap of appropriate size
                if ( !pPart->isGap() ||
                     (mMaxGapSize && pPart->getLength() > mMaxGapSize) ||
                     // or if it fits cleanly on the existing graph
                     pPart->isConformingCapturedGap(MAX_JITTER) )
                { rItr += pPart->getLength();
                    continue; } // just CONTINUE
                // we have a gap of appropriate size, and it has new, interesting kmers
                // that might change the graph structure
                auto itr = rItr - 1;
                BRQ_Kmer kkk(itr); itr += K;
                update(kkk,KMerContext::initialContext(*itr));
                auto last = itr + pPart->getLength();
                while ( itr != last )
                { unsigned char predCode = kkk.front();
                    kkk.toSuccessor(*itr);
                    KMerContext kc(predCode,*++itr);
                    *oItr++ = kkk.isRev()?BRQ_Entry(BRQ_Kmer(kkk).rc(),kc.rc()) : BRQ_Entry(kkk,kc); }
                KMerContext kc = KMerContext::finalContext(kkk.front());
                kkk.toSuccessor(*itr);
                update(kkk,kc);
                rItr += pPart->getLength(); } }

        void reduce( BRQ_Entry* e1, BRQ_Entry* e2 )
        { summarizeEntries(e1,e2);
            if ( e1->getKDef().getCount() >= mMinFreq )
            { mDict.insertEntry(std::move(*e1)); } }

        BRQ_Entry* overflow( BRQ_Entry* e1, BRQ_Entry* e2 )
        { if ( e2-e1 > 1 ) summarizeEntries(e1,e2); return e1+1; }

    private:
        void update( BRQ_Kmer kmer, KMerContext kc )
        { if ( kmer.isRev() ) { kmer.rc(); kc=kc.rc(); }
            mDict.applyCanonical(kmer,
                                 [kc](BRQ_Entry const& ent)
                                 { KDef& kd = const_cast<KDef&>(ent.getKDef());
                                     kd.setContext(kc|kd.getContext()); }); }

        static unsigned const MAX_JITTER = 1;
        vecbvec const& mReads;
        BRQ_Dict& mDict;
        BRQ_Pather mPather;
        unsigned mMaxGapSize;
        unsigned mMinFreq;
    };
    typedef MapReduceEngine<GapFiller,BRQ_Entry,BRQ_Kmer::Hasher> GFMRE;

    void fillGaps( vecbvec const& reads, unsigned maxGapSize, unsigned minFreq,
                   vecbvec* pEdges, BRQ_Dict* pDict )
    {
        //std::cout << Date() << ": filling gaps." << std::endl;
        GapFiller gf(reads,*pEdges,maxGapSize,minFreq,pDict);
        GFMRE mre(gf);
        if ( !mre.run(5*reads.size(),0ul,reads.size()) )
            FatalErr("Map/Reduce operation failed in gap filling.");

        //std::cout << "Now the dictionary has " << pDict->size()
        //                << " entries." << std::endl;

        //std::cout << Date() << ": finding edge sequences again." << std::endl;
        pDict->nullEntries();
        pDict->recomputeAdjacencies();
        pEdges->clear();
        buildEdges(*pDict,pEdges);
    }

    class BRQ_Join
    {
    public:
        BRQ_Join( EdgeLoc const& edgeLoc1, EdgeLoc const& edgeLoc2, unsigned overlap )
                : mEdgeLoc1(edgeLoc1), mEdgeLoc2(edgeLoc2), mOverlap(overlap)
        {}

        EdgeLoc const& getEdgeLoc1() const { return mEdgeLoc1; }
        EdgeLoc const& getEdgeLoc2() const { return mEdgeLoc2; }
        unsigned getOverlap() const { return mOverlap; }

        size_t hash() const
        { return 47 * (47*mEdgeLoc1.hash()+mEdgeLoc2.hash()) + mOverlap; }

        struct Hasher
        { size_t operator()( BRQ_Join const& join ) { return join.hash(); } };

        friend bool operator<( BRQ_Join const& join1, BRQ_Join const& join2 )
        { if ( join1.mEdgeLoc1 < join2.mEdgeLoc1 ) return true;
            if ( join2.mEdgeLoc1 < join1.mEdgeLoc1 ) return false;
            if ( join1.mEdgeLoc2 < join2.mEdgeLoc2 ) return true;
            if ( join2.mEdgeLoc2 < join1.mEdgeLoc2 ) return false;
            return join1.mOverlap < join2.mOverlap; }

        friend std::ostream& operator<<( std::ostream& os, BRQ_Join const& join )
        { return os << join.mEdgeLoc1 << "<-" << join.mOverlap << "->"
                 << join.mEdgeLoc2; }

    private:
        EdgeLoc mEdgeLoc1;
        EdgeLoc mEdgeLoc2;
        unsigned mOverlap;
    };

    class BRQ_Joiner
    {
    public:
        BRQ_Joiner( vecbvec const& reads, vecbvec const& edges, BRQ_Dict const& dict,
                unsigned maxGapSize, unsigned minFreq,
                vecbvec* pFakeReads )
                : mReads(reads), mEdges(edges), mPather(dict,edges),
                  mMaxGapSize(maxGapSize), mMinFreq(minFreq), mFakeReads(*pFakeReads)
        { ForceAssertLt(maxGapSize,K-1); }

        template <class OItr>
        void map( size_t readId, OItr oItr )
        { bvec const& read = mReads[readId];
            std::vector<PathPart> const& parts = mPather.path(read);
            if ( parts.size() < 3 ) return; // EARLY RETURN!
            auto end = parts.end()-1;
            for ( auto pPart=parts.begin()+1; pPart != end; ++pPart )
            { // if it's a captured gap of appropriate size
                if ( pPart->isGap() && pPart->getLength() <= mMaxGapSize )
                { PathPart const& prev = pPart[-1];
                    PathPart const& next = pPart[1];
                    unsigned overlap = K-pPart->getLength()-1;
                    if ( next.getEdgeID() < prev.getEdgeID() )
                        *oItr++ = BRQ_Join(next.rc().lastLoc(),prev.rc().firstLoc(),overlap);
                    else
                        *oItr++ = BRQ_Join(prev.lastLoc(),next.firstLoc(),overlap); } } }

        void reduce( BRQ_Join* pJoin1, BRQ_Join* pJoin2 )
        { Assert(validOverlap(*pJoin1));
            if ( pJoin2-pJoin1 >= mMinFreq )
            { //std::cout << "Joining: " << *pJoin1 << std::endl;
                bvec bv; bv.reserve(2*K);
                append(pJoin1->getEdgeLoc1(),0,&bv);
                append(pJoin1->getEdgeLoc2(),pJoin1->getOverlap(),&bv);
                add(bv); } }

        BRQ_Join* overflow( BRQ_Join* pJoin1, BRQ_Join* pJoin2 )
        { return pJoin1+std::min(unsigned(pJoin2-pJoin1),mMinFreq); }

    private:
        bool validOverlap( BRQ_Join const& join )
        { EdgeLoc const& el1 = join.getEdgeLoc1();
            bvec const& bv1 = mEdges[el1.getEdgeID().val()];
            EdgeLoc const& el2 = join.getEdgeLoc2();
            bvec const& bv2 = mEdges[el2.getEdgeID().val()];
            bool result;
            if ( el1.isRC() )
            { auto end = bv1.rcbegin(el1.getEdgeOffset()+K);
                auto beg = end - join.getOverlap();
                if ( el2.isRC() )
                    result = std::equal(beg,end,bv2.rcbegin(el2.getEdgeOffset()));
                else
                    result = std::equal(beg,end,bv2.begin(el2.getEdgeOffset())); }
            else
            { auto end = bv1.begin(el1.getEdgeOffset()+K);
                auto beg = end - join.getOverlap();
                if ( el2.isRC() )
                    result = std::equal(beg,end,bv2.rcbegin(el2.getEdgeOffset()));
                else
                    result = std::equal(beg,end,bv2.begin(el2.getEdgeOffset())); }
            return result; }

        void append( EdgeLoc const& el, unsigned indent, bvec* pBV )
        { bvec const& edge = mEdges[el.getEdgeID().val()];
            unsigned offset = el.getEdgeOffset()+indent;
            unsigned len = K-indent;
            if ( el.isRC() )
            { auto itr = edge.rcbegin(offset);
                pBV->append(itr,itr+len); }
            else
            { auto itr = edge.begin(offset);
                pBV->append(itr,itr+len); } }

        void add( bvec const& bv )
        { static SpinLockedData gLock;
            if ( true )
            { SpinLocker lock(gLock);
                mFakeReads.push_back(bv); } }

        vecbvec const& mReads;
        vecbvec const& mEdges;
        BRQ_Pather mPather;
        unsigned mMaxGapSize;
        unsigned mMinFreq;
        vecbvec& mFakeReads;
    };
    typedef MapReduceEngine<BRQ_Joiner,BRQ_Join,BRQ_Join::Hasher> JMRE;

    void joinOverlaps( vecbvec const& reads, unsigned maxGapSize, unsigned minFreq,
                       vecbvec* pEdges, BRQ_Dict* pDict )
    {
        //std::cout << Date() << ": joining overlaps." << std::endl;
        vecbvec fakeReads;
        fakeReads.reserve(pEdges->size()/10);
        BRQ_Joiner joiner(reads,*pEdges,*pDict,maxGapSize,minFreq,&fakeReads);
        JMRE mre(joiner);
        if ( !mre.run(reads.size(),0ul,reads.size()) )
            FatalErr("Map/Reduce operation failed when joining overlaps.");

        if ( fakeReads.size() )
        {
            //std::cout << "Found " << fakeReads.size() << " joins." << std::endl;
            pDict->nullEntries();
            pDict->process(fakeReads);
            //std::cout << "Now the dictionary has " << pDict->size()
            //                    << " entries." << std::endl;
            //std::cout << Date() << ": finding edge sequences again." << std::endl;
            pEdges->clear();
            buildEdges(*pDict,pEdges);
        }
    }


    inline size_t pathPartToEdgeID( PathPart const& part, std::vector<int> const& mFwdEdgeXlat, std::vector<int> const& mRevEdgeXlat )
    {
        ForceAssert(!part.isGap());
        size_t idx = part.getEdgeID().val();
        return ( part.isRC() ? mRevEdgeXlat[idx] : mFwdEdgeXlat[idx] );
    }

    inline void pathPartsToReadPath( std::vector<PathPart> const& parts, ReadPath& path, std::vector<int> const& mFwdEdgeXlat, std::vector<int> const& mRevEdgeXlat )
    {
        path.clear();
        PathPart const *pLast = 0;
        for (PathPart const &part : parts) {
            if (part.isGap()) continue;
            if (pLast && pLast->isSameEdge(part)) continue;
//        if ( part.getEdgeOffset() == 0 && part.getLength() < 5) continue;
            size_t idx = part.getEdgeID().val();
            path.push_back(part.isRC() ? mRevEdgeXlat[idx] : mFwdEdgeXlat[idx]);
            pLast = &part;
        }
        if (path.empty()) { path.setOffset(0); /*path.setLastSkip(0);*/ }
        else {
            PathPart const &firstPart = parts.front();
            if (!firstPart.isGap())
                path.setOffset(firstPart.getEdgeOffset());
            else {
                int eo1 = parts[1].getEdgeOffset();
                int eo2 = firstPart.getLength();
                path.setOffset(eo1 - eo2);
            }
        }
    }

    void path_reads_OMP( vecbvec const& reads, VecPQVec const& quals, BRQ_Dict const& dict, vecbvec const& edges,
            HyperBasevector const& hbv, std::vector<int> const& fwdEdgeXlat, std::vector<int> const& revEdgeXlat,
                     ReadPathVec* pPaths) {
        static unsigned const MAX_JITTER = 3;
        #pragma omp parallel
        {

            vec<int> toLeft,toRight;
            hbv.ToLeft(toLeft);
            hbv.ToLeft(toRight);
            BRQ_Pather mPather(dict,edges);
            ReadPath mPath;
            ExtendReadPath mExtender(hbv,&toLeft,&toRight);
            qvec mQV;

            #pragma omp for
            for (auto readId=0;readId<reads.size();++readId){
                std::vector<PathPart> parts = mPather.path(reads[readId]);     // needs to become a forward_list

                // convert any seeds on hanging edges to gaps
                std::vector<PathPart> new_parts;
                for ( auto part : parts ) {
                    if ( !part.isGap() ) {
                        size_t edge_id = pathPartToEdgeID(part,fwdEdgeXlat,revEdgeXlat);
                        size_t vleft  = toLeft[edge_id];
                        size_t vright = toRight[edge_id];
                        if ( hbv.ToSize(vleft) == 0
                             && hbv.ToSize(vright) > 1
                             && hbv.FromSize(vright) > 0
                             && part.getEdgeLen() <= 100 ) {
                            // delete a seed on a hanging edge
                            part = PathPart(part.getLength() );
                        }
                    }

                    // either add to an existing gap or create a new part (gap or not)
                    if ( part.isGap() && new_parts.size() && new_parts.back().isGap() )
                        new_parts.back().incrLength( part.getLength() );
                    else
                        new_parts.push_back(part);
                }
                std::swap(parts,new_parts);

                // if a gap is captured between two edges and things don't make sense, delete the
                // seed, if there are more than one, and subsequent seeds.  This hopefully avoids
                // the loss of sensitivity of just dropping the seeds.
                if ( parts.size() >= 3 ) {
                    size_t seeds = ( parts.begin()->isGap()) ? 0u : 1u;
                    for ( auto pPart = parts.begin()+1; pPart != parts.end()-1; ++pPart ) {
                        if ( !pPart->isGap() ) { seeds++; continue; }
                        if ( !pPart->isConformingCapturedGap(MAX_JITTER) ||
                             !mPather.isJoinable(pPart[-1],pPart[1]) ) {


                            if ( seeds > 1 ) {
                                PathPart tmpPart(pPart[-1].getLength());
                                for ( auto pPart2 = pPart; pPart2 != parts.end(); ++ pPart2 )
                                    tmpPart.incrLength( pPart2->getLength() );
                                parts.erase( pPart-1, parts.end() );
                                parts.push_back(tmpPart);
                            } else {
                                for ( auto pPart2 = pPart+1; pPart2 != parts.end(); ++ pPart2 )
                                    pPart->incrLength( pPart2->getLength() );
                                parts.erase( pPart+1, parts.end() );
                            }

                            break;
                        }
                    }
                }


                // if a seed has only taken us <= 5 kmers onto an edge, we back off
                // because we probably don't have enough evidence to commit.  Extension
                // code later will be more careful.
                if ( parts.back().isGap() && parts.size() > 1 ) {
                    auto const& last2 = parts[parts.size()-2];
                    if ( last2.getEdgeOffset() == 0 && last2.getLength() <= 5 ) {
                        auto last = parts.back();
                        last.incrLength(last2.getLength());
                        parts.pop_back();
                        parts.pop_back();
                        parts.push_back(last);
                    }
                } else if ( !parts.back().isGap() ) {
                    auto& last = parts.back();
                    if ( last.getEdgeOffset()==0 && last.getLength() <= 5 )  {
                        last = PathPart( last.getLength() );
                    }
                }

                pathPartsToReadPath(parts, mPath,fwdEdgeXlat,revEdgeXlat);

                quals[readId].unpack(&mQV);
                mExtender.attemptLeftRightExtension(mPath,reads[readId],mQV);

                (*pPaths)[readId] = mPath;
            }
        }

    }

} // end of anonymous namespace


inline void combine_Entries( BRQ_Entry & dest, BRQ_Entry & source )
{
    KDef& kDef = dest.getKDef();
    KMerContext kc=dest.getKDef().getContext();
    kc|=source.getKDef().getContext();
    kDef.setContext(kc);
    kDef.setCount(kDef.getCount()+source.getKDef().getCount());
}

void print_results_status(uint nthreads, bool results[]){
#pragma omp critical
    {
        std::cout<< "Results status: [ ";
        for (auto i=0;i<nthreads;++i) std::cout<<results[i]<<" ";
        std::cout<<" ]"<<std::endl;
    }
}

void merge_kmers_into_master(std::vector<BRQ_Entry> &kmer_list,unsigned workers, std::vector<BRQ_Entry> * partial_results,bool last_status[]) {
    //merge from the bottom

    //count the total number of partial results to process (ease out the looping).
    uint64_t results_left=0;
    bool active[workers];
    BRQ_Entry next_value[workers];
    for (auto i=0;i<workers;i++) {
        if (last_status[i]) {
            results_left += partial_results[i].size();
            active[i] = (partial_results[i].size() > 0);
        }
        else active[i]=false;
        if (active[i]) next_value[i]=partial_results[i].back();
    }


    if (results_left==0) return;
    //create all iterators and such
#pragma omp critical
    std::cout << Date() << ": Master merging " << results_left << " results, current master kmer list size "<< kmer_list.size()<< std::endl;
    kmer_list.reserve(kmer_list.size()+results_left);//is this wise?
    auto master_read=kmer_list.rbegin();

    kmer_list.resize(kmer_list.size()+results_left);//is this wise?

    auto master_write=kmer_list.rbegin();
    bool first=true;
    //while any partial results left


    while (results_left) {
        //find the max current partial result
        int max=-1;
        for (auto i=0;i<workers;i++)
            if (active[i]){
                if (-1==max or next_value[i]>next_value[max])
                    max=i;
            }
        //while the master has larger entries, insert them:
        while (master_read!=kmer_list.rend() and *master_read>=next_value[max]){
            if (first) first=false;
            else ++master_write;
            *master_write=*master_read;
            ++master_read;
        }
        //is the max > is the minimum = current last insertion?
        if (!first and partial_results[max].back()==*master_write) {
            combine_Entries(*master_write,next_value[max]);
            partial_results[max].pop_back();
        } else {
            if (first) first=false;
            else ++master_write;
            *master_write=next_value[max];
            partial_results[max].pop_back();
        }
        if (partial_results[max].size()==0) active[max]=false;
        else next_value[max]=partial_results[max].back();
        --results_left;
    }

    while (master_read!=kmer_list.rend()){//whatever is left from the master
        ++master_write;
        *master_write=*master_read;
        ++master_read;
    }
    //move to the top
    int64_t new_size=master_write+1-kmer_list.rbegin();
    std::cout<<"new master kmer list size:"<<new_size<<std::endl;
    kmer_list.assign(master_write.base()-1,kmer_list.end());
    kmer_list.resize(new_size);
}



BRQ_Dict * createDictOMP(vecbvec const& reads, VecPQVec const& quals, unsigned minQual, unsigned minFreq){
    //Count kmers per read, probably stupidly overkill
    std::vector<unsigned> goodLens(reads.size(),0);
    const uint64_t batch_size=100000000;
    //parallelForBatch(0ul,reads.size(),100000,
    //                 GoodLenTailFinder(quals,minQual,&goodLens));
    std::atomic_uint nKmers;
#pragma omp parallel
    {
        qvec uq;
        auto itr = uq.end();
        auto beg = uq.begin();
        unsigned good = 0;
#pragma omp for
        for (auto i = 0; i < reads.size(); ++i) {
            quals[i].unpack(&uq);
            itr = uq.end();
            beg = uq.begin();
            good = 0;
            while (itr != beg) {
                if (*--itr < minQual) good = 0;
                else if (++good == K) {
                    goodLens[i] = (itr - beg) + K;
                    break;
                }
            }
            nKmers += goodLens[i];

        }
    }

    uint64_t nThreads = omp_get_max_threads();

    std::atomic_uint blocks_to_do(0);
    std::atomic_bool unprocessed_results[nThreads];
    for (auto i=0;i<nThreads;i++) unprocessed_results[i]=false;
//    print_results_status(nThreads,unprocessed_results);
    std::vector<BRQ_Entry> partial_results[nThreads];
    std::vector<BRQ_Entry> kmer_list;
    auto start = 0;
    uint64_t task_kmers = 0;
    std::atomic_bool all_tasks_created(false);
    BRQ_Dict *pDict;
    #pragma omp parallel
    {
#pragma omp master
        {//The master thread, creates tasks first
#pragma omp critical
            std::cout << Date() << ": Master entering the main list merging loop with " << blocks_to_do
                      << " blocks to do " << std::endl;

            while (blocks_to_do or not all_tasks_created) {
                //copy the current status of the unprocessed_results so it is STABLE (we fix the groups to process as those in this iteration)
                bool last_status[nThreads];
                unsigned results_to_process = 0;//TODO: active wait is ugly
                for (auto i = 0; i < nThreads; ++i) {
                    last_status[i] = unprocessed_results[i];
                    if (last_status[i]) results_to_process++;
                }
                //TODO: replace this by a nice in-place merge sort.
                if (results_to_process) {
#pragma omp critical
                    std::cout << Date() << ": Master merging " << results_to_process << " batches into main list"
                              << std::endl;
                    merge_kmers_into_master(kmer_list, nThreads, partial_results, last_status);
#pragma omp critical
                    std::cout << Date() << ": Master merging done" << std::endl;

                    //decrease the blocks_to_do counter
                    for (auto i = 0; i < nThreads; ++i) if (last_status[i]) unprocessed_results[i] = false;
                    blocks_to_do -= results_to_process;
                    std::cout << Date() << ": Master still waiting for " << blocks_to_do << "blocks" << std::endl;
                }
            }
            std::cout << Date() << ": Master processed all blocks and tasks, creating dict... "<<std::endl;
            //now do the insertion on the dict, with no locking or anything, just a fast one, only if count>min_freq
            pDict = new BRQ_Dict(kmer_list.size());
            uint64_t used=0;
            for (auto &kent:kmer_list) {
                if (kent.getKDef().getCount() >= minFreq) {
                    pDict->insertEntry(kent);
                    used++;
                }
            }
            std::cout<< Date()<<": "<<used<<" out of "<<kmer_list.size()<<" kmers with Freq >= "<<minFreq<<std::endl;
            pDict->recomputeAdjacencies(); //does this mean we should not have saved adjacencies in the first place?

        }

#pragma omp single
        {
            for (auto i = 0; i <= reads.size(); ++i) {
                if (task_kmers + goodLens[i] > batch_size or i == reads.size()) {
                    //create the task to count the kmers
                    auto end = i;
                    //increase the task counter.
                    ++blocks_to_do;
    #pragma omp task shared(goodLens,reads,unprocessed_results,partial_results) firstprivate(start,end)
                    {
    #pragma omp critical
                        std::cout << Date() << ": Thread " << omp_get_thread_num() << " started a batch from " << start
                                  << " to " << end << std::endl;
                        //create the kmer count in a local vector
                        std::vector<BRQ_Entry> local_results;
                        local_results.reserve(task_kmers);
                        for (auto readId = start; readId < end; ++readId) {
                            unsigned len = goodLens[readId];
                            if (len > K) {
                                auto beg = reads[readId].begin(), itr = beg + K, last = beg + (len - 1);
                                BRQ_Kmer kkk(beg);
                                KMerContext kc = KMerContext::initialContext(*itr);

                                do {
                                    local_results.push_back(
                                            kkk.isRev() ? BRQ_Entry(BRQ_Kmer(kkk).rc(), kc.rc()) : BRQ_Entry(kkk, kc));
                                    unsigned char pred = kkk.front();
                                    kkk.toSuccessor(*itr);
                                    ++itr;
                                    kc = KMerContext(pred, *itr);

                                } while (itr <= last);

                            }
                        }
                        //sort and collapse
                        std::sort(local_results.begin(), local_results.end());
                        KMerContext kc;
                        auto okItr = local_results.begin();
                        for (auto kItr = local_results.begin(); kItr < local_results.end(); ++kItr) {
                            kc = kItr->getKDef().getContext();
                            uint count = 0;
                            while (kItr < local_results.end() - 1 and *kItr == *(kItr + 1)) {
                                KDef const &kDef = (kItr + 1)->getKDef();
                                kc |= kDef.getContext();
                                count += std::max(kDef.getCount(), 1ul);
                                ++kItr;
                            }
                            okItr = kItr;
                            KDef &kDef = okItr->getKDef();
                            kDef.setContext(kc);
                            kDef.setCount(count);
                            ++okItr;
                        }
                        local_results.resize(okItr - local_results.begin());

                        //wait for the std::atomic_bool "unprocessed results" flag to be false for this particular block
#pragma omp critical
                        std::cout << Date() << ": Thread " << omp_get_thread_num() << " waiting to send..."<< std::endl;
                        while (unprocessed_results[omp_get_thread_num()]);//TODO: active wait is ugly
                        //transfer results
                        partial_results[omp_get_thread_num()] = local_results;
                        unprocessed_results[omp_get_thread_num()] = true;
                        #pragma omp critical
                        std::cout << Date() << ": Thread " << omp_get_thread_num() << " sent batch to master" << std::endl;

                    }
                    //std::cout<<"Task created!!!"<<std::endl;
                    task_kmers = 0;
                    start = i;
                }

                if (i < reads.size()) task_kmers += goodLens[i];
            }
            all_tasks_created = true;
        }



    }


    return pDict;

}

void buildReadQGraph( vecbvec const& reads, VecPQVec const& quals,
                      bool doFillGaps, bool doJoinOverlaps,
                      unsigned minQual, unsigned minFreq,
                      double minFreq2Fract, unsigned maxGapSize,
                      HyperBasevector* pHBV, ReadPathVec* pPaths, int _K)
{
    std::cout << Date() << ": loading reads (createdictOMP!!!)." << std::endl;
    BRQ_Dict* pDict = createDictOMP(reads,quals,minQual,minFreq);

    std::cout << Date() << ": finding edge sequences." << std::endl;
    // figure out the complete base sequence of each edge
    vecbvec edges;
    edges.reserve(pDict->size()/100); //TODO: this is probably WAY too much in most scenarios
    buildEdges(*pDict,&edges);

    unsigned minFreq2 = std::max(2u,unsigned(minFreq2Fract*minFreq+.5));

    if ( doFillGaps ) { // Off by default
        std::cout << Date() << ": filling gaps." << std::endl;
        fillGaps(reads, maxGapSize, minFreq2, &edges, pDict);
    }

    if ( doJoinOverlaps ) { // Off by default
        std::cout << Date() << ": joining Overlaps." << std::endl;
        joinOverlaps(reads, _K / 2, minFreq2, &edges, pDict);
    }

    std::vector<int> fwdEdgeXlat;
    std::vector<int> revEdgeXlat;
    if ( !pPaths )
    {
        delete pDict;
        std::cout << Date() << ": building HBV." << std::endl;
        buildHBVFromEdges(edges,_K,pHBV,fwdEdgeXlat,revEdgeXlat);
    }
    else
    {
        std::cout << Date() << ": building HBV." << std::endl;
        buildHBVFromEdges(edges,K,pHBV,fwdEdgeXlat,revEdgeXlat);

        for (auto i=0;i<edges.size();++i) if (fwdEdgeXlat[i]<0 or fwdEdgeXlat[i]>pHBV->EdgeObjectCount()-1 or revEdgeXlat[i]<0 or revEdgeXlat[i]>pHBV->EdgeObjectCount()-1)
                std::cout<<"Edge "<<i<<" improper translation: "<<fwdEdgeXlat[i]<<" "<<revEdgeXlat[i]<<std::endl;
        std::cout << Date() << ": creating Read Paths." << std::endl;

        pPaths->clear().resize(reads.size());
        path_reads_OMP(reads, quals, *pDict, edges, *pHBV, fwdEdgeXlat, revEdgeXlat, pPaths);

        //parallelForBatch(0ul,reads.size(),100000,pather);

        delete pDict;
    }
    std::cout << Date() << ": graph created." << std::endl;
}