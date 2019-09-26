/*
 * BuildReadQGraph.cc
 *
 *  Created on: Jan 22, 2014
 *      Author: tsharpe
 */

#include "paths/long/BuildReadQGraph.h"
#include "Basevector.h"
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
#include <util/OutputLog.h>
#include "paths/HyperBasevector.h"
#include "paths/long/ExtendReadPath.h"

namespace
{

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

    class EdgeBuilder {
    public:
        EdgeBuilder(BRQ_Dict const &dict, vecbvec *pEdges)
                : mDict(dict), mEdges(*pEdges) {}

        void buildEdge(BRQ_Entry const &entry) {
            if (isPalindrome(entry))
                make1KmerEdge(entry);
            else if (upstreamExtensionPossible(entry)) {
                if (downstreamExtensionPossible(entry))
                    return;
                extendUpstream(entry);
            } else if (downstreamExtensionPossible(entry))
                extendDownstream(entry);
            else
                make1KmerEdge(entry);
        }

        bool isEnd(BRQ_Entry const &entry) {
            if (isPalindrome(entry))
                return true;
            else if (upstreamExtensionPossible(entry) and downstreamExtensionPossible(entry))
                    return false;
            return true;
        }

        // not thread-safe
        void simpleCircle(BRQ_Entry const &entry) {
            BRQ_Entry const *pFirstEntry = &entry;
            mEdgeSeq.assign(entry.begin(), entry.end());
            mEdgeEntries.push_back(pFirstEntry);
            KMerContext context = entry.getKDef().getContext();
            BRQ_Kmer kmer(entry);
            while (true) {
                ForceAssertEq(context.getPredecessorCount(), 1u);
                ForceAssertEq(context.getSuccessorCount(), 1u);
                unsigned char succCode = context.getSingleSuccessor();
                kmer.toSuccessor(succCode);
                BRQ_Entry const *pEntry = lookup(kmer, &context);
                if (pEntry == pFirstEntry)
                    break;
                if (!pEntry->getKDef().isNull()) {
                    std::cout << "Failed to close circle.\n";
                    for (auto beg = mEdgeEntries.begin(), end = mEdgeEntries.end();
                         beg != end; ++beg)
                        std::cout << *beg << ' ' << **beg << '\n';
                    std::cout << pEntry << ' ' << *pEntry << std::endl;
                    CRD::exit(1);
                }
                mEdgeSeq.push_back(succCode);
                mEdgeEntries.push_back(pEntry);
            }
            canonicalizeCircle();
            addEdge();
        }

    private:
        void canonicalizeCircle() {
            auto itr = std::min_element(mEdgeEntries.begin(), mEdgeEntries.end(),
                                        [](BRQ_Kmer const *pEnt1, BRQ_Kmer const *pEnt2) { return *pEnt1 < *pEnt2; });
            BRQ_Kmer const &minKmer = **itr;
            size_t idx = itr - mEdgeEntries.begin();
            if (CF<K>::getForm(mEdgeSeq.begin(idx)) == CanonicalForm::REV) {
                mEdgeSeq.ReverseComplement();
                std::reverse(mEdgeEntries.begin(), mEdgeEntries.end());
                idx = mEdgeSeq.size() - idx - K;
                Assert(std::equal(minKmer.begin(), minKmer.end(), mEdgeSeq.begin(idx)));
            }
            if (!idx)
                return;
            bvec bv;
            bv.reserve(mEdgeSeq.size());
            bv.assign(mEdgeSeq.begin(idx), mEdgeSeq.end());
            bv.append(mEdgeSeq.begin(K - 1), mEdgeSeq.begin(K + idx - 1));
            Assert(std::equal(mEdgeSeq.begin(), mEdgeSeq.begin(K - 1),
                              mEdgeSeq.end() - (K - 1)));
            mEdgeSeq = bv;
            mEdgeEntries.reserve(mEdgeEntries.size() + idx);
            std::copy(mEdgeEntries.begin(), mEdgeEntries.begin() + idx,
                      std::back_inserter(mEdgeEntries));
            mEdgeEntries.erase(mEdgeEntries.begin(), mEdgeEntries.begin() + idx);
        }

        bool isPalindrome(BRQ_Kmer const &kmer) {
            if (!(K & 1))
                return kmer.isPalindrome();
            BRQ_SubKmer subKmer(kmer);
            if (subKmer.isPalindrome())
                return true;
            subKmer.toSuccessor(kmer.back());
            return subKmer.isPalindrome();
        }

        bool upstreamExtensionPossible(BRQ_Entry const &entry) {
            KMerContext context = entry.getKDef().getContext();
            if (context.getPredecessorCount() != 1)
                return false;
            BRQ_Kmer pred(entry);
            pred.toPredecessor(context.getSinglePredecessor());
            if (isPalindrome(pred))
                return false;
            lookup(pred, &context);
            return context.getSuccessorCount() == 1;
        }

        bool downstreamExtensionPossible(BRQ_Entry const &entry) {
            KMerContext context = entry.getKDef().getContext();
            if (context.getSuccessorCount() != 1)
                return false;
            BRQ_Kmer succ(entry);
            succ.toSuccessor(context.getSingleSuccessor());
            if (isPalindrome(succ))
                return false;
            lookup(succ, &context);
            return context.getPredecessorCount() == 1;
        }

        void make1KmerEdge(BRQ_Entry const &entry) {
            mEdgeSeq.assign(entry.begin(), entry.end());
            mEdgeEntries.push_back(&entry);
            addEdge();
        }

        void extendUpstream(BRQ_Entry const &entry) {
            mEdgeSeq.assign(entry.rcbegin(), entry.rcend());
            mEdgeEntries.push_back(&entry);
            extend(BRQ_Kmer(entry).rc(), entry.getKDef().getContext().rc());
        }

        void extendDownstream(BRQ_Entry const &entry) {
            mEdgeSeq.assign(entry.begin(), entry.end());
            mEdgeEntries.push_back(&entry);
            extend(entry, entry.getKDef().getContext());
        }

        void extend(BRQ_Kmer const &kmer, KMerContext context) {
            BRQ_Kmer next(kmer);
            while (context.getSuccessorCount() == 1) {
                unsigned char succCode = context.getSingleSuccessor();
                next.toSuccessor(succCode);
                if (isPalindrome(next))
                    break;
                BRQ_Entry const *pEntry = lookup(next, &context);
                if (context.getPredecessorCount() != 1)
                    break;
                mEdgeSeq.push_back(succCode);
                mEdgeEntries.push_back(pEntry);
            }
            switch (mEdgeSeq.getCanonicalForm()) {
                case CanonicalForm::PALINDROME:
                    ForceAssertEq(mEdgeSeq.size(), K);
                    // allow flow-through to FWD case
                case CanonicalForm::FWD:
                    addEdge();
                    break;
                case CanonicalForm::REV:
                    mEdgeSeq.clear();
                    mEdgeEntries.clear();
                    break;
            }
        }

        BRQ_Entry const *lookup(BRQ_Kmer const &kmer, KMerContext *pContext) {
            BRQ_Entry const *result;
            if (kmer.isRev()) {
                result = mDict.findEntryCanonical(BRQ_Kmer(kmer).rc());
                ForceAssert(result);
                *pContext = result->getKDef().getContext().rc();
            } else {
                result = mDict.findEntryCanonical(kmer);
                ForceAssert(result);
                *pContext = result->getKDef().getContext();
            }
            return result;
        }

        void addEdge() {
            static SpinLockedData gLock;
            EdgeID edgeID;
            if (mEdgeSeq.getCanonicalForm() == CanonicalForm::REV) {
                mEdgeSeq.ReverseComplement();
                std::reverse(mEdgeEntries.begin(), mEdgeEntries.end());
            }
            if (true) {
                SpinLocker lock(gLock);
                edgeID.setVal(mEdges.size());
                mEdges.push_back(mEdgeSeq);
            }
            unsigned offset = 0;
            bool err = false;
            for (BRQ_Entry const *pEnt : mEdgeEntries) {
                KDef const &kDef = pEnt->getKDef();
                if (kDef.isNull()) {
                    Assert(std::equal(pEnt->begin(), pEnt->end(), mEdgeSeq.begin(offset)) ||
                           std::equal(pEnt->rcbegin(), pEnt->rcend(), mEdgeSeq.begin(offset)));
                    const_cast<KDef &>(kDef).set(edgeID, offset++);
                } else {
                    std::cout << edgeID << ':' << offset++ << ' ' << pEnt
                              << " Already occupied as " << kDef.getEdgeID()
                              << ':' << kDef.getEdgeOffset() << std::endl;
                    err = true;
                }
            }
            if (err)
                FatalErr("Having trouble with preoccupied kmers.");
            mEdgeSeq.clear();
            mEdgeEntries.clear();
        }

        BRQ_Dict const &mDict;
        vecbvec &mEdges;
        std::vector<BRQ_Entry const *> mEdgeEntries;
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
        OutputLog(2) << nRegularEdges
                  << " edges of total length " << nTotalLength << std::endl;

        // if a kmer isn't marked as being on an edge, it must be a part of a smooth
        // circle: add those edges, too.  simpleCircle method isn't thread-safe, so
        // this part is single-threaded.
        //std::cout << Date() << ": finding smooth circles." << std::endl;
        for ( auto const& hhs : dict )
            for ( auto const& entry : hhs )
                if ( entry.getKDef().isNull() )
                    eb.simpleCircle(entry);
        OutputLog(2) << pEdges->size()-nRegularEdges
                  << " circular edges of total length "
                  << pEdges->SizeSum()-nTotalLength << std::endl;
    }

    template <class Itr1, class Itr2>
    size_t matchLen( Itr1 itr1, Itr1 end1, Itr2 itr2, Itr2 end2 )
    {
        size_t result = 0;

        while ( itr1 != end1 && itr2 != end2 && *itr1 == *itr2 )
        { ++result; ++itr1; ++itr2; }

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
        {
            mPathParts.clear();
            if (read.size() < K) {
                mPathParts.emplace_back(read.size());
                return mPathParts;
            } // EARLY RETURN!

            auto itr = read.begin();
            auto end = read.end() - K + 1;
            while (itr != end) {
                BRQ_Kmer kmer(itr);
                BRQ_Entry const *pEnt = mDict.findEntry(kmer);
                if (!pEnt) {
                    unsigned gapLen = 1u;
                    auto itr2 = itr + K;
                    ++itr;
                    auto end2 = read.end();
                    while (itr2 != end2) {
                        kmer.toSuccessor(*itr2);
                        ++itr2;
                        if ((pEnt = mDict.findEntry(kmer)))
                            break;
                        ++gapLen;
                        ++itr;
                    }
                    mPathParts.emplace_back(gapLen);
                }
                if (pEnt) {
                    KDef const &kDef = pEnt->getKDef();
                    bvec const &edge = mEdges[kDef.getEdgeID().val()];
                    int offset = kDef.getEdgeOffset();
                    auto eBeg(edge.begin(offset));
                    size_t len = 1u;
                    bool rc = CF<K>::isRC(itr, eBeg);
                    if (!rc)
                        len += matchLen(itr + K, read.end(), eBeg + K, edge.end());
                    else {
                        offset = edge.size() - offset;
                        auto eBegRC(edge.rcbegin(offset));
                        len += matchLen(itr + K, read.end(), eBegRC, edge.rcend());
                        offset = offset - K;
                    }
                    unsigned edgeKmers = edge.size() - K + 1;

                    mPathParts.emplace_back(kDef.getEdgeID(), rc, offset, len, edgeKmers);
                    itr += len;
                }
            }
            return mPathParts;
        }

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

    void path_reads_OMP( vecbvec const& reads, BRQ_Dict const& dict, vecbvec const& edges,
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

            #pragma omp for
            for (size_t readId=0;readId<reads.size();++readId){
                //this brings parts as sections of perfect matches to an edge. gaps are represented as 0 with len 0 on edge, and gap size on read
                std::vector<PathPart> parts = mPather.path(reads[readId]);     // needs to become a forward_list

                // convert any seeds on TIPS smaller than 100bp to gaps, then add them to previous gaps if needed.
                //TODO: XXX shouldnt the tip size be read sie (i.e. 250)
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

                (*pPaths)[readId].insert((*pPaths)[readId].end(),mPath.begin(),mPath.end());
            }
        }

    }

    void improve_read_paths_OMP( vecbvec const& reads, VecPQVec const& quals, HyperBasevector const& hbv, ReadPathVec & pPaths) {
#pragma omp parallel
        {
            vec<int> toLeft,toRight;
            hbv.ToLeft(toLeft);
            hbv.ToLeft(toRight);
            QualVec mQV;
            ExtendReadPath mExtender(hbv,&toLeft,&toRight);
#pragma omp for
            for (size_t readId=0;readId<reads.size();++readId){
                quals[readId].unpack(&mQV);
                mExtender.attemptLeftRightExtension(pPaths[readId],reads[readId],mQV);
            }
        }

    }

} // end of anonymous namespace




void create_read_lengths(std::vector<uint16_t> & rlen, VecPQVec const& quals, unsigned minQual){

    uint64_t qsize=quals.size();
    rlen.resize(qsize);
    #pragma omp parallel shared(rlen,quals,qsize)
    {
        QualVec uq;
        #pragma omp for
        for (uint64_t i = 0; i < qsize; ++i) {
            quals[i].unpack(&uq);
            auto itr = uq.end();
            auto beg = uq.begin();
            uint16_t good = 0;
            while (itr != beg) {
                if (*--itr < minQual) good = 0;
                else if (++good == K) {
                    rlen[i] = (itr - beg) + K;
                    break;
                }
            }
        }
    }
    OutputLog(2) << "Read lengths created"<<std::endl;
}

//a kmerlist class, replaces an std::vector with faster methods for the merge-sort
KmerList::~KmerList() {
    clear();
}

void KmerList::clear() {
    if (nullptr != kmers and 0 != size) {
        free(kmers);
        size = 0;
        kmers = nullptr;
    }
}

void KmerList::merge(KmerList &other) {
    //std::cout <<"merging batches of size "<<size<<" and "<<other.size<<std::endl;
    auto itr1 = kmers, end1 = kmers + size;
    auto itr2 = other.kmers, end2 = other.kmers + other.size, itr2w = other.kmers;
    //std::cout<<"merging, accumulating on counts1"<<std::endl;
    while (itr2 != end2) {
        while (itr1 != end1 and *itr1 < *itr2) ++itr1;
        if (itr1 != end1 and *itr1 == *itr2) {
            //combine_Entries(*itr1,*itr2);
            itr1->combine(*itr2);
            ++itr1;
            ++itr2;
        }
        while (itr2 != end2 and (itr1 == end1 or *itr2 < *itr1)) {
            *itr2w = *itr2;
            ++itr2w;
            ++itr2;
        }
    }
    //shrink counts2 to the size of its remaining elements
    //std::cout<<" equal elements merge done, now resizing and merging "<< (itr2w - other.kmers) << " new elements "<<std::endl;
    //other.resize(itr2w - other.kmers);
    other.size=itr2w - other.kmers;
    //counts2->shrink_to_fit();
    //expand counts1 to allow the insertion of the unique values on counts2
    //std::cout<<"merging, resizing counts1 to " << counts1.size()+counts2.size() << std::endl;
    auto old_size=size;
    resize(size + other.size);
    auto ritr1 = kmers + old_size - 1;

    //merge-sort from the bottom into count1.
    //std::cout<<"merging, final merging"<<std::endl;
    auto writr1 = kmers + size - 1, rend1 = kmers;
    auto ritr2 = other.kmers + other.size - 1, rend2 = other.kmers;

    while (writr1 >= rend1) {
        if (ritr2 >= rend2 and (ritr1 < rend1 or *ritr2 > *ritr1)) {
            memcpy(writr1, ritr2, sizeof(KMerNodeFreq_s));
            --ritr2;
        } else {
            memcpy(writr1, ritr1, sizeof(KMerNodeFreq_s));
            --ritr1;
        }
        --writr1;
    }
    other.clear();
    //std::cout<<"merge done, new size "<<size<<std::endl;
}


void KmerList::sort() {
    std::sort(kmers,kmers+size);
}

void KmerList::uniq() {
    auto wptr = kmers;
    auto endptr = kmers + size;
    for (auto rptr = kmers; rptr < endptr; ++wptr) {
        *wptr = *rptr;
        ++rptr;
        while (rptr < endptr and *wptr == *rptr) {
            wptr->combine(*rptr);
            ++rptr;
        }
    }
    resize(wptr - kmers);
}

void KmerList::dump(std::string filename) {
    std::ofstream batch_file(filename, std::ios::out | std::ios::trunc | std::ios::binary);
    uint64_t total_kmers = size;
    batch_file.write((const char *) &total_kmers, sizeof(uint64_t));
    batch_file.write((const char *) kmers, sizeof(KMerNodeFreq_s) * size);
    batch_file.close();
}

void KmerList::load(std::string filename) {
    std::ifstream batch_file(filename, std::ios::in | std::ios::binary);
    uint64_t total_kmers;
    batch_file.read((char *) &total_kmers, sizeof(uint64_t));
    resize(total_kmers);
    OutputLog(3)<<"Reading "<<size<<" kmers"<<std::endl;
    batch_file.read((char *) kmers, sizeof(KMerNodeFreq_s) * size);
    batch_file.close();
}

void KmerList::resize(size_t new_size) {
    if (size != new_size) {
        //std::cout << " allocating space for "<< new_size <<" elements: " << sizeof(KMerNodeFreq_s) * new_size <<std::endl;
        if (0==size) kmers = (KMerNodeFreq_s *) malloc(sizeof(KMerNodeFreq_s) * new_size);
        else kmers = (KMerNodeFreq_s *) realloc(kmers, sizeof(KMerNodeFreq_s) * new_size);
        //if (new_size>0 and kmers == nullptr) std::cout << " realloc error!!! "<<std::endl;
        size = new_size;
    }
}



std::shared_ptr<KmerList> kmerCountOMP(vecbvec const& reads, std::vector<uint16_t> const &rlen,
                                                               uint64_t gfrom, uint64_t gto, uint64_t batch_size=0) {
    std::atomic<uint64_t> totalKmers(0);
    //Compute how many "batches" will be used. and malloc a structure for them and a bool array to mark them "ready to process".
    //optionally, there could be a total count of "ready kmers" to set a limit for memory usage, if total count is too high, slaves would wait
    if (batch_size == 0) batch_size = (gto - gfrom) / (4 * omp_get_max_threads()) + 1;
    const uint64_t batches = ((gto - gfrom) + batch_size - 1) / batch_size;
    OutputLog(3) << "OMP-merge kmer counting in " << batches << " batches of " << batch_size << " reads" << std::endl;

    std::vector<std::atomic_uint_fast8_t *> batch_status; //level->batch->status

    std::vector<std::vector<std::shared_ptr<KmerList>>> batch_lists; //level->batch-> pointer to batch
    uint16_t levels = 0;
    for (auto elements = batches; elements > 1; elements = (elements + 1) / 2) {
        batch_status.push_back((std::atomic_uint_fast8_t *) calloc(sizeof(std::atomic_uint_fast8_t), elements));
        batch_lists.push_back(std::vector<std::shared_ptr<KmerList > >());
        batch_lists.back().resize(elements);
        OutputLog(4) << "level " << levels << " created with " << elements << " elements " << std::endl;
        ++levels;
    }
    std::atomic_uint_fast8_t level_count[levels]; //level->count
    for (auto &l:level_count)l = 0;

#pragma omp parallel shared(rlen,reads)
    {
#pragma omp for schedule(dynamic)
        for (auto batch = 0; batch < batches; ++batch) {
            //==== Part 1: actually counting ====
            uint64_t from = gfrom + batch * batch_size;
            uint64_t read_count = (batch < batches - 1) ? batch_size : (gto - gfrom) - batch_size * (batches - 1);
            uint64_t to = from + read_count;
            std::shared_ptr<KmerList> local_kmer_list = std::make_shared<KmerList>();
            uint64_t total_good_lenght = std::accumulate(rlen.begin() + from, rlen.begin() + to, 0);
            local_kmer_list->resize(total_good_lenght);
            //Populate the kmer list
            uint64_t last_kmer=0;
            for (auto readId = from; readId < to; ++readId) {
                unsigned len = rlen[readId];
                if (len > K) {
                    auto beg = reads[readId].begin(), itr = beg + K, last = beg + (len - 1);
                    KMerNodeFreq kkk(beg);
                    kkk.hash();
                    kkk.kc = KMerContext::initialContext(*itr);
                    kkk.count = 1;
                    (kkk.isRev() ? KMerNodeFreq(kkk, true) : kkk).to_struct( local_kmer_list->kmers[last_kmer] );
                    ++last_kmer;
                    while (itr != last) {
                        unsigned char pred = kkk.front();
                        kkk.toSuccessor(*itr);
                        ++itr;
                        kkk.kc = KMerContext(pred, *itr);
                        (kkk.isRev() ? KMerNodeFreq(kkk, true) : kkk).to_struct( local_kmer_list->kmers[last_kmer] );
                        ++last_kmer;
                    }
                    kkk.kc = KMerContext::finalContext(kkk.front());
                    kkk.toSuccessor(*last);
                    (kkk.isRev() ? KMerNodeFreq(kkk, true) : kkk).to_struct( local_kmer_list->kmers[last_kmer] );
                    ++last_kmer;
                }
            }
            totalKmers += last_kmer;
            //std::cout<<last_kmer<<" kmers inserted in a batch"<<std::endl;
            local_kmer_list->resize(last_kmer);
            local_kmer_list->sort();
            local_kmer_list->uniq();
            //std::cout<<local_kmer_list->size<<" unique kmers in a batch"<<std::endl;
            //std::sort(local_kmer_list->begin(), local_kmer_list->end());
            //collapse_entries(*local_kmer_list);

            //==== Part 2: merging till no results of same level available ====
            //merge_level=0
            uint16_t merge_level = 0;
            bool just_merged = true;
            //while (just merged)
            while (just_merged) {
                //   insert the list into this merge level's queue
                uint16_t slot = level_count[merge_level]++;
                if (slot == batch_lists[merge_level].size() - 1) {
#pragma omp critical
                    OutputLog(4) << "level " << merge_level << " done." << std::endl;
                }
                //   if insertion number is odd or last batch:
                if (merge_level < levels - 1 and (slot % 2 == 1 or slot == batch_lists[merge_level].size() - 1)) {
                    // mix with the previous even number (using own list as base, this way we preserve locality)
                    if (slot % 2 == 1) {
                        while (batch_status[merge_level][slot - 1] != 1)
                            usleep(10); //wait for the previous batch to be finished.
                        batch_status[merge_level][slot - 1] = 2;
                        batch_status[merge_level][slot] = 2;
                        //inplace_count_merge(local_kmer_list, batch_lists[merge_level][slot - 1]);
                        local_kmer_list->merge(*batch_lists[merge_level][slot - 1]);
                        batch_lists[merge_level][slot - 1].reset();
                    }
                    //      increase level
                    ++merge_level;
                    just_merged = true;
                } else {
                    // no batch available for merging, next thread at this level will pick up
                    batch_lists[merge_level][slot] = local_kmer_list;
                    local_kmer_list.reset();
                    batch_status[merge_level][slot] = 1;
                    just_merged = false;
                }
                //IDEA: fixed partitioning can be used to keep 2*thread-number lists for longer.
            }
        }

    }


    OutputLog(3) << "Top level merge starting" << std::endl;
    //inplace_count_merge(batch_lists.back()[0], batch_lists.back()[1]);
    batch_lists.back()[0]->merge(*batch_lists.back()[1]);
    OutputLog(3) << "Top level merge done" << std::endl;
    std::shared_ptr<KmerList> kmer_list = batch_lists.back()[0];
    batch_lists.back()[0].reset();
    batch_lists.back()[1].reset();
    for (auto &bs:batch_status) free(bs);
    OutputLog(3) << "cleanup done" << std::endl;
    OutputLog(3) << "Total kmers processed " << totalKmers << std::endl;
    return kmer_list;
}




std::shared_ptr<KmerList> kmerCountOMPDiskBased(vecbvec const& reads, std::vector<uint16_t> const &rlen, unsigned minCount,
                                                                   std::string tmpdir, std::string workdir, unsigned char disk_batches=8, uint64_t batch_size=0) {

    //If size larger than batch (or still not enough cpus used, or whatever), Lauch 2 tasks to sort the 2 halves, with minFreq=0

    OutputLog(2) << "disk-based kmer counting with "<<(int) disk_batches<<" batches"<<std::endl;
    uint64_t total_kmers_in_batches=0;
    for (auto batch=0;batch < disk_batches;batch++) {
        uint64_t to = (batch+1) * reads.size()/disk_batches;
        uint64_t from= batch * reads.size()/disk_batches;
        auto kmer_list = kmerCountOMP(reads,rlen,from,to);
        std::ofstream batch_file(tmpdir+"/kmer_count_batch_"+std::to_string((int)batch),std::ios::out | std::ios::trunc | std::ios::binary);
        batch_file.write((const char *)kmer_list->kmers,sizeof(KMerNodeFreq_s)*kmer_list->size);
        batch_file.close();
        OutputLog(2) << "batch "<<(int) batch<<" done and dumped with "<<kmer_list->size<< " kmers" <<std::endl;
        total_kmers_in_batches+=kmer_list->size;
    }

    //now a multi-merge sort between all batch files into the Dict
    OutputLog(2) << "merging from disk"<<std::endl;
    //open all batch files
    std::ifstream dbf[disk_batches];
    bool dbf_active[disk_batches];
    KMerNodeFreq_s next_knf_from_dbf[disk_batches];
    uint finished_files=0;
    for (auto i=0;i<disk_batches;i++){
        dbf[i].open(tmpdir+"/kmer_count_batch_"+std::to_string((int)i),std::ios::in | std::ios::binary);
        dbf[i].read((char *)&next_knf_from_dbf[i],sizeof(KMerNodeFreq_s));
        dbf_active[i]=true;
    }
    //set all finished flags to false
    //TODO: stupid minimum search
    KMerNodeFreq_s current_kmer;
    //while finished_files<batches
    bool first=true;
    uint min=0;
    for (auto i=1;i<disk_batches;++i)
        if (dbf_active[i]){
            if (next_knf_from_dbf[i]<next_knf_from_dbf[min]) min=i;
        }
    current_kmer=next_knf_from_dbf[min];
    current_kmer.count=0;
    uint64_t used = 0,not_used=0;
    uint64_t hist[256];
    for (auto &h:hist) h=0;
    std::shared_ptr<KmerList> kmerlist=std::make_shared<KmerList>();
    //size_t last_kmer=0;
    size_t alloc_block=10000000;
    while (finished_files<disk_batches) {
        //find minimum of the non-finished files
        uint min=0;
        while (!dbf_active[min])++min;
        for (auto i=1;i<disk_batches;++i)
            if (dbf_active[i]){
                if (next_knf_from_dbf[i]<next_knf_from_dbf[min]) min=i;
            }
        //larger than current kmer?
        if (next_knf_from_dbf[min] > current_kmer) {
            ++hist[std::min(255,(int)current_kmer.count)];
            if (current_kmer.count>=minCount) {
                //(*dict)->insertEntryNoLocking(BRQ_Entry((BRQ_Kmer) current_kmer, current_kmer.kc));
                //kmerlist.push_back(current_kmer);
                if (used>=kmerlist->size) kmerlist->resize(kmerlist->size+alloc_block);
                kmerlist->kmers[used]=current_kmer;
                ++used;
            }
            else ++not_used;
            current_kmer=next_knf_from_dbf[min];
        } else {
            current_kmer.combine(next_knf_from_dbf[min]);
        }
        //advance min file
        dbf[min].read((char *)&next_knf_from_dbf[min],sizeof(KMerNodeFreq_s));
        if ( dbf[min].eof() ) {
            dbf_active[min]=false;
            ++finished_files;
        }
    }

    ++hist[(int)current_kmer.count];
    if (current_kmer.count>=minCount) {
        //kmerlist.push_back(current_kmer);
        if (used>=kmerlist->size) kmerlist->resize(kmerlist->size+alloc_block);
        kmerlist->kmers[used]=current_kmer;
        used++;
    }
    else not_used++;
    kmerlist->resize(used);
    for (auto i=0;i<disk_batches;i++) {
        dbf[i].close();
        std::remove((tmpdir + "/kmer_count_batch_" +std::to_string((int)i)).c_str());
    }
    OutputLog(2)<< used << "/" << used+not_used << " kmers with Freq >= " << minCount << std::endl;
    if (""!=workdir) {
        std::ofstream kff(workdir + "/small_K.freqs");
        for (auto i = 1; i < 256; i++) kff << i << ", " << hist[i] << std::endl;
        kff.close();
    }
    return kmerlist;
}


std::shared_ptr<KmerList> buildKMerCount( vecbvec const& reads,
                                                                std::vector<uint16_t> & rlen, unsigned minCount,
                                                                std::string workdir, std::string tmpdir,
                                                                unsigned char disk_batches, uint64_t count_batch_size )
{
    std::shared_ptr<KmerList> spectrum=std::make_shared<KmerList>();
    OutputLog(2) << "creating kmers from reads..." << std::endl;
    if (1 >= disk_batches) {
        uint64_t hist[256]={0};
        uint64_t used=0;
        spectrum=kmerCountOMP(reads, rlen, 0, rlen.size(), count_batch_size);
        auto witr=spectrum->kmers;
        auto send=spectrum->kmers+spectrum->size;
        for (auto itr=spectrum->kmers;itr!=send;++itr){
            ++hist[itr->count];
            if (itr->count>=minCount) {
                (*witr)=(*itr);
                ++witr;
                ++used;
            }

        }
        OutputLog(2)<< used << "/" << spectrum->size << " kmers with Freq >= " << minCount << std::endl;
        spectrum->resize(used);
        std::ofstream kff(workdir + "/small_K.freqs");
        for (auto i = 1; i < 256; i++) kff << i << ", " << hist[i] << std::endl;
        kff.close();
        //todo: filter kmers to minCount
    } else {
        if ("" == tmpdir) tmpdir = workdir;
        spectrum=kmerCountOMPDiskBased(reads, rlen, minCount, tmpdir, workdir, disk_batches, count_batch_size);
    }

    return spectrum;
}

void dumpkmers( std::shared_ptr<std::vector<KMerNodeFreq_s>> const kmercounts, std::string filename) {
    std::ofstream batch_file(filename,std::ios::out | std::ios::trunc | std::ios::binary);
    uint64_t total_kmers=kmercounts->size();
    batch_file.write((const char *)&total_kmers,sizeof(uint64_t));
    batch_file.write((const char *)kmercounts->data(),sizeof(KMerNodeFreq_s)*kmercounts->size());
    batch_file.close();
}

std::shared_ptr<std::vector<KMerNodeFreq_s>> loadkmers( std::string filename) {
    std::shared_ptr<std::vector<KMerNodeFreq_s>> kmercounts;
    kmercounts=std::make_shared<std::vector<KMerNodeFreq_s>>();
    std::ifstream batch_file(filename,std::ios::in | std::ios::binary);
    uint64_t total_kmers;
    batch_file.read((char *)&total_kmers,sizeof(uint64_t));
    kmercounts=std::make_shared<std::vector<KMerNodeFreq_s>>();
    kmercounts->resize(total_kmers);
    batch_file.read(( char *)kmercounts->data(),sizeof(KMerNodeFreq_s)*kmercounts->size());
    batch_file.close();
    return kmercounts;
}

void buildReadQGraph( std::string out_dir,
                      bool doFillGaps, bool doJoinOverlaps,
                      unsigned minFreq, double minFreq2Fract, unsigned maxGapSize,  HyperBasevector* pHBV,
                      ReadPathVec* pPaths, int _K)
{


    uint64_t numKmers(0),usedKmers(0);
    std::FILE* kmers_from_disk;
    kmers_from_disk = std::fopen(std::string(out_dir+"/raw_kmers.data").data(), "rb");
    if (!kmers_from_disk) {
        std::perror("Failed to open raw_kmers.data, ");
    }
    OutputLog(2) << "Creating dict." << std::endl;
    BRQ_Dict * pDict = new BRQ_Dict(numKmers);
    OutputLog(2) << "Filtering kmers into Dict..." << std::endl;
    std::fread(&numKmers, sizeof(numKmers), 1, kmers_from_disk);
    KMerNodeFreq_s reads_kmer;
    //while next(a) or next(b)
    for (auto read_bytes = std::fread(&reads_kmer, sizeof(reads_kmer), 1, kmers_from_disk);
         0<read_bytes;
         read_bytes = std::fread(&reads_kmer, sizeof(reads_kmer), 1, kmers_from_disk)) {
        if (reads_kmer.count>=minFreq) {
            KMerNodeFreq knf(reads_kmer);
            pDict->insertEntryNoLocking(BRQ_Entry((BRQ_Kmer) knf, knf.kc));
            ++usedKmers;
        }

    }
    OutputLog(2) << usedKmers << "/" << numKmers <<" kmers with freq >= "<< minFreq << std::endl;
    pDict->recomputeAdjacencies();
    OutputLog(2) << "finding edges (unique paths)" << std::endl;
    // figure out the complete base sequence of each edge
    vecbvec edges;
    edges.reserve(pDict->size()/100); //TODO: this is probably WAY too much in most scenarios
    buildEdges(*pDict,&edges);
    uint64_t totalk=0;
    for (auto &e:edges) totalk+=e.size()+1-_K;
    OutputLog(2) <<edges.size()<<" edges with "<<totalk<<" "<<_K<<"-mers"<<std::endl;

    unsigned minFreq2 = std::max(2u,unsigned(minFreq2Fract*minFreq+.5));
    vecbvec reads;
    VecPQVec quals;
    if ( doFillGaps ) { // Off by default
        OutputLog(2) << "filling gaps." << std::endl;
        fillGaps(reads, maxGapSize, minFreq2, &edges, pDict);
    }

    if ( doJoinOverlaps ) { // Off by default
        OutputLog(2) << "joining Overlaps." << std::endl;
        joinOverlaps(reads, _K / 2, minFreq2, &edges, pDict);
    }

    std::vector<int> fwdEdgeXlat;
    std::vector<int> revEdgeXlat;
    if ( !pPaths )
    {
        delete pDict;
        OutputLog(2) << "building graph..." << std::endl;
        buildHBVFromEdges(edges,_K,pHBV,fwdEdgeXlat,revEdgeXlat);
        OutputLog(2) << "graph built" << std::endl;

    }
    else
    {
        OutputLog(2) << "building graph..." << std::endl;
        buildHBVFromEdges(edges,_K,pHBV,fwdEdgeXlat,revEdgeXlat);

        if (reads.size()==0) {
            OutputLog(2) << "Loading bases..." << std::endl;
            reads.ReadAll(out_dir + "/pe_data.fastb");
        }


        // TODO: Cleanup tips shorter than K+5, these won't be used for pathing the reads anyway
        {
            vec<int> to_left, to_right;
            pHBV->ToLeft(to_left), pHBV->ToRight(to_right);

            vec<int> dels;

            // Cleanup of "From" edges
            for (int v = 0; v < pHBV->N(); v++) {
                int num_from(pHBV->From(v).size());
                if (num_from > 0) {
                    for (int ei = 0; ei < num_from; ei++) {
                        int e = pHBV->EdgeObjectIndexByIndexFrom(v, ei);
                        if (pHBV->EdgeObject(e).size() < pHBV->K()+5 && pHBV->FromSize(to_right[e]) == 0) {
                            dels.emplace_back(e);
                        }
                    }
                }
            }

            // Cleanup of "To" edges
            for (int v = 0; v < pHBV->N(); v++) {
                int num_to(pHBV->To(v).size());
                if (num_to > 0) {
                    for (int ei = 0; ei < num_to; ei++) {
                        int e = pHBV->EdgeObjectIndexByIndexTo(v, ei);
                        if (pHBV->EdgeObject(e).size() < pHBV->K()+5 && pHBV->ToSize(to_left[e]) == 0) {
                            dels.emplace_back(e);
                        }
                    }
                }
            }

            UniqueSort(dels);
            pHBV->DeleteEdges(dels);
        }

        pHBV->RemoveUnneededVertices();

        // Recalculate pDict!
        {
            std::vector<KMerNodeFreq> kmers;
            kmers.reserve(numKmers);
            for (int edge=0; edge < pHBV->E(); edge++) {
                auto edgeObject(pHBV->EdgeObject(edge));
                unsigned len = edgeObject.size();
                if (len > K) {
                    auto beg = edgeObject.begin(), itr = beg + K, last = beg + (len - 1);
                    KMerNodeFreq kkk(beg);
                    kkk.hash();
                    kkk.kc = KMerContext::initialContext(*itr);
                    kkk.count = 1;
                    (kkk.isRev()) ? kmers.emplace_back(kkk, true): kmers.emplace_back(kkk);

                    while (itr != last) {
                        unsigned char pred = kkk.front();
                        kkk.toSuccessor(*itr);
                        ++itr;
                        kkk.kc = KMerContext(pred, *itr);
                        (kkk.isRev()) ? kmers.emplace_back(kkk, true): kmers.emplace_back(kkk);

                    }
                    kkk.kc = KMerContext::finalContext(kkk.front());
                    kkk.toSuccessor(*last);
                    (kkk.isRev()) ? kmers.emplace_back(kkk, true): kmers.emplace_back(kkk);
                }
            }
            std::sort(kmers.begin(), kmers.end());
            auto iter = std::unique(kmers.begin(), kmers.end());
            // Now v becomes {1 2 3 7 8 10 * * * * * *}
            // * means undefined

            // Resizing the vector so as to remove the undefined terms
            kmers.resize(std::distance(kmers.begin(), iter));

            delete pDict;
            pDict = new BRQ_Dict(kmers.size());
            for (const auto & k : kmers) {
                pDict->insertEntryNoLocking(BRQ_Entry ( (BRQ_Kmer)k, k.kc ));
            }
        }
        // Recalculate edges
        {
            edges.clear();
            for (int e = 0; e < pHBV->E(); e++) {
                edges.push_back(pHBV->EdgeObject(e));
            }
        }

        pPaths->clear();
        pPaths->resize(reads.size());
        OutputLog(2) << "pathing "<<reads.size()<<" reads into graph..." << std::endl;
        path_reads_OMP(reads, *pDict, edges, *pHBV, fwdEdgeXlat, revEdgeXlat, pPaths);
        delete pDict;
        edges.clear(); edges.shrink_to_fit();revEdgeXlat.clear();fwdEdgeXlat.clear();

        if (quals.size()==0) {
            OutputLog(2) << "Loading quals..." << std::endl;
            load_quals(quals, out_dir + "/pe_data.cqual");
        }
        improve_read_paths_OMP(reads,quals,*pHBV,*pPaths);
        OutputLog(2) << "reads pathed"<<std::endl;
        vec<int> to_right, to_left;
        pHBV->ToRight(to_right);
        pHBV->ToLeft(to_left);
        auto multi = 0, fixed = 0;
        for (auto pi = 0; pi < pPaths->size(); ++pi) {
            auto &p = (*pPaths)[pi];
            if (p.size() < 2) continue;
            ++multi;
            for (auto i = 1; i < p.size(); ++i) {
                if (to_right[p[i - 1]] != to_left[p[i]]) {
                    //std::cout<<"Path "<<pi<<" has a false connection "<<p[i-1]<<" -> "<<p[i]<<std::endl;
                    p.resize(i - 1);
                    ++fixed;
                }
            }
        }
        OutputLog(2) << "checking/fixing paths done, " << fixed << "/" << multi << " paths with multiple edges fixed" << std::endl;
    }

}