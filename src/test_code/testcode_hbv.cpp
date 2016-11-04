//
// Created by Gonzalo Garcia (TGAC) on 28/10/2016.
//

#include <iostream>
#include <typeinfo>
#include "testcode_hbv.h"

#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "kmers/BigKMer.h"
#include "feudal/HashSet.h"
#include "kmers/KMerHasher.h"

template <unsigned BIGK>
using BigDict = HashSet<BigKMer<BIGK>,typename BigKMer<BIGK>::hasher>;

template <unsigned BIGK>
class BigKMerizer
{
public:
    typedef BigDict<BIGK> BigKDict;
    BigKMerizer( BigKDict* pDict ) : mDict(*pDict) {}

    void kmerize( bvec const& bv )
    { if ( bv.size() < BIGK ) return;
      KMerHasher<BIGK> hasher;
      auto itr = bv.begin();
      size_t hash = hasher.hash(itr);
      if ( bv.size() == BIGK ) { add(BigKMer<BIGK>(bv,hash)); return; }
      // [GONZA]: Creates Kmer --> first (initial context means there is no predecessor, there is no previous base)
      BigKMer<BIGK> kmer(bv,hash,KMerContext::initialContext(bv[BIGK]));
      add(kmer);
      auto end = bv.end()-BIGK;
      while ( ++itr != end )
      { hash = hasher.stepF(itr);
        kmer.successor(hash,KMerContext(itr[-1],itr[BIGK]));
        add(kmer); }
      hash = hasher.stepF(itr);
      kmer.successor(hash,KMerContext::finalContext(itr[-1]));
      add(kmer); }

    template <class OItr>
    void map( vecbvec::const_iterator iItr, OItr oItr )
    { bvec const& bv = *iItr;
      if ( bv.size() < BIGK ) return;
      KMerHasher<BIGK> hasher;
      auto itr = bv.begin();
      size_t hash = hasher.hash(itr);
      if ( bv.size() == BIGK )
      { BigKMer<BIGK> kmer(bv,hash);
        *oItr = kmer.isRev() ? kmer.rc() : kmer; ++oItr;
        return; }
      BigKMer<BIGK> kmer(bv,hash,KMerContext::initialContext(bv[BIGK]));
      *oItr = kmer.isRev() ? kmer.rc() : kmer; ++oItr;
      auto end = bv.end()-BIGK;
      while ( ++itr != end )
      { hash = hasher.stepF(itr);
        kmer.successor(hash,KMerContext(itr[-1],itr[BIGK]));
        *oItr = kmer.isRev() ? kmer.rc() : kmer; ++oItr; }
      hash = hasher.stepF(itr);
      kmer.successor(hash,KMerContext::finalContext(itr[-1]));
      *oItr = kmer.isRev() ? kmer.rc() : kmer; ++oItr; }

    void reduce( BigKMer<BIGK>* key1, BigKMer<BIGK>* key2 )
    { for ( auto itr=key1+1; itr != key2; ++itr )
        key1->addContext(itr->getContext());
      mDict.insertUniqueValue(*key1); }

    // We build the dictionary with BigKMers pointing to reads.  This updates
    // all the BigKMers in the dictionary to point to edges instead.
    void updateDict( bvec const& edge )
    { KMerHasher<BIGK> hasher;
      auto itr = edge.begin();
      BigKMer<BIGK> kmer(edge,hasher.hash(itr));
      update(kmer.isRev()?kmer.rc():kmer);
      auto end = edge.end()-BIGK+1;
      while ( ++itr != end )
      { kmer.successor(hasher.stepF(itr));
        update(kmer.isRev()?kmer.rc():kmer); } }

private:
    void add( BigKMer<BIGK> const& kmer )
    { canonicalAdd( kmer.isRev() ? kmer.rc() : kmer ); }

    void update( BigKMer<BIGK> const& kmer )
    { const_cast<BigKMer<BIGK>*>(mDict.lookup(kmer))->updateLocation(kmer); }

    void canonicalAdd( BigKMer<BIGK> const& kmer ) const
    { mDict.apply(kmer,
                  [&kmer]( BigKMer<BIGK> const& entry )
                  { const_cast<BigKMer<BIGK>&>(entry).addContext(kmer.getContext()); }); }

    BigKDict& mDict;
};

template <unsigned BIGK>
class Pather
{
public:
    typedef BigDict<BIGK> BigKDict;

    Pather( vecbvec const& reads, BigKDict const& dict, vecbvec const& edges,
            std::vector<int> const& fwdXlat, std::vector<int> const& revXlat,
            ReadPathVec* pReadPaths,
            HyperKmerPath const* pHKP, vecKmerPath* pKmerPaths )
        : mReads(reads), mDict(dict), mEdges(edges),
          mFwdXlat(fwdXlat), mRevXlat(revXlat), mpReadPaths(pReadPaths),
          mpHKP(pHKP), mpKmerPaths(pKmerPaths) {}

    void operator()( size_t readId ) {
      bvec const &read = mReads[readId];
      if (read.size() < BIGK) return;
      KMerHasher<BIGK> hasher;
      BigKMer<BIGK> kmer(read, hasher(read.begin()));
      BigKMer<BIGK> entry = lookup(kmer);
      mTmpReadPath.setFirstSkip(entry.getOffset());
      size_t idx = &entry.getBV() - &mEdges[0];
      size_t edgeId = entry.isRC() ? mRevXlat[idx] : mFwdXlat[idx];
      mTmpReadPath.push_back(edgeId);
      size_t readLenRemaining = read.size();
      size_t edgeLenRemaining = entry.getBV().size() - entry.getOffset();
      while (readLenRemaining > edgeLenRemaining) {
        readLenRemaining = readLenRemaining - edgeLenRemaining + BIGK - 1;
        unsigned readOffset = read.size() - readLenRemaining;
        BigKMer<BIGK> nextKmer(read, hasher(read.begin(readOffset)),
                               KMerContext(), readOffset);
        BigKMer<BIGK> nextEntry = lookup(nextKmer);
        ForceAssertEq(nextEntry.getOffset(), 0u);
        idx = &nextEntry.getBV() - &mEdges[0];
        edgeId = nextEntry.isRC() ? mRevXlat[idx] : mFwdXlat[idx];
        mTmpReadPath.push_back(edgeId);
        edgeLenRemaining = nextEntry.getBV().size();
      }
      edgeLenRemaining -= readLenRemaining;
      /* mTmpReadPath.setLastSkip(edgeLenRemaining); */
      if (mpReadPaths)
        (*mpReadPaths)[readId] = mTmpReadPath;
      if (mpKmerPaths) {
        buildKmerPath(edgeLenRemaining);
        (*mpKmerPaths)[readId] = mTmpKmerPath;
        mTmpKmerPath.clear();
      }
      mTmpReadPath.clear();
    }

private:
    BigKMer<BIGK> lookup(BigKMer<BIGK> const &kmer) {
      if (kmer.isRev()) {
        BigKMer<BIGK> const *pEnt = mDict.lookup(kmer.rc());
        ForceAssert(pEnt);
        return pEnt->rc();
      }
      BigKMer<BIGK> const *pEnt = mDict.lookup(kmer);
      ForceAssert(pEnt);
      return *pEnt;
    }

    void buildKmerPath(size_t edgeLenRemaining) {
      if (mTmpReadPath.size() == 1) {
        KmerPath const &edge = mpHKP->EdgeObject(mTmpReadPath.front());
        mTmpKmerPath.Assign(edge.front().Start() + mTmpReadPath.getFirstSkip(),
                            edge.back().Stop() - edgeLenRemaining);
      } else {
        KmerPath const &edge1 = mpHKP->EdgeObject(mTmpReadPath.front());
        mTmpKmerPath.Assign(edge1.front().Start() + mTmpReadPath.getFirstSkip(),
                            edge1.back().Stop());
        auto end = mTmpReadPath.end() - 1;
        for (auto itr = mTmpReadPath.begin() + 1; itr != end; ++itr) {
          KmerPath const &edge = mpHKP->EdgeObject(*itr);
          SerfVec<KmerPathInterval> &tmp = mTmpKmerPath;
          tmp.append(edge.begin(), edge.end());
        }
        KmerPath const &edgeN = mpHKP->EdgeObject(*end);
        kmer_id_t stop = edgeN.back().Stop() - edgeLenRemaining;
        mTmpKmerPath.Append(edgeN.front().Start(), stop);
      }
    }

    vecbvec const &mReads;
    BigKDict const &mDict;
    vecbvec const &mEdges;
    std::vector<int> const &mFwdXlat;
    std::vector<int> const &mRevXlat;
    ReadPathVec *mpReadPaths;
    HyperKmerPath const *mpHKP;
    vecKmerPath *mpKmerPaths;
    ReadPath mTmpReadPath;
    KmerPath mTmpKmerPath;
};

hbv_explorer::hbv_explorer(std::string hbv_filename){


  // CARGAR LAS LECTURAS
  vecbvec reads;
  reads.ReadAll("/Users/ggarcia/Documents/test_dataset/test/PE1_frag_reads_orig.fastb");
  std::cout << "Loaded reads: " << reads.size() << std::endl;

  HyperBasevector hbv;
  BinaryReader::readFile(hbv_filename, &hbv);

  auto edges = hbv.Edges();
  std::cout << edges.size() << " " << edges[0] << std::endl;

  vec<int> to_left;
  vec<int> to_right;

  hbv.ToLeft(to_left);
  hbv.ToRight(to_right);

  std::cout << "Hbv loaded..." << hbv.EdgeObjectCount() << " edges!"<<std::endl;
//  for (auto a = 0; a<hbv.EdgeObjectCount(); ++a){
//    std::cout<< "Edge #: " << a << " Edge size: " << hbv.EdgeObject(a).isize()
//    << " toleft: " << to_left[a]<< " to rigth: " << to_right[a]  << std::endl;
//  }

  using BigDict = HashSet<BigKMer<200>,typename BigKMer<200>::hasher>;
  BigDict bigDict(3000000);
  BigKMerizer<200> tkmerizer(&bigDict);

  for (auto a = 0; a<hbv.EdgeObjectCount(); ++a){
    tkmerizer.kmerize(edges[a]);
  }

//  for (auto a = 0; a<hbv.EdgeObjectCount(); ++a){
//    // Udatea el diccionario con los ejes del grafo en lugar de con las lecturas
//    tkmerizer.updateDict(edges[a]);
//  }

  for (auto nr = 0; nr<20; ++nr){

    auto read = reads[nr];
    std::set<uint64_t > mapped_edges;

    KMerHasher<200> hasher;
    auto itr = read.begin();
    auto end = read.end()-200;

    size_t hash = hasher.hash(itr);
//    BigKMer<200> kmer(read,hash,KMerContext::initialContext(read[200]));
    BigKMer<200> kmer(read, hasher(itr));

    int cont = 0;
    while (++itr != end){
      hash = hasher.stepF(itr);
//      kmer.successor(hash,KMerContext(itr[-1],itr[200]));
      kmer.successor(hash);

      auto fkmer = bigDict.lookup(kmer);
      if (fkmer != 0){
        std::cout<<"Read #: " << nr << std::endl;
        std::cout << "Contador: " << cont
                  << " Kmer: " << (fkmer->isRC() ? kmer.rc() : kmer)
                  << " Offset: " << fkmer->getOffset()
                  << " Is RC: " << fkmer->isRC()
                  << " Context: " << fkmer->getContext()
                  << " BV size: " << fkmer->getBV().size()
                  << " BV # : " << &fkmer->getBV() - &edges[0]
                  << std::endl<< std::endl
                  << " " << hbv.EdgeObject(&fkmer->getBV() - &edges[0])
                  << std::endl;
        cont ++;
      }
    }
  }
}