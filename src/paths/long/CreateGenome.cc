// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <sys/time.h>

#include "Basevector.h"
#include "ParallelVecUtilities.h"
#include "paths/HyperBasevector.h"
//#include "paths/LongReadTools.h"
#include "paths/long/CreateGenome.h"
#include "paths/long/LongProtoTools.h"


namespace
{
    int KmerId(const basevector& b, const int L, const int p)
    {
        int n = 0;
        for (int l = 0; l < L; l++) {
            n <<= 2;  // same as *= 4
            n += b[p + l];
        }
        return n;
    }

struct WorkItem
{
    WorkItem() : mRefId(0), mStart(0), mEnd(0), mDest(0) {}
    WorkItem( size_t refId, unsigned start, unsigned end, size_t dest )
    : mRefId(refId), mStart(start), mEnd(end), mDest(dest) {}

    // default copying, moving, and destructor are OK

    size_t mRefId;
    unsigned mStart;
    unsigned mEnd;
    size_t mDest;
};

typedef triple<int,int,int> IntTriple;

class KmerizationProc
{
public:
    KmerizationProc( vecbvec const& ref, int K, vec<IntTriple>& out )
    : mpRef(&ref), mpOut(&out), mK(K) {}

    // default copying, moving, and destructor are OK

    void operator()( WorkItem const& item ) const
    { bvec const& tig = (*mpRef)[item.mRefId];
      auto oItr = mpOut->begin()+item.mDest;
      for ( unsigned idx = item.mStart; idx != item.mEnd; ++idx )
      { *oItr = IntTriple(KmerId(tig,mK,idx),item.mRefId,idx); ++oItr; } }

private:
    vecbvec const* mpRef;
    vec<triple<int,int,int>>* mpOut;
    int mK;
};


} // end of anonymous namespace

void CreateGlocs( const vecbasevector& G, unsigned const LG, VecIntPairVec& Glocs )
{
    size_t nKmers = 0;
    for ( auto itr=G.begin(),end=G.end(); itr != end; ++itr )
        if ( itr->size() >= LG )
            nKmers += itr->size()-LG+1;
    vec<IntTriple> out(nKmers);
    if ( true )
    {
        KmerizationProc proc(G,LG,out);
        Worklist<WorkItem,KmerizationProc> wl(proc);
        size_t nTigs = G.size();
        size_t dest = 0;
        unsigned const BATCH_SIZE = 100000;
        ForceAssertGe(BATCH_SIZE,LG);
        for ( size_t idx = 0; idx != nTigs; ++idx )
        {
            unsigned off = 0;
            unsigned sz = G[idx].size();
            if ( sz < LG )
                continue;
            sz = sz - LG + 1;
            while ( sz >= BATCH_SIZE )
            {
                wl.add(WorkItem(idx,off,off+BATCH_SIZE,dest));
                off += BATCH_SIZE;
                dest += BATCH_SIZE;
                sz -= BATCH_SIZE;
            }
            if ( sz )
            {
                wl.add(WorkItem(idx,off,off+sz,dest));
                dest += sz;
            }
        }
    }
    ParallelSort(out);
    Glocs.clear().resize( 1 << (2*LG) );
    auto beg = out.begin(), end = out.end();
    while ( beg != end )
    {
        auto itr = beg;
        int kmerId = beg->first;
        while ( ++itr != end )
            if ( kmerId != itr->first )
                break;
        IntPairVec& locs = Glocs[kmerId];
        locs.reserve(itr-beg);
        while ( beg != itr )
        {
            locs.push_back(std::pair<int,int>(beg->second,beg->third));
            ++beg;
        }
    }
}
