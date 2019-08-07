#include "CoreTools.h"
#include "graph/DigraphTemplate.h"
#include "paths/long/KmerCount.h"

template void RemoveHangingEnds<kmer_count>(digraphE<kmer_count>&, 
     int (kmer_count::*)() const, int, double);

template void RemoveHangingEnds3<kmer_count>(digraphE<kmer_count>&, 
     int (kmer_count::*)() const, int, double, int);

template void digraphE<kmer_count>::Used(vec<Bool>&) const;
template void DistancesToEndArr( const digraphE<kmer_count>&, vec<int> const&, int const, Bool const, vec<int>& );
