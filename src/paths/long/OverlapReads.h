#ifndef OVERLAP_READS_H
#define OVERLAP_READS_H

#include "CoreTools.h"
#include "Basevector.h"
#include "efasta/EfastaTools.h"
// The utilities to overlay a small number of efasta reads ( < 100 ) and create
// a linear assembly. 

class ReadOverlap{
public:
    int rid1, rid2;        // read ids
    int alt1, alt2;        // which alternatives from the two reads
    int overlap;           // how much is the overlap
    ReadOverlap(int r1, int r2, int a1, int a2, int o):
        rid1(r1),rid2(r2),alt1(a1),alt2(a2),overlap(o) {};
    ReadOverlap(){ rid1=-1; rid2=-1;};
    bool operator== ( const ReadOverlap& r ) const {
        return rid1 == r.rid1 && rid2 == r.rid2 &&
            alt1 == r.alt1 && alt2 == r.alt2 && overlap == r.overlap ;
    }
    friend std::ostream& operator<<( std::ostream& out, const ReadOverlap& r);
};

inline std::ostream& operator<< ( std::ostream& out, const ReadOverlap& r) {
    return out << r.rid1 << "_" << r.alt1 << "( " <<  r.overlap <<" )" << r.rid2 << "_" << r.alt2;
}

inline std::ostream& operator<<  ( std::ostream& out, const vec<ReadOverlap>& rs) {
    out << "[";
    for ( size_t i = 0; i < rs.size(); ++i ) out << rs[i] << ", ";
    out << "]";
    return out;
}




#endif
