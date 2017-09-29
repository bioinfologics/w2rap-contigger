//
// Created by Luis Yanes (EI) on 29/09/2017.
//

#include "BuildReadQGraph.h"

std::ostream& operator<<(std::ostream& os, const KMerNodeFreq_s& kmer) {
    os << kmer.kdata[0] << '-' << kmer.kdata[1];
    os << "\t" << (int)kmer.count << ' ' << std::bitset<8>(kmer.kc);
    return os;
}

std::istream& operator>>(std::istream& is, KMerNodeFreq_s& kmer) {
    is.read((char*)&kmer, sizeof(kmer));
    return is;
}
