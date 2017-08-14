//
// Created by Luis Yanes (EI) on 09/08/2017.
//

#include "paths/long/BuildReadQGraph.h"

std::ostream& operator<<(std::ostream& os, const KMerNodeFreq_s& kmer)
{
    os << kmer.kdata[0] << '-' << kmer.kdata[1] << "\n";
    os << kmer.count << ' ' << kmer.kc << "\n";
    return os;
}

std::istream& operator>> (std::istream& is, KMerNodeFreq_s& kmer)
{
    is.read((char*)&kmer, sizeof(kmer));
    return is;
}


int main(int argc, char **argv) {
    std::ifstream inf(argv[1], std::ios_base::binary|std::ios_base::in);
    if (!inf) {
        return -1;
    }
    uint64_t totalKmers(0), fileKmers(0);
    inf.read((char *) &totalKmers, sizeof(totalKmers));
    std::cout << "Final count " << totalKmers << "\n";
    KMerNodeFreq_s mer, prev;
    inf >> mer;
    std::cout << mer;
    prev = mer;
    while (!inf.eof()) {
        inf >> mer;
        //std::cout << mer;
        if (prev > mer) {
            std::cerr << "Input error at " << inf.tellg() << " prev > current\n" << "Prev: "<< prev << "Cur: " << mer;
            std::cerr << "File kmer = " << fileKmers << std::endl;
            return 2;
        }
        prev = mer;
        fileKmers++;
    }
    std::cout << mer;
    if (totalKmers!=fileKmers) {
        std::cerr << "Total kmers (" << totalKmers << ") is different from read kmers " << fileKmers << std::endl;
        return 1;
    }

    std::cout << totalKmers << " = " << fileKmers << std::endl;

    return 0;
}