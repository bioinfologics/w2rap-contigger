//
// Created by Luis Yanes (EI) on 2019-04-29.
//

#ifndef W2RAP_CONTIGGER_SPECTRACN_HPP
#define W2RAP_CONTIGGER_SPECTRACN_HPP


#include <paths/HyperBasevector.h>
#include <paths/long/ReadPath.h>
#include <paths/long/BuildReadQGraph.h>
#include <util/kseq.hpp>
#include "Vec.h"

/**
 * @brief
 * Utility to calculate spectra-cn from a fasta or an hbv, works with gzipped files.
 */
class SpectraCN {
    static uint8_t findCountBin(std::ifstream &file, KMerNodeFreq_s &target, size_t lo, size_t hi);
    static uint64_t get_total_fasta_length(const std::string &fasta_path) {
        int l1(0);
        kseq seq1;
        int c1(0);
        int finished_early(0);
        uint64_t total_good_length(0);
        gzFile fp1 = gzopen(fasta_path.c_str(), "r");
        FunctorZlib gzr1;
        kstream<gzFile, FunctorZlib> ks1(fp1, gzr1);
        l1 = ks1.read(seq1);
        while (l1 >= 0) {
            if (seq1.seq.empty()) {
                std::cout << "Error " << std::string(fasta_path.c_str()) << " on read " << c1 << " is invalid"
                          << std::endl;
                finished_early = 1;
                break;
            }

            // Process this sequence!
            total_good_length += seq1.seq.size();


            c1++;
            l1 = ks1.read(seq1);
        }

        if (!finished_early) {
            if (l1 == -2) {
                std::cout << std::string(fasta_path.c_str()) << " is invalid at " << c1 << std::endl;
            }
        }
        gzclose(fp1);
        return total_good_length;
    }
    static void generateAssemblyKmers(const HyperBasevector &hb, const vec<int> &inv, std::shared_ptr<KmerList> &graph_kmer_freqs);

    static void
    writeHistogram(const String &dir, const String &name, const std::map<std::vector<uint64_t>, uint64_t> &totals_by_freq);

    static std::map<std::vector< uint64_t>, uint64_t>
    calculateHistogram(const std::string &dir, const std::shared_ptr<KmerList> &graph_kmer_freqs);

    static std::shared_ptr<KmerList> &
    generateFastaKmers(const std::string &fasta_path, std::shared_ptr<KmerList> &local_kmer_list);
public:
    static void DumpSpectraCN(const HyperBasevector &hb, const vec<int> &inv, const String &dir, const String &name);
    static void DumpSpectraCN(const std::string &fasta_path, const std::string &dir, const std::string &name);

};


#endif //W2RAP_CONTIGGER_SPECTRACN_HPP
