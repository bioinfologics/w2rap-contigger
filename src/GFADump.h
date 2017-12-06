//
// Created by Bernardo Clavijo (TGAC) on 07/07/2016.
//

#ifndef W2RAP_CONTIGGER_GFADUMP_H
#define W2RAP_CONTIGGER_GFADUMP_H

void GFADumpLines(std::string filename, const HyperBasevector &hb, const vec<int> &inv, const ReadPathVec &paths,
                  const int MAX_CELL_PATHS, const int MAX_DEPTH, bool find_lines,
                  const std::vector<int> &marked_edges = std::vector<int>());

// GFA Abyss equivalent, needs to get the KC (kmer coverage) from the raw_kmers.data file
// The kmer.data should probably be local to this function call
void GFADumpAbyss (std::string filename, const HyperBasevector &hb, const vec<int> &inv,
    const ReadPathVec &paths, const int MAX_CELL_PATHS, const int MAX_DEPTH, bool find_lines,
    const std::vector<int> & marked_edges = std::vector<int>());

void GFADumpRaw(std::string filename, const HyperBasevector &hb, const vec<int> &inv,
                const std::vector<int> &marked_edges = std::vector<int>(),
                const std::vector<int> &marked_vertices = std::vector<int>());
bool check_from_to(const HyperBasevector &hb);

#endif //W2RAP_CONTIGGER_GFADUMP_H
