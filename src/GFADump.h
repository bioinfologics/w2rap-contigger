//
// Created by Bernardo Clavijo (TGAC) on 07/07/2016.
//

#ifndef W2RAP_CONTIGGER_GFADUMP_H
#define W2RAP_CONTIGGER_GFADUMP_H

void GFADump (std::string filename, const HyperBasevector &hb, const vec<int> &inv, const
ReadPathVec &paths, const int MAX_CELL_PATHS, const int MAX_DEPTH, bool find_lines);

#endif //W2RAP_CONTIGGER_GFADUMP_H
