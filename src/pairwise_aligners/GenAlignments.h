// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
//


#ifndef GENALIGNMENTS_H
#define GENALIGNMENTS_H

#include "Basevector.h"
#include "pairwise_aligners/MakeAlignsMethod.h"
#include "pairwise_aligners/MutmerGraph.h"
#include "dvString.h"

template<int I, int k, int BLOCKS_PER_NODE> void
GenAlignments(
    const vecbasevector& EE, mutmer_graph<I, BLOCKS_PER_NODE>& M, int N0,
    String& aligns_file, ostream* log, int min_mutmer, int max_alignments,
    makealigns_method *method, int offset_A, int offset_B );
#endif
