#ifndef EXTRACT_READS_H
#define EXTRACT_READS_H

#include "Basevector.h"
#include "CoreTools.h"
#include "feudal/PQVec.h"

void ExtractPairedReads( std::string r1file, std::string r2file, vecbvec & pReads, VecPQVec & quals );

#endif
