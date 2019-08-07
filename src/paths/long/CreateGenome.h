#ifndef CREATE_GENOME_H
#define CREATE_GENOME_H

#include "Basevector.h"
#include "CoreTools.h"
#include "IntPairVec.h"
#include "paths/HyperBasevector.h"
#include "paths/long/Logging.h"
#include "paths/long/DiscovarTools.h"

class ref_data {

     public:

     ref_data( ) : LG(12) { }

     vecbasevector G, G3, G3plus;
     vec<int> Gplus_ext;
     vec<HyperBasevector> GH;
     vec<double> ploidy;
     vec<bool> is_circular;
     VecIntPairVec Glocs, G3locs, G3pluslocs;
     int LG;
};

#endif
