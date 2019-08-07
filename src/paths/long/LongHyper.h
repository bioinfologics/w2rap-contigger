// Make the corrected reads into a HyperKmerPath assembly.

#ifndef LONG_HYPER_H
#define LONG_HYPER_H

#include "CoreTools.h"
#include "efasta/EfastaTools.h"
#include "paths/long/LongProtoTools.h"
#include "paths/long/PairInfo.h"
#include "paths/long/SupportedHyperBasevector.h"

Bool LongHyper(const VecEFasta &correctede, const vec<pairing_info> &cpartner, SupportedHyperBasevector &shb,
               const long_heuristics &heur, const long_logging_control &log_control, const long_logging &logc,
               bool useOldLRPMethod);

#endif
