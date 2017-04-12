///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "ParallelVecUtilities.h"
#include "efasta/EfastaTools.h"
#include "kmers/KmerRecord.h"
#include "paths/long/CreateGenome.h"
#include "paths/long/ultra/GetFriendsAndAlignsInitial.h"
#include "paths/long/LargeKDispatcher.h"

void ReportPeakMem( const String msg )
{    std::cout << Date( ) << ": " << msg << ( msg != "" ? ", " : "" ) 
          << "peak mem = " << std::setiosflags(std::ios::fixed)
          << std::setprecision(1) << PeakMemUsageBytes( ) / 1000000000.0
          << std::resetiosflags(std::ios::fixed) << " GB" << std::endl;    }

int SelectK2(const VecEFasta &corrected, const double K2frac,
             const long_logging &logc, const long_heuristics &heur) {
    int K2 = -1;
    double hclock = WallClockTime();
    vec<int> lens;
    for (size_t i = 0; i < corrected.size(); i++) {
        if (corrected[i].size() == 0) continue;
        lens.push_back(corrected[i].Length1());
    }
    Sort(lens);
    if (lens.empty()) DiscovarTools::ExitNoCorrectedReads();
    int med = Median(lens);
    double target_length = K2frac * double(med);
    if (logc.MIN_LOGGING) {
        std::cout << Date() << ": " << lens.size()
                  << " corrected/closed pairs, having median length " << med
                  << std::endl;
    }
    double min_err = 1000000;
    // iterate over allowable K values
    for (auto itr = BigK::begin(), end = BigK::end(); itr != end; ++itr) {
        double err = Abs(target_length - *itr);
        if (err < min_err) {
            min_err = err;
            K2 = *itr;
        }
    }
    if (K2 == -1)
        FatalErr("Unable to identify a suitable value for K2.");
    if (logc.STATUS_LOGGING) std::cout << Date() << ": using K2 = " << K2 << std::endl;
    REPORT_TIME(hclock, "used choosing K2");
    return K2;
}



