// Discovar messages are tagged "DISCOVAR MESSAGE".

#include <omp.h>
#include "paths/long/DiscovarTools.h"
#include "system/System.h"

namespace DiscovarTools {

    void ExitAssemblyEmpty() {    // DISCOVAR MESSAGE
        std::cout << "\nDear user, we are sad to report that your assembly is now "
                  << "empty.\nThis could be because your coverage is too low, or "
                  << "that something\nelse is wrong with your data.  Assembly of very "
                  << "small regions\n(e.g. ~1 kb) can also cause empty assemblies.\n"
                  << "Discovar complete.  Have a nice day.\n\n";
        Scram();
    }

    void ExitPathsEmpty() {    // DISCOVAR MESSAGE
        std::cout << "\nDear user, we are sad to report the current assembly is not supported by reads.\n"
                  << "This could happen if the reads do not align well with the reference (if supplied), or "
                  << "that something\nelse is wrong with your data.\n"
                  << "Discovar complete.  Have a nice day.\n\n";
        Scram();
    }


    void ExitNoCorrectedReads() {    // DISCOVAR MESSAGE
        std::cout << "\nDear user, we are sad to report that after error correction, "
                  << "no reads remain.\n"
                  << "Further assembly is thus impossible.\n"
                  << "Discovar complete.  Have a nice day.\n\n";
        Scram();
    }


    void ExitShortReads(const String &additional_info) {
        // DISCOVAR MESSAGE
        std::cout << "\nDiscovar has found that all reads in your input are too short";
        if (additional_info != "")
            std::cout << ":\n" << additional_info << "\n";
        std::cout << "\nDiscovar complete.  Have a nice day.\n\n";
        Scram();
    }
}

void SetThreads(uint &NUM_THREADS, bool malloc_per_thread_check) {
    if (malloc_per_thread_check) {
        static bool malloc_warning = false;
        String mtt = Getenv("MALLOC_PER_THREAD", "undefined");
        if (mtt != "1" && !malloc_warning) {
            // DISCOVAR MALLOC WARNING
            std::cout << Date() << ": Warning: recommend doing "
                      << "'setenv MALLOC_PER_THREAD 1'\n"
                      << Date() << ": before Discovar, to improve "
                      << "computational performance." << std::endl;
            malloc_warning = true;
        }
    }
    NUM_THREADS = configNumThreads(NUM_THREADS);
    omp_set_num_threads(NUM_THREADS);
}
