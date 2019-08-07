#ifndef GENERIC_H_
#define GENERIC_H_

#include "system/ThreadsafeIO.h"
#include <cstddef>
#include <fstream>

class UtilizationReporter
{
public:
    void report( void const* addr, size_t bits, size_t siz, size_t cap,
                    char const* type )
    { if ( mpOS ) writeLine(addr,bits,siz,cap,type); }

    static UtilizationReporter gInstance;

private:
    UtilizationReporter();
    ~UtilizationReporter();
    UtilizationReporter( UtilizationReporter const& )=delete;
    UtilizationReporter& operator=( UtilizationReporter const& )=delete;

    void writeLine( void const* addr, size_t bits, size_t siz, size_t cap,
                    char const* type );

    std::ofstream mLog;
    ThreadsafeOStream* mpOS;
};

#endif // GENERIC_H_
