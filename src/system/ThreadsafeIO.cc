/*
 * \file ThreadsafeIO.cc
 * \author tsharpe
 * \date May 14, 2009
 *
 * \brief Redefine cout, cerr in a threadsafe manner.
 */
#include "system/ThreadsafeIO.h"

SpinLockedData ThreadsafeStreambuf::gLock;

ThreadsafeStreambuf::int_type ThreadsafeStreambuf::overflow( int_type ch )
{
    if ( ch != traits_type::eof() )
    {
        *pptr() = ch;
        pbump(1);
    }
    return sync() ? traits_type::eof() : ch;
}

int ThreadsafeStreambuf::sync()
{
    SpinLocker locker(gLock);
    char* buf = pbase();
    mOS.write(buf,pptr() - buf);
    setp(buf,epptr());
    return mOS.fail(); // i.e., 0 if we wrote everything, 1 if we didn't.
}
