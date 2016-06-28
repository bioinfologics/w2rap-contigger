// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
//


#include "Kclock.h"
#include <sys/resource.h>
#include <sys/time.h>
#include "system/Assert.h"

kclock::kclock() {

    struct timeval tp;
    struct timezone tzp;
    struct rusage r_usage ;

    int ftime_ret = gettimeofday( &tp, &tzp ) ;
    Assert( ftime_ret == 0 ) ;

    wtime_start = (double) tp.tv_sec + (double) tp.tv_usec / 1000000.0 ;

    int getrusage_ret = getrusage( RUSAGE_SELF, &r_usage ) ;
    Assert( getrusage_ret == 0 ) ;
    stime_start = (double) r_usage.ru_stime.tv_sec +
                  (double) r_usage.ru_stime.tv_usec / 1000000.0 ;
    utime_start = (double) r_usage.ru_utime.tv_sec +
                  (double) r_usage.ru_utime.tv_usec / 1000000.0 ;


#if 0
    cout << " gettimeof day returned: " <<  tp.tv_sec <<
         " and " <<  tp.tv_usec << std::endl ;
    ftime_ret = gettimeofday( &tp, &tzp ) ;
    cout << " gettimeof day returned: " <<  tp.tv_sec <<
         " and " <<  tp.tv_usec << std::endl ;

    cout << " wtime_start = " << wtime_start << std::endl ;
    cout << " utime_start = " << utime_start << std::endl ;
    cout << " stime_start = " << stime_start << std::endl ;

#endif
}

double kclock::WallClock() {
    struct timeval tp;
    struct timezone tzp;
    int ftime_ret = gettimeofday( &tp, &tzp ) ;
    Assert( ftime_ret == 0 ) ;

    double wtime_now = (double) tp.tv_sec + (double) tp.tv_usec / 1000000.0 ;

    return wtime_now - wtime_start ;
}

double kclock::SystemTime() {
    struct rusage r_usage ;
    int getrusage_ret = getrusage( RUSAGE_SELF, &r_usage ) ;
    Assert( getrusage_ret == 0 ) ;

    double stime_now = (double) r_usage.ru_stime.tv_sec + (double) r_usage.ru_stime.tv_usec / 1000000.0 ;

    return stime_now - stime_start ;
}

double kclock::UserTime() {
    struct rusage r_usage ;
    int getrusage_ret = getrusage( RUSAGE_SELF, &r_usage ) ;
    Assert( getrusage_ret == 0 ) ;

    double utime_now = (double) r_usage.ru_utime.tv_sec + (double) r_usage.ru_utime.tv_usec / 1000000.0 ;

    return utime_now - utime_start ;
}

ostream& operator<<(ostream &o, kclock &c) {
    struct timeval tp;
    struct timezone tzp;
    struct rusage r_usage ;
    double wtime_now ;
    double stime_now ;
    double utime_now ;

    int ftime_ret = gettimeofday( &tp, &tzp ) ;
    Assert( ftime_ret == 0 ) ;

    wtime_now = (double) tp.tv_sec + (double) tp.tv_usec / 1000000.0 ;

    int getrusage_ret = getrusage( RUSAGE_SELF, &r_usage ) ;
    Assert( getrusage_ret == 0 ) ;

    stime_now = (double) r_usage.ru_stime.tv_sec + (double) r_usage.ru_stime.tv_usec / 1000000.0 ;
    utime_now = (double) r_usage.ru_utime.tv_sec + (double) r_usage.ru_utime.tv_usec / 1000000.0 ;

#ifdef SHOW_MICROSECONDS
    o << "## WALLCLOCK= " <<  ( wtime_now - c.wtime_start )
      << " SYS= " <<   ( stime_now - c.stime_start )
      << " USER= " <<  ( utime_now - c.utime_start )
      << " ##" ;
#else
    o << "## WALLCLOCK= " << (int) ( wtime_now - c.wtime_start )
      << " SYS= " <<  (int) ( stime_now - c.stime_start )
      << " USER= " << (int) ( utime_now - c.utime_start )
      << " ##" ;
#endif
    return o ;
}
