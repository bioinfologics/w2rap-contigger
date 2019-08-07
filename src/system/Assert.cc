///\file Assert.cc
/// Assert.cc provides some of the guts of the Assert macros
#include "system/Assert.h"
#include "system/Exit.h"
#include <iostream>

void Assert::reportVals( char const* loc, char const* func, char const* vals )
{
    std::cout << loc << " failed in function\n" << func << "\n";

    if ( vals )
        std::cout << "with values " << vals;

    std::cout << std::endl;
}

void Assert::reportValsAndDie( char const* loc, char const* func, char const* vals )
{
    reportVals(loc,func,vals);
    CRD::exit(1);
}
