///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef	SYSTEM_ASSERT_H
#define SYSTEM_ASSERT_H

#include <sstream>
#include <execinfo.h>
#include <cxxabi.h>

namespace Assert
{

    /** Print a demangled stack backtrace of the caller function to FILE* out. */
    static inline void print_stacktrace(FILE *out = stderr, unsigned int max_frames = 63)
    {
        fprintf(out, "stack trace:\n");

        // storage array for stack trace address data
        void* addrlist[max_frames+1];

        // retrieve current stack addresses
        int addrlen = backtrace(addrlist, sizeof(addrlist) / sizeof(void*));

        if (addrlen == 0) {
            fprintf(out, "  <empty, possibly corrupt>\n");
            return;
        }

        // resolve addresses into strings containing "filename(function+address)",
        // this array must be free()-ed
        char** symbollist = backtrace_symbols(addrlist, addrlen);

        // allocate string which will be filled with the demangled function name
        size_t funcnamesize = 256;
        char* funcname = (char*)malloc(funcnamesize);

        // iterate over the returned symbol lines. skip the first, it is the
        // address of this function.
        for (int i = 1; i < addrlen; i++)
        {
            char *begin_name = 0, *begin_offset = 0, *end_offset = 0;

            // find parentheses and +address offset surrounding the mangled name:
            // ./module(function+0x15c) [0x8048a6d]
            for (char *p = symbollist[i]; *p; ++p)
            {
                if (*p == '(')
                    begin_name = p;
                else if (*p == '+')
                    begin_offset = p;
                else if (*p == ')' && begin_offset) {
                    end_offset = p;
                    break;
                }
            }

            if (begin_name && begin_offset && end_offset
                && begin_name < begin_offset)
            {
                *begin_name++ = '\0';
                *begin_offset++ = '\0';
                *end_offset = '\0';

                // mangled name is now in [begin_name, begin_offset) and caller
                // offset in [begin_offset, end_offset). now apply
                // __cxa_demangle():

                int status;
                char* ret = abi::__cxa_demangle(begin_name,
                                                funcname, &funcnamesize, &status);
                if (status == 0) {
                    funcname = ret; // use possibly realloc()-ed string
                    fprintf(out, "  %s : %s+%s\n",
                            symbollist[i], funcname, begin_offset);
                }
                else {
                    // demangling failed. Output function name as a C function with
                    // no arguments.
                    fprintf(out, "  %s : %s()+%s\n",
                            symbollist[i], begin_name, begin_offset);
                }
            }
            else
            {
                // couldn't parse the line? print the whole line.
                fprintf(out, "  %s\n", symbollist[i]);
            }
        }

        free(funcname);
        free(symbollist);
    }

    void reportVals( char const* loc, char const* func, char const* vals );
void reportValsAndDie( char const* loc, char const* func, char const* vals )
            __attribute__((__noreturn__));

template <class T, class U>
void reportAndDie( T const& t, U const& u, char const* loc, char const* func )
{

    print_stacktrace();
    std::ostringstream oss;
    oss << "arg1 = " << t << " and arg2 = " << u;
  reportValsAndDie(loc,func,oss.str().c_str()); }

template <class T, class U>
inline void eq( T const& t, U const& u, char const* loc, char const* func )
{ if ( !(t == u) ) reportAndDie(t,u,loc,func); }

template <class T, class U>
inline void ne( T const& t, U const& u, char const* loc, char const* func )
{ if ( !(t != u) ) reportAndDie(t,u,loc,func); }

template <class T, class U>
inline void gt( T const& t, U const& u, char const* loc, char const* func )
{ if ( !(t > u) ) reportAndDie(t,u,loc,func); }

template <class T, class U>
inline void ge( T const& t, U const& u, char const* loc, char const* func )
{ if ( !(t >= u) ) reportAndDie(t,u,loc,func); }

template <class T, class U>
inline void lt( T const& t, U const& u, char const* loc, char const* func )
{ if ( !(t < u) ) reportAndDie(t,u,loc,func); }

template <class T, class U>
inline void le( T const& t, U const& u, char const* loc, char const* func )
{ if ( !(t <= u) ) reportAndDie(t,u,loc,func); }

inline void yes( bool t, char const* loc, char const* func )
{ if ( !t ) reportValsAndDie(loc,func,0); }

inline void no( bool t, char const* loc, char const* func )
{ if ( t ) reportValsAndDie(loc,func,0); }

template <class T, class U>
void report( T const& t, U const& u, char const* loc, char const* func )
{ std::ostringstream oss; oss << "arg1 = " << t << " and arg2 = " << u;
  reportVals(loc,func,oss.str().c_str()); }

template <class T, class U>
inline void eqTest( T const& t, U const& u, char const* loc, char const* func )
{ if ( !(t == u) ) report(t,u,loc,func); }

template <class T, class U>
inline void neTest( T const& t, U const& u, char const* loc, char const* func )
{ if ( !(t != u) ) report(t,u,loc,func); }

template <class T, class U>
inline void gtTest( T const& t, U const& u, char const* loc, char const* func )
{ if ( !(t > u) ) report(t,u,loc,func); }

template <class T, class U>
inline void geTest( T const& t, U const& u, char const* loc, char const* func )
{ if ( !(t >= u) ) report(t,u,loc,func); }

template <class T, class U>
inline void ltTest( T const& t, U const& u, char const* loc, char const* func )
{ if ( !(t < u) ) report(t,u,loc,func); }

template <class T, class U>
inline void leTest( T const& t, U const& u, char const* loc, char const* func )
{ if ( !(t <= u) ) report(t,u,loc,func); }

inline void yesTest( bool t, char const* loc, char const* func )
{ if ( !t ) reportVals(loc,func,0); }

inline void noTest( bool t, char const* loc, char const* func )
{ if ( t ) reportVals(loc,func,0); }

}

#define ASSERT_ARGS_HELP(x) #x
#define ASSERT_ARGS_HELP2(x) ASSERT_ARGS_HELP(x)
#define ASSERT_ARGS(macro,x)\
    #macro "(" #x ") at " __FILE__ ":" ASSERT_ARGS_HELP2(__LINE__),\
    __PRETTY_FUNCTION__
#define ASSERT_ARGS2(macro,x,y)\
    #macro "(" #x "," #y ") at " __FILE__ ":" ASSERT_ARGS_HELP2(__LINE__),\
    __PRETTY_FUNCTION__

// this set of macros halt execution if the condition is false
#define ForceAssertEq(x,y) Assert::eq(x,y,ASSERT_ARGS2(ForceAssertEq,x,y))
#define ForceAssertNe(x,y) Assert::ne(x,y,ASSERT_ARGS2(ForceAssertNe,x,y))
#define ForceAssertGt(x,y) Assert::gt(x,y,ASSERT_ARGS2(ForceAssertGt,x,y))
#define ForceAssertGe(x,y) Assert::ge(x,y,ASSERT_ARGS2(ForceAssertGe,x,y))
#define ForceAssertLt(x,y) Assert::lt(x,y,ASSERT_ARGS2(ForceAssertLt,x,y))
#define ForceAssertLe(x,y) Assert::le(x,y,ASSERT_ARGS2(ForceAssertLe,x,y))
#define ForceAssert(x)     Assert::yes(x,ASSERT_ARGS(ForceAssert,x))
#define ForceAssertNot(x)  Assert::no(x,ASSERT_ARGS(ForceAssertNot,x))

// this set of macros just print a message if the condition is false
#define TestAssertEq(x,y) Assert::eqTest(x,y,ASSERT_ARGS2(AssertEq,x,y))
#define TestAssertNe(x,y) Assert::neTest(x,y,ASSERT_ARGS2(AssertNe,x,y))
#define TestAssertGt(x,y) Assert::gtTest(x,y,ASSERT_ARGS2(AssertGt,x,y))
#define TestAssertGe(x,y) Assert::geTest(x,y,ASSERT_ARGS2(AssertGe,x,y))
#define TestAssertLt(x,y) Assert::ltTest(x,y,ASSERT_ARGS2(AssertLt,x,y))
#define TestAssertLe(x,y) Assert::leTest(x,y,ASSERT_ARGS2(AssertLe,x,y))
#define TestAssert(x)     Assert::yesTest(x,ASSERT_ARGS(Assert,x))
#define TestAssertNot(x)  Assert::noTest(x,ASSERT_ARGS(AssertNot,x))

// this set of macros does nothing if NDEBUG is defined
// they halt execution when the condition is false if NDEBUG isn't defined
#ifdef NDEBUG
#define AssertEq(x,y) ((void)(sizeof(x)+sizeof(y)))
#define AssertNe(x,y) ((void)(sizeof(x)+sizeof(y)))
#define AssertGt(x,y) ((void)(sizeof(x)+sizeof(y)))
#define AssertGe(x,y) ((void)(sizeof(x)+sizeof(y)))
#define AssertLt(x,y) ((void)(sizeof(x)+sizeof(y)))
#define AssertLe(x,y) ((void)(sizeof(x)+sizeof(y)))
#define Assert(x)     ((void)sizeof(x))
#define AssertNot(x)  ((void)sizeof(x))
#else
#define AssertEq(x,y) Assert::eq(x,y,ASSERT_ARGS2(AssertEq,x,y))
#define AssertNe(x,y) Assert::ne(x,y,ASSERT_ARGS2(AssertNe,x,y))
#define AssertGt(x,y) Assert::gt(x,y,ASSERT_ARGS2(AssertGt,x,y))
#define AssertGe(x,y) Assert::ge(x,y,ASSERT_ARGS2(AssertGe,x,y))
#define AssertLt(x,y) Assert::lt(x,y,ASSERT_ARGS2(AssertLt,x,y))
#define AssertLe(x,y) Assert::le(x,y,ASSERT_ARGS2(AssertLe,x,y))
#define Assert(x)     Assert::yes(x,ASSERT_ARGS(Assert,x))
#define AssertNot(x)  Assert::no(x,ASSERT_ARGS(AssertNot,x))
#endif

#endif // SYSTEM_ASSERT_H
