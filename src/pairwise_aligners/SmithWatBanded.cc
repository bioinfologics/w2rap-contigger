// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 


// SmithWatBanded( S, T, offset, bandwidth )
//
// Let S and T be basevectors.  Return the best (lowest) score of an alignment of S 
// with T, relative to the following rules:
//
// (a) a mismatch scores +1.0;
// (b) each blank in a gap scores +1.5;
// (c) gaps on either end of either read do not count.

// Only search for alignments with the given offset, plus or minus the given
// bandwidth.  The meaning of "offset" is explained by:
//
// example: offset = 2
// S -----------------
// T   ---------------------

// example: offset = -2
// S   ---------------------
// T -----------------

#include "Basevector.h"
#include "math/Functions.h"
#include "ShortVector.h"
#include "pairwise_aligners/SmithWatBanded.h"
#include "Vec.h"

float SmithWatBanded( const basevector& S, const basevector& T, 
     int offset, int bandwidth )
{
     if ( offset > (int) S.size( ) || offset < - (int) T.size( ) )
     {    std::cout << "Warning: SmithWatBanded passed nonsense arguments:" << std::endl;
          PRINT2( S.size( ), T.size( ) );
          PRINT2( offset, bandwidth );    }

     int left = offset - bandwidth, right = offset + bandwidth;

     const unsigned int mismatch_penalty = 2;
     const unsigned int gap_penalty = 3;
     const float divider = 2.0;

     const unsigned int Infinity = 1000000000;

     unsigned int n = S.size( ), N = T.size( );

     static vec<char> s;
     s.resize(n);
     for ( unsigned int i = 0; i < n; i++ )
          s[i] = S[i];
     unsigned int best_score = Infinity;
     static vec<unsigned int> x;
     x.resize(n+1);
     int istart = 0, istop = 0;

     for ( unsigned int i = 0; i <= n; i++ )
     {    x[i] = 0;
          if ( !( left <= (int) i && (int) i <= right ) ) x[i] = Infinity;    }

     for ( unsigned int j = Max( 0, -right-1 ); j <= Min( N-1, n-left ); j++ )
     {    x[0] = 0;
          if ( !( left <= -(int) j && -(int) j <= right ) ) x[0] = Infinity;
          unsigned int lastx = 0;
          if ( !( left <= -(int) (j-1) && -(int) (j-1) <= right ) ) lastx = Infinity;
          char* sp = &s[0];

          istart = Max( 0, left + (int) j - 1 );
          istop = Min( (int) n - 1, right + (int) j + 1 );
          unsigned int* xp = &x[0] + istart;
          if ( istart > 0 )
          {    lastx = *xp;
               *xp = Infinity;    }

          #define SWCORE(J)                                                        \
               else if ( T[j] == J )                                               \
               {    if ( istart <= istop )                                         \
                    {    unsigned int a                                            \
                              = lastx + mismatch_penalty * (sp[istart] != J);      \
                         unsigned int b = *xp + gap_penalty;                       \
                         ++xp;                                                     \
                         unsigned int c = *xp + gap_penalty;                       \
                         if ( !( left <= (int) (istart-1) - (int) (j-1)            \
                              && (int) (istart-1) - (int) (j-1) <= right ) )       \
                              a = Infinity;                                        \
                         if ( !( (int) istart - (int) (j-1) <= right ) )           \
                              b = Infinity;                                        \
                         if ( !( left <= (int) (istart-1) - (int) j ) )            \
                              c = Infinity;                                        \
                         lastx = *xp;                                              \
                         *xp = Min( Min( a, b ), c );    }                         \
                    if ( istart + 1 <= istop )                                     \
                    {    unsigned int a                                            \
                              = lastx + mismatch_penalty * (sp[istart+1] != J);    \
                         unsigned int b = *xp + gap_penalty;                       \
                         ++xp;                                                     \
                         unsigned int c = *xp + gap_penalty;                       \
                         if ( !( left <= (int) istart - (int) (j-1)                \
                              && (int) istart - (int) (j-1) <= right ) )           \
                              a = Infinity;                                        \
                         if ( !( (int) istart+1 - (int) (j-1) <= right ) )         \
                              b = Infinity;                                        \
                         if ( !( left <= (int) istart - (int) j ) )                \
                              c = Infinity;                                        \
                         lastx = *xp;                                              \
                         *xp = Min( Min( a, b ), c );    }                         \
                    for ( int i = istart + 2; i <= istop - 2; i++ )                \
                    {    unsigned int a = lastx + mismatch_penalty * (sp[i] != J); \
                         unsigned int b = *xp + gap_penalty;                       \
                         ++xp;                                                     \
                         unsigned int c = *xp + gap_penalty;                       \
                         lastx = *xp;                                              \
                         *xp = Min( Min( a, b ), c );    }                         \
                    if ( istop - 1 >= istart + 2 )                                 \
                    {    unsigned int a                                            \
                              = lastx + mismatch_penalty * (sp[istop-1] != J);     \
                         unsigned int b = *xp + gap_penalty;                       \
                         ++xp;                                                     \
                         unsigned int c = *xp + gap_penalty;                       \
                         if ( !( istop - (int) j <= right ) ) b = Infinity;        \
                         lastx = *xp;                                              \
                         *xp = Min( Min( a, b ), c );    }                         \
                    if ( istop >= istart + 2 )                                     \
                    {    unsigned int a                                            \
                              = lastx + mismatch_penalty * (sp[istop] != J);       \
                         unsigned int b = *xp + gap_penalty;                       \
                         ++xp;                                                     \
                         unsigned int c = *xp + gap_penalty;                       \
                         if ( !( istop - (int) j <= right ) ) a = Infinity;        \
                         if ( !( istop - (int) (j-1) <= right ) ) b = Infinity;    \
                         lastx = *xp;                                              \
                         *xp = Min( Min( a, b ), c );    }    }

          if ( 0 == 1 );
          SWCORE(0)
	  SWCORE(1)
	  SWCORE(2)
	  SWCORE(3)

          if ( istop < (int) n - 1 )
          {    ++xp;
               *xp = Infinity;    }

          if ( istop == (int) n - 1 ) best_score = Min( best_score, x[n] );    }

     for ( int i = istart; i <= istop; i++ )
          best_score = Min( x[i], best_score );
     return float(best_score)/divider;    }
