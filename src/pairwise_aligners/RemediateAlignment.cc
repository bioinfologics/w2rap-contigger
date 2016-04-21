// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 


// RemediateAlignment: given two reads (and quality scores), plus an alignment 
// between them (whose score has already been computed), try to replace the 
// alignment by a better one.

#include "Basevector.h"
#include "math/Functions.h"
#include "PackAlign.h"
#include "Qualvector.h"
#include "pairwise_aligners/RemediateAlignment.h"
#include "ScoreAlignment.h"
#include "pairwise_aligners/SmithWatBanded.h"
#include "pairwise_aligners/SmithWatBandedA.h"
#include "system/System.h"
#include "system/Types.h"
#include "Vec.h"

Bool RemediateAlignment( const basevector& rd1, const qualvector& q1,
     const basevector& rd2, const qualvector& q2, align& a, int& errors, 
     float& score )
{    if ( a.Nblocks( ) == 1 && errors < 10 ) return false;
     int pos1 = a.pos1( ), pos2 = a.pos2( );
     int bandwidth = Bandwidth(a);
     float actual = float( ActualErrors( rd1, rd2, a, 2, 3 ) ) / 2.0;
     float revised = SmithWatBanded(rd1, rd2, pos1 - pos2, bandwidth);
     if ( revised > actual ) PRINT2( revised, actual ); // would be bug!
     if ( revised >= actual ) return false;
     SmithWatBandedA( rd1, rd2, pos1 - pos2, bandwidth, a, errors );
     Regap( a, rd1, q1, rd2, q2 );
     score = ScoreAlignment( a, rd1, q1, rd2, q2 );
     return true;    }
