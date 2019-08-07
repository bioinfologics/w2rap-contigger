#ifndef SCOREALIGNMENTONE
#define SCOREALIGNMENTONE

#include "Alignment.h"
//#include "math/Arith.h"
#include "Basevector.h"
#include "Qualvector.h"

float ScoreAlignment( const align& a, const basevector& rd1, 
     const QualVec& scores1, const basevector& rd2,
     const QualVec& scores2 = QualVec(0), int start1 = 0,
     int stop1 = -1, int start2 = 0, int stop2 = -1, Bool ignore_gaps = False );

float ScoreAlignment( Bool rd2_is_rc, const align& a, const basevector& rd1, 
     const QualVec& scores1, const basevector& rd2,
     const QualVec& scores2 = QualVec(0), int start1 = 0,
     int stop1 = -1, int start2 = 0, int stop2 = -1, Bool ignore_gaps = False );

int ScoreAlignmentPoly( const align& a, const basevector& rd1, 
     const QualVec& scores1, const basevector& rd2,
     const QualVec& scores2 = QualVec(0), int start1 = 0,
     int stop1 = -1, int start2 = 0, int stop2 = -1 );

int ScoreAlignmentPoly( Bool rd2_is_rc, const align& a, const basevector& rd1, 
     const QualVec& scores1, const basevector& rd2,
     const QualVec& scores2 = QualVec(0), int start1 = 0,
     int stop1 = -1, int start2 = 0, int stop2 = -1 );

void Regap( align& a, 
	    const basevector& rd1, const QualVec& scores1,
	    const basevector& rd2, const QualVec& scores2 );

void Regap( Bool rd2_is_rc, align& a, 
	    const basevector& rd1, const QualVec& scores1,
	    const basevector& rd2, const QualVec& scores2 );

inline 
void Regap( alignment& a, 
	    const basevector& rd1, const QualVec& scores1,
	    const basevector& rd2, const QualVec& scores2 )
{    
  align al = align(a);
  Regap( al, rd1, scores1, rd2, scores2 );
  a.Set( packalign(al), a.Errors( ) );    
}

#endif
