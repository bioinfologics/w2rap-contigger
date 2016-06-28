// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
//


//  Id: MainArrays.h,v 1.26 2000/07/17 15:30:22 serafim Exp $
#ifndef MAIN_ARRAYS_H
#define MAIN_ARRAYS_H

#define cken cnull
#define cser cnull
#define cdav cnull
#ifdef KEN
#undef cken
#define cken cloga
#endif
#ifdef KEN
#define KENSPRINTVAL(NAME) cken << #NAME << " = " << NAME << std::endl ;
#else
#define KENSPRINTVAL(NAME)
#endif
#ifdef DAVE
#undef cdav
#define cdav cloga
#endif
#ifdef SERAFIM
#undef cser
#define cser cloga
#endif
#include <fstream>
#include <list>
#include <set>
#include <deque>
#include <algorithm>
#include <fstream>

#include "Kclock.h"
#include "math/Functions.h"


#include <ctime>
const int ios_binary =  ios::binary;


enum genome_mark { /* This code for marking the genome is used
		      only in the DisplayGenomeAssembly function which
		      provides a display of the layout, after ./Arachne
		      is run. That function is not used too much lately
		      (i.e. as of May 25, 2001). */
    SC_START               =   1,
    SC_STOP                =   2,
    FILLED_GAP             =   4,
    UNFILLED_GAP           =   8,
    CLEAN_SEQUENCE         =  16,
    FEW_MISPLACED          =  32, // < 1%
    SOME_MISPLACED         =  64, // 1% to 10%
    MANY_MISPLACED         = 128, // 10% to 50%
    MOST_MISPLACED         = 256, // 50% to 90%
    ALL_MISPLACED          = 512, // 90% to 100%
    COVERED_WITH_READ      = 1024,
    COVERED_WITH_TWO_READS = 2048
};

const int ReadNotSet = - 1;
extern int rcnt;      // Total number of reads;
extern int rcnt_main; /* Number of main reads, excluding auxiliary fake reads
			 that may be added to provide our knowledge of pairing
			 information say, such as in comparative assembly projects.
			 Only the main reads are counted when computing the
			 logodd ratio of a contig being unique. Also, only the
			 main reads give coordinates for determining the actual
			 location of a contig. The main reads are assumed to be
			 the first rcnt_main reads in the read array.
			 -- Serafim, July 12, 2001. */

class rawgraph;

// Sante --- Thu Oct 11 14:25:26 EDT 2001
extern longlong genome_length ;
// extern unsigned int genome_length ;

extern ofstream cnull ;
extern ofstream cloga;
extern ofstream clogs;
extern ofstream permanent_file;

#endif
