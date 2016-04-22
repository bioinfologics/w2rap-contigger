// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 


#ifndef COMMON
#define COMMON

#include <iosfwd>
#include <fstream>


//
// The following are our standard output streams
//
extern std::ofstream cnull; 
extern std::ofstream cloga;
extern std::ofstream clogs;
extern std::ofstream permanent_file;


void coutDot();                     /* Outputs a dot in std::cout */

void coutDot(int inti, int intmod); /* Outputs a dot in std::cout if
				       inti%intmod == 0 */

void clogaDot();                    /* Outputs a dot in cloga */

void clogaDot(int, int);            /* Outputs a dot in cloga if
				       inti%intmod == 0 */

#endif
