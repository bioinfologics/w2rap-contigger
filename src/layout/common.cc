// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 


#include <iostream>
#include <fstream>
#include "layout/common.h"

ofstream cnull; 
ofstream cloga;
ofstream clogs;
ofstream permanent_file;

void coutDot() {
  cout << ".";
  flush(cout);
}

void coutDot(int inti, int intmod) {
  if (!intmod || !(inti%intmod)) coutDot();
}

void clogaDot() {
  cloga << ".";
  flush(cloga);
}

void clogaDot(int inti, int intmod) {
  if (!intmod || !(inti%intmod)) clogaDot();
}





