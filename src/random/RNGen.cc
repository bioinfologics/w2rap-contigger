/*
 * \file RNGen.cc
 * \author tsharpe
 * \date Nov 17, 2011
 *
 * \brief
 */
#include "random/RNGen.h"

RNGen RNGen::gDflt(1u);
RNGen RNGen::gSystem;
SpinLockedData RNGen::gLock;
