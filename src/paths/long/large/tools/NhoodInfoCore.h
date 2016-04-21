///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef NHOOD_INFO_CORE_H
#define NHOOD_INFO_CORE_H

#include "CoreTools.h"
#include "TokenizeString.h"
#include "lookup/LookAlign.h"
#include "paths/RemodelGapTools.h"
#include "paths/long/CleanEfasta.h"
#include "paths/long/DisplayTools.h"
#include "paths/long/large/Lines.h"
#include "paths/long/large/tools/NhoodInfoState.h"
#include "paths/long/large/tools/NhoodInfoStuff.h"


void NhoodInfoCore( int argc, char *argv[] );

inline void NhoodInfoCore( const String& args )
{    vec<String> tokens;
     Tokenize( "NhoodInfo " + args, tokens );
     int argc = tokens.size( );
     vec<char*> argv(argc);
     for ( int i = 0; i < argc; i++ )
          argv[i] = const_cast<char*>(tokens[i].c_str());

     NhoodInfoCore( argc, &argv[0] );    }

#endif
