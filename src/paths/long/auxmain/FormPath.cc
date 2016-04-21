///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// FormPath.  Form a sequence of edges into a fasta record.  This has probably
// been written before.

#include "Basevector.h"
#include "ParseSet.h"
#include "TokenizeString.h"
#include "paths/long/SupportedHyperBasevector.h"

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String_Doc(IN, "hbv or shbv file for assembly");
     CommandArgument_String_Doc(SEQ, "comma-separated sequence of edges; multiple "
          "sequences may be specified using semicolon delimiters");
     EndCommandArguments;

     HyperBasevector hb;
     HyperBasevectorX hbx;
     bool use_hbx = false;
     if ( IN.Contains( ".shbv", -1 ) )
     {    SupportedHyperBasevector shb;
          BinaryReader::readFile( IN, &shb );
          hb = shb;    }
     else if ( IN.Contains( ".hbv", -1 ) )
     {    BinaryReader::readFile( IN, &hb );    }
     else if ( IN.Contains( ".hbx", -1 ) )
     {    use_hbx = true;
	  BinaryReader::readFile( IN, &hbx );    }
     else
     {    cout << "Illegal suffix for IN." << std::endl;
          Scram(1);    }

     vec<int> to_left, to_right;
     if (use_hbx)  {  // hbx
	 to_left = hbx.to_left();
	 to_right = hbx.to_right();
     }  else  // hbv or shbv
	 hb.ToLeft(to_left), hb.ToRight(to_right);

     vec<String> seq;
     Tokenize( SEQ, ';', seq );
     for ( int i = 0; i < seq.isize( ); i++ )
     {    vec<int> x;
          ParseIntSet( "{" + seq[i] + "}", x, false );
          for ( int j = 1; j < x.isize( ); j++ )
          {    if ( to_left[ x[j] ] != to_right[ x[j-1] ] )
               {    cout << "Warning: edges " << x[j-1] << " and " << x[j] 
                         << " are not adjacent.  Doesn't make sense." 
                         << std::endl;    }    }
          basevector b = (use_hbx ? hbx.Cat(x) :  hb.Cat(x)) ;
          b.Print( cout, seq[i] );    }    }
