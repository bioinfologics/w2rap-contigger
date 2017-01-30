///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "Qualvector.h"

void Print( std::ostream &out, const QualVec &q, const String &name,
            const int scores_per_line )
{
    out << '>' << name;
    for ( QualVec::size_type i = 0; i < q.size(); ++i )
    {
        if (i % scores_per_line)
        {
            out << ' ';
        }
        else
        {
            out << '\n';
        }
        out << static_cast<unsigned int>(q[i]);
    }
    out << '\n';
}

std::pair <String, String> Stacked( const QualVec& quals) {
  uint read_length = quals.size();
  String line1(read_length, '9'), line2(read_length, '9'); // Max value displayed is 99, no real quality score should exceed this.
  for (uint i = 0; i < read_length; i++ ) {
    String qual_str = ToString(static_cast<unsigned int>(quals[i]));
    if (qual_str.size() == 1){
      line1[i] = ' ';
      line2[i] = qual_str[0];
    } else if (qual_str.size() == 2) {
      line1[i] = qual_str[0];
      line2[i] = qual_str[1];
    } // values > 99 are set to 99.
  }
  return std::make_pair(line1, line2);
}


void PrintStacked( std::ostream &out ,const QualVec& quals) {
  std::pair <String, String> qualities = Stacked(quals);
  out << qualities.first << std::endl << qualities.second << std::endl;
}

