#ifndef EFASTA_TOOLS_H
#define EFASTA_TOOLS_H

#include "Basevector.h"
#include "Bitvector.h"
#include "CoreTools.h"
#include "Fastavector.h"
#include "feudal/MasterVec.h"

struct Ambiguity
{
  size_t start;
  int size;
  String replace;
};

// Class efasta represents a single efasta record, exclusive of the header line.

class efasta;
typedef MasterVec<efasta> VecEFasta;

class efasta : public String {

     public:

     efasta( ) { }
     efasta( String::Allocator const& alloc ) : String(alloc) {}
     efasta( const String& s ) : String(s) { }
     efasta( const fastavector& v );  
     efasta( const basevector& v );
     efasta( const vec<basevector>& x );
     // simular to efasta( const vec<basevector>& x ), but zip
     // two sequences by more costly pairwise alignment
     efasta( const basevector& seq1, const basevector& seq2, const String& marking=String("") );
     efasta( const VecEFasta& x );

     // Length1 returns the length in bases of the record, assuming that the
     // first choice is made for all brackets.

     int Length1( Bool count_Ns = false ) const;

     int AmbEventCount( ) const;


     void FlattenTo( basevector& b ) const;

     // ExpandTo.  Convert to a list of fastavectors.  This will find ambiguous
     // base codes.   If max_count is specified and the number of fastavectors would 
     // exceed it, leave v empty and return False.

     Bool ExpandTo( vec<basevector>& v, const int max_count = -1 ) const;

     // ExpandToPaths.  Convert to a graph (see MakeGraph) and associated paths
     // encoded as graph branches. 

     Bool ExpandTo( vec< vec<basevector> >& G, vec< vec<int> >& paths, 
          const int max_count = -1 ) const;

     // Compute reverse-complement.

     void ReverseComplement( );

     // Compute blocks.

     void GetBlocks( vec< vec<String> >& blocks) const;

     // Erase the segment [start,stop) where start and stop are determined by 
     // the first choice in brackets.

     void Erase1( const int start, const int stop );
     
     // MakeGraph.  Convert to a graph.  Note that in the basevector case, this 
     // will convert Ns to As.

     void MakeGraph( vec< vec<basevector> >& G ) const;
     void MakeGraph( vec< vec<fastavector> >& G ) const;

      // Prints in a fasta format: "><string_id>\n" followed by the full base
      // sequence stored in the basevector; breaks
      // long sequences nicely into 80-character lines

      void Print( std::ostream& out, const String& id = "" ) const;


      // Write a efasta-format file which contains the bases in the
      // input fasta scaffolded together (with gaps, etc.) as defined
      // by the scaffolds.  If supplied, rc defines which contigs are
      // reverse-complement.


      friend void swap( efasta& e1, efasta& e2 )
      { e1.swap(e2); }
};

template<> struct Serializability<efasta>
{ typedef SelfSerializable type; };

// ExpandAmbCode.  Expand ambiguous base codes.  For example, M is converted to
// {A,C}, and N is converted to {A,C,G,T}.  This also converts n to N.

String ExpandAmbCode( char x );

// ValidateEfastaRecord.  Test a group of lines to see if they represent a valid 
// efasta record, exclusive of the header line.  This doesn't check for duplicates 
// in choose expressions.

void ValidateEfastaRecord( const vec<String>& lines, const String& msg = "",
     const Bool allow_empty_record = False );

// AllPlus.  Return the concatenation of some lines, including newlines.

String AllPlus( const vec<String>& lines );

#endif
