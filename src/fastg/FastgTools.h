///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef FASTG_TOOLS_H
#define FASTG_TOOLS_H


#include "CoreTools.h"
#include "efasta/EfastaTools.h"

// class specifying all FASTG specifications
class fastg_meta {

  public:
    static size_t GetDefaultGapSize(); //the number of 'N' in a gap's canonical sequence
    static String GetVersion();
    static String GetFileHeader(const String assembly_name);
    static String GetFileFooter();

//   For FASTG<->EFASTA conversion, 3 modes are required due to legacy reasons
//   Mode 0 yields to correct value, Mode 1/2 are used to accommodate old code
//   MODE_0 == 0: variation in length according to actual FASTG spec
//   MODE_1 == 1: the variation of 'gaps' are ignore;
//               instead, the gap equals to max( average gap width, 1)
//   MODE_2 == 2: the variation of 'gaps' are ignore;
//               instead, the gap equals to max( average gap width,0)
//
// Future work should stick with MODE_0
//

    enum efasta_mode {MODE_0,MODE_1,MODE_2,N_MODE};

  private:
    //these are implmenetation dependent and are defined in the implementation file FastgTools.cc
    static const String sVersion_;
    static const size_t uDefaultGapSize_; //the number of 'N' in the canonical sequence
};


class recfastg;

// Class basefastg represents bases record in fastg string.

class basefastg : public String {

  public:

    basefastg( );
    basefastg( const String& s );
    basefastg( const efasta& ef );
    basefastg( const int& sep, const int& dev );
    basefastg( const superb& s, const VecEFasta& econtigs);
    basefastg( const superb& s, const vec<basefastg>& fcontigs);
    basefastg( const superb& s, const vec<recfastg>& fcontigs);

    void AsSections( vec<String>& sections ) const;

    int Length1(const Bool count_Ns = false) const;



//   Three modes are required due to legacy reasons
//   Mode 0 yields to correct value, Mode 1/2 are used to accommodate old code
//   fastg_meta::MODE_0 == 0: variation in length according to actual FASTG spec
//   fastg_meta::MODE_1 == 1: the variation of 'gaps' are ignore;
//               instead, the gap equals to max( average gap width, 1)
//   fastg_meta::MODE_2 == 2: the variation of 'gaps' are ignore;
//               instead, the gap equals to max( average gap width,0)
    int MinLength(const fastg_meta::efasta_mode uMode ) const;
    int MaxLength(const fastg_meta::efasta_mode uMode ) const;

    Bool IsGapless() const;
};


class headfastg : public String {

  public:

    headfastg( );
    headfastg( const String& id );
    headfastg( const String& id, const vec<String>& next_ids );
    headfastg( const String& id, const vec<int>& next_ids );

};

class recfastg {

  private:
    headfastg header_;
    basefastg bases_;

  public:
    recfastg();
    recfastg( const headfastg& header, const basefastg& bases);


    void Set( const headfastg& header, const basefastg& bases);

    int Length1() const;
//   Three modes are required due to legacy reasons
//   Mode 0 yields to correct value, Mode 1/2 are used to accommodate old code
//   fastg_meta::MODE_0 == 0: variation in length according to actual FASTG spec
//   fastg_meta::MODE_1 == 1: the variation of 'gaps' are ignore;
//               instead, the gap equals to max( average gap width, 1)
//   fastg_meta::MODE_2 == 2: the variation of 'gaps' are ignore;
//               instead, the gap equals to max( average gap width,0)
    int MinLength(const fastg_meta::efasta_mode uMode ) const;
    int MaxLength(const fastg_meta::efasta_mode uMode ) const;
    Bool IsGapless() const;

    Bool ReadRecord( std::ifstream& in );

    const basefastg& bases() const;
    const headfastg& header() const;

    // Prints in a following format: "><header_>;\n" followed by the full bases
    // sequence and ambiguity information; breaks
    // long sequences nicely into 80-character lines

    void Print( std::ostream& out ) const;

    void Print( std::ostream& out, const String& id ) const;

    void AsScaffold( superb& s, vec<recfastg>& fcontigs, int& lastTid ) const;

    void AsFasta( fastavector& fa ) const;
    void AsFastb( basevector& fb ) const;

    // convert to EFASTA, MODE_0==MODE_1, with gap=max(1,avg-gap), and MODE_2 for gap=max(0,avg-gap)
    void AsEfasta( efasta& fe , const fastg_meta::efasta_mode uMode) const;

};


void LoadFastg( const String& fn, vec<recfastg>& records );

void WriteFastg( const String& fn, const vec<recfastg>& records );

// take a vec of recfastg and convert and create FASTA and EFASTA files
void WriteFastaEfasta( const String& fasta_file, const String& efasta_file
                       , const vec<recfastg>& records, const Bool ncbi_format=false);




#endif
