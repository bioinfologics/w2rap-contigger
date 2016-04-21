///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// ============================================================================
//
/// CmpSeq: compare two sets of DNA sequences, in fasta or fastb format.
/// \file CmpSeq.cc
/// Sequences are defined by FILE1 and FILE2.  Minimal calling syntax:
///
///      CmpSeq FILE1=file1 FILE2=file2
///
/// where the base for FILE1 and FILE2 is defined by the environment variable
/// ARACHNE_PRE.(You can set this on an individual run using PRE and/or BASE.)
///
/// See below:
///
/// 1.  Miscellaneous options.
/// 2.  Control of alignment acceptance.
/// 3.  Examples.
/// 4.  Many undocumented options, after "BeginCommandArguments".
///
///
/// Miscellaneous options:
///
/// OUTPUT=output-file sends all output to the given file; otherwise std::cout is used.
///
/// Note: Reads or contigs or whatever in FILE1 and FILE2 are broken at long
/// ( >= 20 base ) sequences of ambiguous bases, unless you set NOBREAK=True.
///
/// COUNT1 (if set) is the number of sequences from FILE1 to be used.
/// COUNT2 (if set) is the number of sequences from FILE2 to be used.
/// START2 (if set) is the first sequence in FILE2 to be used.
///
/// QUAL1 (if set) is the name of a quality score file, in qualb format,
/// matching FILE1.  Ditto for QUAL2.  Setting only QUAL2 is not allowed.
///
/// If OMIT_PROPER_ALIGNERS=True, and a sequences from FILE1 properly aligns
/// with a sequences from FILE2, it is completely ignored.
///
/// If MIN_COVERAGE1 > 0, then a sequences from FILE1 is ignored unless at
/// least the fraction MIN_COVERAGE1 (a number between 0 and 1) of its bases
/// are covered  by bases from FILE2 reads.
///
/// If NE=True, don't align sequences from FILE1 and FILE2 with the same indices.
/// If EQ=True, don't align sequences from FILE1 and FILE2 with different indices.
/// If LT=True, don't align sequences from FILE1 and FILE2 unless the index of the
/// first sequences is less than the index of the second.  Normally this option is
/// only useful if FILE1 = FILE2.
///
/// If VISUAL_COVERAGE=True, for each accepted proper alignment, mark the
/// corresponding positions on the FILE2 sequence, as "hit right" or "hit wrong".
/// Cumulate, keeping track of all positions.  Then output the FILE2 sequence
/// positions, as not hit at all, hit wrong, or hit right.  This gives a picture
/// of how diverged one set of sequences is from another, although it will
/// significantly underrepresent this divergence.
///
/// ==========================================================================
///
///
/// ===========================================================================
///
/// Control of alignment acceptance.
///
/// If TRIM_ALIGNMENT_ENDS is set, then alignments are cleaned up by removing
/// poorly aligning parts at the ends of the aligning region.  This influences
/// the affect of the following options.
///
/// Alignments for which the overlap is < MIN_OVERLAP (default 50) are ignored.
///
/// If MAX_ERROR_RATE is set, then alignments whose error rate exceeds
/// MAX_ERROR_RATE are ignored.
///
/// If a value for MAX_QUAL_SCORE is given (in which case QUAL1 must be set), then
/// quality scores (defined by QUAL1 and QUAL2 if it is set) are used to grade
/// alignments, and those with scores > MAX_QUAL_SCORE are ignored.
///
/// If MIN_PERFECT_MATCH is set, then alignments which do not subsume a perfect match
/// of length MIN_PERFECT_MATCH are ignored.
///
/// ==========================================================================
///
///
/// =======================================================================
///
/// Examples.
///
/// 1. Compare all the trimmed reads in a given project with each other, using
/// quality scores to evaluate their alignments. (This assumes the existence of
/// a small project "projects/L3191", with run directory "reads",
/// that you have done
/// FetchAndTrim in, e.g. creatable by
/// "Assemble DATA=projects/L3191 RUN=reads STOP=Fetch" ).
///
///      CmpSeq BASE=projects/L3191/reads
///             FILE1=reads.fastb QUAL1=reads.qualb
///             FILE2=reads.fastb QUAL2=reads.qualb
///             LT=True
///
/// 2. Similar, but use final contigs instead of trimmed reads.
/// (Remove STOP=Fetch
/// if you want to generate an appropriate data set via the above Assemble.)
///
///      CmpSeq BASE=projects/L3191/reads
///             FILE1=mergedcontigs.fasta QUAL1=mergedcontigs.qualb
///             FILE2=mergedcontigs.fasta QUAL2=mergedcontigs.qualb
///             LT=True
///
/// 3. And now, trimmed reads vs final contigs:
///
///      CmpSeq BASE=projects/L3191/reads
///             FILE1=reads.fastb QUAL1=reads.qualb
///             FILE2=mergedcontigs.fasta QUAL2=mergedcontigs.qualb
///
/// 4. Trimmed reads from one project versus final contigs on another, showing
/// only those reads (none) which align to high stringency:
///
///      CmpSeq BASE=projects
///             FILE1=L653/reads/reads.fastb QUAL1=L653/reads/reads.qualb
///             FILE2=L3191/reads/mergedcontigs.fasta
///             QUAL2=L3191/reads/mergedcontigs.qualb
///             MAX_QUAL_SCORE=100
///
/// =======================================================================

#include "Alignment.h"
#include "math/Arith.h"
#include "Badness.h"
#include "Basevector.h"
#include "kmers/KmerShape.h"
#include "MainTools.h"
#include "FetchReads.h"
#include "math/Functions.h"
#include "math/HoInterval.h"
#include "pairwise_aligners/MakeAligns.h"
#include "NQS.h"
#include "Overlap.h"
#include "pairwise_aligners/PerfectAlignment.h"
#include "PrintAlignment.h"
#include "Quality.h"
#include "Qualvector.h"
#include "pairwise_aligners/RemediateAlignment.h"
#include "Rmr.h"
#include "ScoreAlignment.h"
#include "Set.h"
#include "ShortVector.h"
#include "system/ErrNo.h"
#include "pairwise_aligners/SmithWatBandedA.h"
#include "TrimAlignmentEnds.h"
#include "VecString.h"

const unsigned int undefined = 1000000000;

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String_OrDefault(BASE, "");

     // Data for first sequences:

     CommandArgument_String_OrDefault(FILE1, "");
     CommandArgument_Bool_OrDefault(FILE1_FASTB, False);
     CommandArgument_UnsignedInt_OrDefault(START1, 0);
     CommandArgument_UnsignedInt_OrDefault(COUNT1, 0);
     CommandArgument_String_OrDefault(NAME1, "");
     CommandArgument_String_OrDefault(QUAL1, "");
     CommandArgument_String_OrDefault(HEAD1, "");

     // Data for second sequences:

     CommandArgument_String_OrDefault(FILE2, "");
     CommandArgument_Bool_OrDefault(FILE2_FASTB, False);
     CommandArgument_UnsignedInt_OrDefault(START2, 0);
     CommandArgument_UnsignedInt_OrDefault(COUNT2, 0);
     CommandArgument_String_OrDefault(NAME2, "");
     CommandArgument_String_OrDefault(QUAL2, "");
     CommandArgument_String_OrDefault(HEAD2, "");

     // Other options:

     CommandArgument_String_OrDefault(OUTPUT, "");
     CommandArgument_String_OrDefault(OMIT_IF_PROPER_OVERLAP_EXISTS, "False");
     CommandArgument_UnsignedInt_OrDefault(PASSES, 100);
     CommandArgument_UnsignedInt_OrDefault(ACTUAL_PASSES, 100);
     CommandArgument_UnsignedInt_OrDefault(K, 24);
     CommandArgument_UnsignedInt_OrDefault(MIN_OVERLAP, 50);
     CommandArgument_Double_OrDefault(MIN_COVERAGE1, 0.0);
     CommandArgument_Bool_OrDefault(USE_ONLY_SEMIPROPER_ALIGNMENTS, False);
     CommandArgument_Bool_OrDefault(OMIT_PROPER_ALIGNERS, False);
     CommandArgument_UnsignedInt_OrDefault(MIN_SIZE1, 0);
     CommandArgument_Double_OrDefault(MAX_ERROR_RATE, 1.0);
     CommandArgument_UnsignedInt_OrDefault(MAXCLIQ, 1000);
     CommandArgument_UnsignedInt_OrDefault(MAX_ALIGNS, 10000);
     CommandArgument_UnsignedInt_OrDefault(MIN_MUTMER, 0);
     CommandArgument_UnsignedInt_OrDefault(BLOCK_SIZE, 500);
     CommandArgument_Bool_OrDefault(AVOID_PROMISCUOUS_KMERS, False);
     CommandArgument_Bool_OrDefault(TRIM_ALIGNMENT_ENDS, False);
     CommandArgument_UnsignedInt_OrDefault(TRIM_PAR1, 3);
     CommandArgument_UnsignedInt_OrDefault(TRIM_PAR2, 5);
     CommandArgument_UnsignedInt_OrDefault(TRIM_PAR3, 2);
     CommandArgument_UnsignedInt_OrDefault(MIN_PERFECT_MATCH, 0);
     CommandArgument_Double_OrDefault(MAX_QUAL_SCORE, 0.0);
     CommandArgument_Double_OrDefault(MAX_QUAL_SCORE_FRACT, 0.0);
     CommandArgument_Bool_OrDefault(PERFECT_ALIGNMENTS, False);
     CommandArgument_Bool_OrDefault(SMITH_WATERMAN, False);
     CommandArgument_Bool_OrDefault(FIRST_COMPLETE, False);
     CommandArgument_Bool_OrDefault(SECOND_COMPLETE, False);
     CommandArgument_Bool_OrDefault(SILENT, False);
     CommandArgument_Bool_OrDefault(ALT_METHOD, False);
     CommandArgument_Bool_OrDefault(SW_GAP_METHOD, False);
     CommandArgument_UnsignedInt_OrDefault(ALT_METHOD_BANDWIDTH, 0);
     CommandArgument_UnsignedInt_OrDefault(MAKEALIGNS_MAXERRS, 200);
     CommandArgument_UnsignedInt_OrDefault(MAKEALIGNS_MAXERRS_LOCAL, 50);
     CommandArgument_UnsignedInt_OrDefault(MAKEALIGNS_MAXERRS_LOCAL_DONE, 0);
     CommandArgument_UnsignedInt_OrDefault(MAKEALIGNS_NSTRETCH, 1);
     CommandArgument_UnsignedInt_OrDefault(MAKEALIGNS_STRETCH, 3);
     CommandArgument_UnsignedInt_OrDefault(MAKEALIGNS_END_STRETCH, 4);
     CommandArgument_UnsignedInt_OrDefault(SW_GAP_MIN_MAX_MUTMER, 200);
     CommandArgument_UnsignedInt_OrDefault(SW_GAP_MAX_OFFSET_DIFF, 1000);
     CommandArgument_UnsignedInt_OrDefault(SW_GAP_MAX_GAP, 1000);
     CommandArgument_UnsignedInt_OrDefault(SW_GAP_MIN_PROG_LENGTH, 0);
     CommandArgument_Double_OrDefault(SW_GAP_MIN_PROG_RATIO, 0.8);
     CommandArgument_Double_OrDefault(SW_GAP_MIN_OVERLAP_FRAC, 0.2);
     CommandArgument_UnsignedInt_OrDefault(SW_GAP_IGNORE_OVERLAP_FRAC_LENGTH, INT_MAX);
     CommandArgument_Bool_OrDefault(SW_GAP_AFFINE_PENALTY, False);
     CommandArgument_Bool_OrDefault(SW_GAP_VERBOSE, False);
     CommandArgument_Bool_OrDefault(POLY_SCORE, False);
     CommandArgument_String_OrDefault(TEMP_ALIGNS_FILE, "");
     CommandArgument_String_OrDefault(TEMP_DIR, "");
     CommandArgument_String_OrDefault(LOG_FILE, "/dev/null");
     CommandArgument_Bool_OrDefault(COUNT_NQS, False);
     CommandArgument_UnsignedInt_OrDefault(PRINT_PERFECT_INTERVALS2, undefined);
     CommandArgument_Bool_OrDefault(PRINT_MISMATCH_POSITIONS, False);
     CommandArgument_Bool_OrDefault(RMR, False);
     CommandArgument_Double_OrDefault(MAX_RMR, undefined);
     CommandArgument_UnsignedInt_OrDefault(MAX_MISMATCHES, undefined);

     CommandArgument_Bool_OrDefault(FW_ONLY, False);

     // Output options:

     CommandArgument_Bool_OrDefault(OUT_PROPER_MATCHES, True);
     CommandArgument_Bool_OrDefault(OUT_IMPROPER_MATCHES, False);
     CommandArgument_Bool_OrDefault(OUT_MISSING1, False);
     CommandArgument_Bool_OrDefault(OUT_COVERAGE2, False);
     CommandArgument_Bool_OrDefault(OUT_READS1, False);
     CommandArgument_Bool_OrDefault(CALL_IMPROPER_PROPER, False);
     CommandArgument_Bool_OrDefault(LT, False);
     CommandArgument_Bool_OrDefault(EQ, False);
     CommandArgument_Bool_OrDefault(NE, False);
     CommandArgument_Bool_OrDefault(VISUAL_COVERAGE, False);
     CommandArgument_Bool_OrDefault(VISUAL_COVERAGE_IMPROPER, False);
     CommandArgument_UnsignedInt_OrDefault(VISUAL_COVERAGE_QUAL, 0);
     CommandArgument_String_OrDefault(BINARY_ALIGNMENTS_FILE, "");
     CommandArgument_Bool_OrDefault(BINARY_ALIGNMENTS_FILE_APPEND, False);
     CommandArgument_Bool_OrDefault(ABBREVIATE_ALIGNMENTS, True);
     CommandArgument_Bool_OrDefault(ABBREVIATE_GOOD, False);
     CommandArgument_Bool_OrDefault(SUMMARY_ALIGNMENTS, False);
     CommandArgument_Bool_OrDefault(SUMMARY_ALIGNMENTS_BRIEF, False);
     CommandArgument_Bool_OrDefault(SUMMARY_ALIGNMENTS_BRIEF_PLUS, False);
     CommandArgument_String_OrDefault(SUMMARY_ALIGNMENTS_BRIEF_PLUS_FILE, "");

     EndCommandArguments;

     longlong NQS_look = 0, NQS_see = 0;

     String PREBASE = PRE + "/" + BASE + "/";

     if ( HEAD1 != "" )
     {    if ( FILE1 != "" || QUAL1 != "" )
          {    FatalErr( "HEAD1 use precludes use of FILE1 or QUAL1." );    }
          String fasta = PREBASE + "/" + HEAD1 + ".fasta";
          String fastb = PREBASE + "/" + HEAD1 + ".fastb";
          Bool fasta_exists = IsRegularFile(fasta);
          Bool fastb_exists = IsRegularFile(fastb);
          if ( fasta_exists && fastb_exists )
          {    FatalErr( "Confused by option HEAD1: both " << fasta << " and "
                    << fastb << " exist." );    }
          if ( !fasta_exists && !fastb_exists )
          {    FatalErr( "Confused by option HEAD1: neither " << fasta << " nor "
                    << fastb << " exist." );    }
          if (fasta_exists) FILE1 = HEAD1 + ".fasta";
          else FILE1 = HEAD1 + ".fastb";
          String qualb = PREBASE + "/" + HEAD1 + ".qualb";
          if ( IsRegularFile(qualb) ) QUAL1 =  HEAD1 + ".qualb";
          String names = PREBASE + "/" + HEAD1 + ".names";
          if ( IsRegularFile(names) ) NAME1 = HEAD1 + ".names";
     }
     if ( HEAD2 != "" )
     {    if ( FILE2 != "" || QUAL2 != "" )
          {    FatalErr( "HEAD2 use precludes use of FILE2 or QUAL2." );    }
          String fasta = PREBASE + "/" + HEAD2 + ".fasta";
          String fastb = PREBASE + "/" + HEAD2 + ".fastb";
          Bool fasta_exists = IsRegularFile(fasta);
          Bool fastb_exists = IsRegularFile(fastb);
          if ( fasta_exists && fastb_exists )
          {    FatalErr( "Confused by option HEAD2: both " << fasta << " and "
                    << fastb << " exist." );    }
          if ( !fasta_exists && !fastb_exists )
          {    FatalErr( "Confused by option HEAD2: neither " << fasta << " nor "
                    << fastb << " exist." );    }
          if (fasta_exists) FILE2 = HEAD2 + ".fasta";
          else FILE2 = HEAD2 + ".fastb";
          String qualb = PREBASE + "/" + HEAD2 + ".qualb";
          if ( IsRegularFile(qualb) && QUAL1 != "" ) QUAL2 =  HEAD2 + ".qualb";
          String names = PREBASE + "/" + HEAD2 + ".names";
          if ( IsRegularFile(names) ) NAME2 = HEAD2 + ".names";
     }
     if ( FILE1 == "" )
     {    FatalErr( "Either FILE1 or HEAD1 must be specified." );    }

     ostream* out_ptr;
     String output = PRE + "/" + OUTPUT;
     if ( OUTPUT != "" ) out_ptr = new ofstream( output.c_str( ) );
     else out_ptr = &cout;
     ostream& out = *out_ptr;

     ostream* bout_ptr = 0;
     if (SUMMARY_ALIGNMENTS_BRIEF_PLUS)
          bout_ptr = new ofstream(
               SUMMARY_ALIGNMENTS_BRIEF_PLUS_FILE.c_str( ), ios::app );

     Bool self = ( FILE2 == "" );
     if (self)
     {    ForceAssert( QUAL2 == "" );
          ForceAssert( START1 == 0 && START2 == 0 );
          ForceAssert( !OUT_COVERAGE2 );
          ForceAssert( !VISUAL_COVERAGE );
          FILE2 = FILE1;
          QUAL2 = QUAL1;
          LT = True;
          if ( !SILENT ) std::cout << "doing self-comparison\n\n";    }

     ForceAssert( K == 8 || K == 12 || K == 16 || K == 24 || K == 48 );
     ForceAssert( MAX_QUAL_SCORE == 0.0 || QUAL1 != "" );
     ForceAssert( MAX_QUAL_SCORE_FRACT == 0.0 || QUAL1 != "" );
     ForceAssert( QUAL2 == "" || QUAL1 != "" );
     ForceAssert( !PERFECT_ALIGNMENTS || ( QUAL1 != "" && QUAL2 != "" ) );
     ForceAssert( TEMP_ALIGNS_FILE == "" || TEMP_DIR == "" );

     if ( ACTUAL_PASSES > PASSES )
       ACTUAL_PASSES = PASSES;

     if ( BINARY_ALIGNMENTS_FILE != "" && !BINARY_ALIGNMENTS_FILE_APPEND )
     {    Remove( PRE + "/" + BINARY_ALIGNMENTS_FILE );
          Remove( PRE + "/" + BINARY_ALIGNMENTS_FILE + ".gz" );   }

     ForceAssert( ! ( ALT_METHOD && SW_GAP_METHOD ) );
     makealigns_method *method_ptr;
     if ( ALT_METHOD )
     {
       makealigns_alt_method *alt_method_ptr = new makealigns_alt_method;

       alt_method_ptr->SetMaxErrs( MAKEALIGNS_MAXERRS );
       alt_method_ptr->SetBandwidth( ALT_METHOD_BANDWIDTH );

       method_ptr = alt_method_ptr;
     }
     else if ( SW_GAP_METHOD )
     {
       makealigns_sw_gap_method *sw_gap_method_ptr = new makealigns_sw_gap_method;

       sw_gap_method_ptr->SetMaxErrs( MAKEALIGNS_MAXERRS );
       sw_gap_method_ptr->SetEndStretch( MAKEALIGNS_END_STRETCH );
       sw_gap_method_ptr->SetMinMaxMutmerLength( SW_GAP_MIN_MAX_MUTMER );
       sw_gap_method_ptr->SetMaxMutmerOffsetDiff( SW_GAP_MAX_OFFSET_DIFF );
       sw_gap_method_ptr->SetMaxGap( SW_GAP_MAX_GAP );
       sw_gap_method_ptr->SetMinProgressionLength( SW_GAP_MIN_PROG_LENGTH );
       sw_gap_method_ptr->SetMinProgressionRatio( SW_GAP_MIN_PROG_RATIO );
       sw_gap_method_ptr->SetMinOverlapFraction( SW_GAP_MIN_OVERLAP_FRAC );
       sw_gap_method_ptr->SetIgnoreOverlapFractionLength( SW_GAP_IGNORE_OVERLAP_FRAC_LENGTH );
       sw_gap_method_ptr->SetAffinePenalties( SW_GAP_AFFINE_PENALTY );

       sw_gap_method_ptr->SetVerbose( SW_GAP_VERBOSE );

       method_ptr = sw_gap_method_ptr;
     }
     else
     {
       makealigns_orig_method *orig_method_ptr = new makealigns_orig_method;

       orig_method_ptr->SetMaxBadness( MAKEALIGNS_MAXERRS );
       orig_method_ptr->SetMaxErrs( MAKEALIGNS_MAXERRS );
       orig_method_ptr->SetLocalMaxErrs( MAKEALIGNS_MAXERRS_LOCAL );
       orig_method_ptr->SetLocalMaxErrsDone( MAKEALIGNS_MAXERRS_LOCAL_DONE );
       orig_method_ptr->SetStretch( MAKEALIGNS_STRETCH );
       orig_method_ptr->SetEndStretch( MAKEALIGNS_END_STRETCH );
       orig_method_ptr->SetNStretch( MAKEALIGNS_NSTRETCH );
       orig_method_ptr->SetCl( 20 );
       orig_method_ptr->SetMaxAlignCtorCalls( 1000000 );

       method_ptr = orig_method_ptr;
     }

     String contigs1_file = PREBASE + FILE1, contigs2_file = PREBASE + FILE2;

     // Read in the contigs, including their id's.

     vecbasevector EE, EErc;
     int N0, N1;
     {    if ( !FILE1_FASTB && !contigs1_file.Contains( ".fastb", -1 ) )
	  {    vecbasevector EE1;
               FetchReads( EE1, 0, contigs1_file, 0 );
	       if ( START1 == 0 )
	       {    EE = EE1;
	            if ( COUNT1 != 0 ) EE.resize( COUNT1 );    }
	       else
	       {    int EE1_rawsize = 0;
	            for ( int i = 0; i < (int) COUNT1; ++i )
		         EE1_rawsize += EE1[ START1 + i ].size();
	            EE.Reserve( EE1_rawsize, COUNT1 );
	            for ( int i = 0; i < (int) COUNT1; ++i )
		         EE.push_back( EE1[ START1 + i ] );    }    }
	  else
	  {    int num_contigs1 = MastervecFileObjectCount(contigs1_file);
	       if ( COUNT1 == 0 ) COUNT1 = num_contigs1 - START1;
	       if ( (int) START1 + (int) COUNT1 > num_contigs1 )
	            FatalErr( "START1=" << START1 << " and COUNT1=" << COUNT1
			 << " are inconsistent with the number of objects in "
			 << contigs1_file << ": " << num_contigs1 << "." );
	       EE.ReadRange( contigs1_file, START1, START1 + COUNT1, 0 );    }
          N0 = EE.size( );
          if ( !self && !FILE2_FASTB && !contigs2_file.Contains( ".fastb", -1 ) )
	  {    vecbasevector cons2;
	       FetchReads( cons2, 0, contigs2_file, 0 );
	       if ( COUNT2 == 0 ) COUNT2 = cons2.size();
	       if ( cons2.size( ) < START2 )
	            FatalErr( "START2 value is too large" );
	       if ( cons2.size( ) < START2 + COUNT2 )
	            FatalErr( "COUNT2 value is too large" );
	       int EE2_rawsize = 0;
	       for ( int i = 0; i < (int)COUNT2; i++ )
	            EE2_rawsize += cons2[ START2 + i ].size();
	       EE.reserve( EE.size() + COUNT2 );
	       for ( int i = 0; i < (int)COUNT2; i++ )
	            EE.push_back( cons2[ START2 + i ] );    }
	  else if ( !self )
	  {    int num_contigs2 = MastervecFileObjectCount(contigs2_file);
	       if ( COUNT2 == 0 ) COUNT2 = num_contigs2 - START2;
	       if ( (int) START2 + (int) COUNT2 > num_contigs2 )
	            FatalErr( "START2=" << START2 << " and COUNT2=" << COUNT2
		 	 << " are inconsistent with the number of objects in "
			 << contigs2_file << ": " << num_contigs2 << "." );
	       EE.ReadRange( contigs2_file, START2, START2 + COUNT2, 0 );    }
	  N1 = EE.size() - N0;    }

     EErc = EE;
     for ( size_t i = 0; i < EE.size( ) ; i++ )
          EErc[i].ReverseComplement( );

     longlong EE1_total = 0, EE2_total = 0;
     vec<int> EE1_length(N0), EE2_length(N1);
     for ( int i = 0; i < N0; i++ )
     {    EE1_total += EE[i].size( );
          EE1_length[i] = EE[i].size( );    }
     for ( int i = 0; i < N1; i++ )
     {    EE2_total += EE[ N0 + i ].size( );
          EE2_length[i] = EE[ N0 + i ].size( );    }

     vecqualvector Q;

     longlong Q_total
          = ( QUAL1.empty( ) ? 0 : EE1_total )
               + ( ( self || QUAL2.empty( ) ) ? 0 : EE2_total );
     Q.Reserve( Q_total, EE.size() );

     if ( QUAL1 != "" ) Q.ReadRange( PREBASE + QUAL1, START1, START1 + N0, 0 );
     if ( !self && QUAL2 != "" )
          Q.ReadRange( PREBASE + QUAL2, START2, START2 + N1, 0 );

     for ( vecqvec::size_type rd_idx = 0; rd_idx < Q.size( ); ++rd_idx )
          if ( Q[rd_idx].size() != EE[rd_idx].size( ) )
          {    FatalErr( "The length of sequence " << rd_idx
                    << " is " << EE[rd_idx].size( ) << ", but the length "
                    << "of the corresponding quality vector is "
                    << Q[rd_idx].size( ) << "." );    }

     vecString ids1, ids2;
     String line;

     if ( !FILE1_FASTB && !contigs1_file.Contains( ".fastb", -1 ) )
     {    Ifstream( c1, contigs1_file );
          while(1)
          {    getline( c1, line );
               if ( !c1 ) break;
	       if ( ids1.size() == START1 + N0 ) break;
               if ( !line.Contains( ">", 0 ) ) continue;
               while(1)
               {    if ( line.Contains( " ", -1 ) ) line.erase(line.size( ) - 1, 1);
                    else break;    }
               ids1.push_back( line.After( ">" ) );    }    }
     else if ( NAME1 != "" )
     {    if ( IsAsciiVec( PREBASE + NAME1 ) )
          {    READ( PREBASE + NAME1, vec<String>, slowids1 );
               ids1.assign( slowids1.begin(), slowids1.end() );    }
          else ids1.ReadAll( PREBASE + NAME1 );
          ids1.resize( START1 + N0);    }
     else
     {    int max_id1 = START1 + N0;
          ids1.resize( max_id1 );
          for ( int i = START1; i < max_id1; i++ )
               ids1[i] = ToString(i);    }

     ForceAssertEq( (int) ids1.size( ), (int) START1 + N0 );

     if ( !self && !FILE2_FASTB && !contigs2_file.Contains( ".fastb", -1 ) )
     {    Ifstream( c2, contigs2_file );
          while(1)
          {    getline( c2, line );
               if ( !c2 ) break;
	       if ( ids2.size() == START2 + N1 ) break;
               if ( !line.Contains( ">", 0 ) ) continue;
               while(1)
               {    if ( line.Contains( " ", -1 ) ) line.erase(line.size( ) - 1, 1);
                    else break;    }
               ids2.push_back( line.After( ">" ) );    }    }
     else if ( !self && NAME2 != "" )
     {    if ( IsAsciiVec( PREBASE + NAME2 ) )
          {    READ( PREBASE + NAME2, vec<String>, slowids2 );
               ids2.assign( slowids2.begin(), slowids2.end() );    }
          else ids2.ReadAll( PREBASE + NAME2 );
          ids2.resize( START2 + N1 );    }
     else if ( !self )
     {    int max_id2 = START2 + N1;
          ids2.resize( max_id2 );
          for ( int i = START2; i < max_id2; i++ )
               ids2[i] = ToString( i );    }
     else ids2 = ids1;

     if ( !self ) ForceAssertEq( (int) ids2.size( ), (int) START2 + N1 );

     int block_size = BLOCK_SIZE;

     String aligns_file;
     if ( TEMP_ALIGNS_FILE.empty() )
     {
       if ( TEMP_DIR.empty() )
         aligns_file = "/tmp";
       else
         aligns_file = TEMP_DIR;
       aligns_file += "/CmpSeq_tmp.XXXXXX";
       aligns_file = temp_file::generateName(aligns_file.c_str());
     }
     else
       aligns_file = TEMP_ALIGNS_FILE;


     Ofstream( log, LOG_FILE );

     vector<augmented_alignment> nobbits;
     int passes = (N0+block_size-1)/block_size;
     if (self) passes = 1;
     vec<Bool> good_found(N0, False), something_found(N0, false);

     vec< vec<int> > proper_second_genome_hits(N1);
     vec< StdSet<int> > proper_hits(N1);
     if ( OUT_COVERAGE2 )
       for ( unsigned int i = 0; i < proper_second_genome_hits.size( ); i++ )
       {    proper_second_genome_hits[i].resize( EE[ i + N0 ].size( ) );
            for ( unsigned int j = 0; j < EE[ i + N0 ].size( ); j++ )
	         proper_second_genome_hits[i][j] = 0;    }

     vec< vec<int> > vis;
     if (VISUAL_COVERAGE||VISUAL_COVERAGE_IMPROPER)
     {    vis.resize(N1);
          for ( int i = 0; i < N1; i++ )
               vis[i].resize( EE[ N0 + i ].size( ), 0 );    }

     const int max_bin_aligns_in_memory = 10000;
     vec<alignment_plus> bin_aligns;
     if ( BINARY_ALIGNMENTS_FILE != "" )
          bin_aligns.reserve(max_bin_aligns_in_memory);

     if ( !SILENT ) std::cout << "running " << passes << " passes" << std::endl;
     for ( int pass = 0; pass < passes; pass++ )
     {
          nobbits.clear();

          int start, EE1part;
          if ( !self )
          {    start = block_size * pass;
               EE1part = Min( block_size, N0 - start );    }
          else
          {    start = 0;
               EE1part = N0;    }
          vecbasevector EEX( EE1part + N1 );
          for ( int i = 0; i < EE1part; i++ )
               EEX[i] = EE[ start + i ];
          for ( int i = 0; i < N1; i++ )
               EEX[ i + EE1part ] = EE[ i + N0 ];

          to_compare what_to_compare( FIRST_VS_SECOND, EE1part );
          if (self) what_to_compare = to_compare( ALL_VS_ALL );

	  #define CALL_MAKE_ALIGNS(_K)                                               \
               MakeAligns<2, _K, 50>( PASSES, ACTUAL_PASSES, EEX, what_to_compare,   \
				      method_ptr, MAXCLIQ,                           \
				      aligns_file, log, MAX_ALIGNS, MIN_MUTMER,      \
				      False, AVOID_PROMISCUOUS_KMERS )
	  DISPATCH_ON_K( K, CALL_MAKE_ALIGNS );

          int aligns_length;
          ifstream aligns_in( aligns_file.c_str( ) );
          int id1, id2;
          Bool rc = False;

	  while (1)
          {    BinRead( aligns_in, aligns_length );
	       if ( !aligns_in ) break;
               for ( int ll = 0; ll < aligns_length; ll++ )
               {    alignment this_align;
	            this_align.Read( aligns_in, id1, id2, rc );
                    if ( rc && FW_ONLY ) continue;

                    int length1 = EEX[id1].size( ), length2 = EEX[id2].size( );
                    int pos1 = this_align.pos1( ), pos2 = this_align.pos2( );
                    int Pos1 = this_align.Pos1( ), Pos2 = this_align.Pos2( );

                    if ( id1 < id2 )
                    {    if ( !rc )
                         {    nobbits.push_back( augmented_alignment( this_align, rc,
                                   length1, length2, id1 + start,
                                   id2 - EE1part + N0, pos1, pos2,
                                   this_align.Errors( ), Pos1, Pos2 ) );    }
                         else
                         {    alignment r = this_align;
                              r.ReverseThis( length1, length2 );
                              nobbits.push_back( augmented_alignment( r, rc, length1,
                                   length2, id1 + start, id2 - EE1part + N0,
                                   r.pos1( ), r.pos2( ),
                                   r.Errors( ), r.Pos1( ), r.Pos2( ) ) );    }    }
                    else
                    {    alignment a = this_align;
                         a.Flip( );
                         nobbits.push_back( augmented_alignment( a, rc, length2,
                              length1, id2 + start, id1 - EE1part + N0, pos2, pos1,
                              a.Errors( ), Pos2, Pos1 ) );    }    }    }

          for ( unsigned int i = 0; i < nobbits.size( ); i++ )
          {    if ( EstimatedOverlap( nobbits[i].a,
 		         EE[nobbits[i].id1], EE[nobbits[i].id2] ) > 0 )
                    good_found[nobbits[i].id1] = True;    }

          set<int> ignored_reads;
          for ( unsigned int i = 0; i < nobbits.size( ); i++ )
          {    int id1 = nobbits[i].id1, id2 = nobbits[i].id2;
               alignment a = nobbits[i].a;
               if ( OMIT_PROPER_ALIGNERS
                    && EstimatedOverlap( a, EE[id1], EE[id2] ) > 0 )
                    ignored_reads.insert(id1-start);
               if ( EE[id1].size( ) < MIN_SIZE1 )
                    ignored_reads.insert(id1-start);    }

          if ( MIN_COVERAGE1 > 0 )
          {    vec< vec<Bool> > hit(EE1part);
               for ( int i = 0; i < EE1part; i++ )
               {    hit[i].resize( EEX[i].size( ) );
                    for ( unsigned int j = 0; j < hit[i].size( ); j++ )
                         hit[i][j] = False;    }
               for ( unsigned int i = 0; i < nobbits.size( ); i++ )
               {    alignment a = nobbits[i].a;
                    int id1 = nobbits[i].id1;
                    int pos1 = nobbits[i].pos1, Pos1 = nobbits[i].Pos1;

                    if (USE_ONLY_SEMIPROPER_ALIGNMENTS)
                    {    if ( pos1 > 1 && Pos1 < (int) EE[id1].size( ) - 1 )
                              continue;    }
                    if ( pos1 == Pos1 || double(nobbits[i].score)/double(Pos1-pos1)
                         > MAX_ERROR_RATE )
                         continue;

                    for ( int j = pos1; j < Pos1; j++ )
                         hit[id1-start][j] = True;    }
               for ( int i = 0; i < EE1part; i++ )
               {    int hits = 0;
                    for ( unsigned int j = 0; j < hit[i].size( ); j++ )
                         if ( hit[i][j] ) ++hits;
                    if ( hit[i].size( ) == 0
                         || double(hits) / double( hit[i].size( ) ) < MIN_COVERAGE1 )
                         ignored_reads.insert(i);    }    }

          for ( unsigned int i = 0; i < nobbits.size( ); i++ )
          {    alignment a = nobbits[i].a;
               static align aa, aar;
               aa.UnpackFrom(a);
               aar = aa;
               aar.ReverseThis( nobbits[i].length1, nobbits[i].length2 );
               // TODO: potentially dangerous truncation of index values by all these int vars
               int RC = nobbits[i].RC, id1 = nobbits[i].id1, id2 = nobbits[i].id2;
               if ( Member( ignored_reads, id1-start ) ) continue;
               int pos1 = nobbits[i].pos1, Pos1 = nobbits[i].Pos1,
                    pos2 = nobbits[i].pos2, Pos2 = nobbits[i].Pos2;
               float score = nobbits[i].score;
               something_found[id1] = True;
               if ( OMIT_IF_PROPER_OVERLAP_EXISTS == "True" && good_found[id1] )
                    continue;

               if ( Pos2 - pos2 < (int) MIN_OVERLAP ) continue;

               if (FIRST_COMPLETE)
               {    if ( pos1 > 0 || Pos1 < (int) EE[id1].size( ) ) continue;    }
               if (SECOND_COMPLETE)
               {    if ( pos2 > 0 || Pos2 < (int) EE[id2].size( ) ) continue;    }

               if (SMITH_WATERMAN)
               {    if ( !RC )
                    {    SmithWatBandedA( EE[id1], EE[id2], a.Offset( ),
                              Bandwidth(a), a );
                         score = ActualErrors( EE[id1], EE[id2], a );    }
                    else
                    {    SmithWatBandedA( EErc[id1], EE[id2], a.Offset( ),
                              Bandwidth(a), a );
                         score = ActualErrors( EErc[id1], EE[id2], a );    }
                    aa.UnpackFrom(a);
                    aar = aa;
                    aar.ReverseThis( nobbits[i].length1, nobbits[i].length2 );
                    pos1 = a.pos1( );
                    pos2 = a.pos2( );
                    Pos1 = a.Pos1( );
                    Pos2 = a.Pos2( );    }

               if (PERFECT_ALIGNMENTS)
               {    if ( !RC )
                         PerfectAlignment( EE[id1], Q[id1], EE[id2], Q[id2], aa,
                              score );
                    else
                    {    static qualvector Qrcid1;
                         Qrcid1.SetToReverseOf( Q[id1] );
                         PerfectAlignment( EErc[id1], Qrcid1, EE[id2], Q[id2], aa,
                              score );    }
                    aa.Compactify( EE[id1].size( ), EE[id2].size( ) );
                    a.Set(aa);
                    pos1 = aa.pos1( );
                    pos2 = aa.pos2( );
                    Pos1 = aa.Pos1( );
                    Pos2 = aa.Pos2( );    }

               if (TRIM_ALIGNMENT_ENDS)
               {    TrimAlignmentEnds( aa, RC, EE[id1], EE[id2], True, True,
                         TRIM_PAR1, TRIM_PAR2, TRIM_PAR3 );
                    a.Set(aa);
                    nobbits[i].pos1 = a.pos1( );
                    nobbits[i].pos2 = a.pos2( );
                    nobbits[i].Pos1 = a.Pos1( );
                    nobbits[i].Pos2 = a.Pos2( );
                    pos1 = nobbits[i].pos1;
                    pos2 = nobbits[i].pos2;
                    Pos1 = nobbits[i].Pos1;
                    Pos2 = nobbits[i].Pos2;
                    if ( !RC ) score = ActualErrors( EE[id1], EE[id2], a );
                    else score = ActualErrors( EErc[id1], EE[id2], a );    }

               if ( Pos2 - pos2 < (int) MIN_OVERLAP ) continue;

               if ( MIN_PERFECT_MATCH > 0 )
               {    if ( !RC && MaxPerfectMatch( False, a, EE[id1], EE[id2] )
                         < (int) MIN_PERFECT_MATCH ) continue;
                    if ( RC && MaxPerfectMatch( False, a, EErc[id1], EE[id2] )
                         < (int) MIN_PERFECT_MATCH ) continue;    }

               if ( MAX_MISMATCHES != undefined )
               {    static vector<int> mgg;
                    if ( !RC ) mgg = a.MutationsGap1Gap2( EE[id1], EE[id2] );
                    else mgg = a.MutationsGap1Gap2( EErc[id1], EE[id2] );
                    if ( mgg[0] > (int) MAX_MISMATCHES ) continue;    }

               if (USE_ONLY_SEMIPROPER_ALIGNMENTS)
               {    if ( pos1 > 1 && Pos1 < (int) EE[id1].size( ) - 1 )
                         continue;    }
               if ( pos1 == Pos1 || double(nobbits[i].score)/double(Pos1-pos1)
                    > MAX_ERROR_RATE )
                    continue;

               if ( MAX_RMR != (double) undefined )
               {    const basevector &rd1 = ( !RC ? EE[id1] : EErc[id1] );
                    const basevector &rd2 = EE[id2];
                    if ( Float(100) * ReciprocalMatchRate( aa, rd1, rd2 )
                         > Float(MAX_RMR) )
                         continue;    }

	       int rel_id2 = id2 - N0;
	       int CONTIG1 = START1 + id1;
               int CONTIG2;
               if ( !self ) CONTIG2 = START2 + rel_id2;
               else CONTIG2 = id2;
               if ( LT && CONTIG1 >= CONTIG2 ) continue;
               if ( EQ && CONTIG1 != CONTIG2 ) continue;
               if ( NE && CONTIG1 == CONTIG2 ) continue;
               int est_over = EstimatedOverlap( a, EE[id1], EE[id2] );

               float qual_score = -1.0, mismatch_rate = -1.0;
               if ( QUAL1 != "" )
               {    if ( QUAL2 == "" )
                    {
                         const basevector &rd1 = ( !RC ? EE[id1] : EErc[id1] );
                         const basevector &rd2 = EE[id2];

                         int mismatches = 0, good_bases = 0;
                         int j, p1 = aa.pos1( ), p2 = aa.pos2( );
                         int nblocks = aa.Nblocks( );
                         for ( j = 0; j < nblocks; j++ )
                         {    if ( aa.Gaps(j) > 0 ) p2 += aa.Gaps(j);
                              if ( aa.Gaps(j) < 0 ) p1 -= aa.Gaps(j);
                              for ( int x = 0; x < aa.Lengths(j); x++ )
                              {    int q1p1
                                        = ( !RC ? Q[id1][p1]
                                             : Q[id1][ Q[id1].size( ) - p1 - 1 ] );
                                   if ( rd1[p1] != rd2[p2] && q1p1 >= 50 )
                                        ++mismatches;
                                   if ( q1p1 >= 50 ) ++good_bases;
                                   ++p1;
                                   ++p2;    }    }
                         if ( good_bases != 0 )
                              mismatch_rate = float(mismatches)/float(good_bases);

                         if ( !POLY_SCORE )
                         {    if ( !RC )
                                   qual_score = ScoreAlignment( aa, EE[id1], Q[id1],
                                        EE[id2] );
                              else qual_score = ScoreAlignment( True, aar, EE[id1],
                                        Q[id1], EE[id2] );
                              if (MAX_QUAL_SCORE > 0 && qual_score > MAX_QUAL_SCORE)
                                   continue;
                              if ( MAX_QUAL_SCORE_FRACT > 0
                                   && float(qual_score)/float(Pos2 - pos2)
                                        > MAX_QUAL_SCORE_FRACT )
                              {    continue;    }    }
                         else
                         {    if ( !RC )
                                   qual_score = ScoreAlignmentPoly( aa, EE[id1],
                                        Q[id1], EE[id2] );
                              else qual_score = ScoreAlignmentPoly( True, aar,
                                        EE[id1], Q[id1], EE[id2] );
                              if (MAX_QUAL_SCORE > 0 && qual_score > MAX_QUAL_SCORE)
                                   continue;
                              if ( MAX_QUAL_SCORE_FRACT > 0
                                   && float(qual_score)/float(Pos2 - pos2)
                                        > MAX_QUAL_SCORE_FRACT )
                              {    continue;    }    }    }
                    else
                    {
                         if (COUNT_NQS)
                         {    const basevector &rd1 = ( !RC ? EE[id1] : EErc[id1] );
                              const basevector &rd2 = EE[id2];
                              int look, see;
                              static qualvector Qrcid1;
                              if (RC) Qrcid1.SetToReverseOf( Q[id1] );
                              const qualvector& q1 = ( !RC ? Q[id1] : Qrcid1 );
                              const qualvector& q2 = Q[id2];
                              CountNQS( a, rd1, rd2, q1, q2, 30, 25, 1, look, see );
                              NQS_look += look;
                              NQS_see += see;    }

                         if ( !POLY_SCORE )
                         {    if ( !RC )
                              {    Regap( aa, EE[id1], Q[id1], EE[id2], Q[id2] );
                                   qual_score = ScoreAlignment( aa, EE[id1], Q[id1],
                                        EE[id2], Q[id2] );    }
                              else
                              {    Regap( True, aar, EE[id1], Q[id1], EE[id2],
                                        Q[id2] );
                                   qual_score = ScoreAlignment( True, aar, EE[id1],
                                        Q[id1], EE[id2], Q[id2] );
                                   aa = aar;
                                   aa.ReverseThis( nobbits[i].length1,
                                        nobbits[i].length2 );    }    }
                         else
                         {    if ( !RC )
                              {    Regap( aa, EE[id1], Q[id1], EE[id2], Q[id2] );
                                   qual_score = ScoreAlignmentPoly( aa, EE[id1],
                                        Q[id1], EE[id2], Q[id2] );    }
                              else
                              {    Regap( True, aar, EE[id1], Q[id1], EE[id2],
                                        Q[id2] );
                                   qual_score = ScoreAlignmentPoly( True, aar,
                                        EE[id1], Q[id1], EE[id2], Q[id2] );
                                   aa = aar;
                                   aa.ReverseThis( nobbits[i].length1,
                                        nobbits[i].length2 );    }    }
                         if ( MAX_QUAL_SCORE > 0 && qual_score > MAX_QUAL_SCORE )
                              continue;
                         if ( MAX_QUAL_SCORE_FRACT > 0
                              && float(qual_score)/float(Pos2 - pos2)
                                   > MAX_QUAL_SCORE_FRACT )
                         {    continue;    }    }    }

               if ( VISUAL_COVERAGE_IMPROPER || VISUAL_COVERAGE && est_over > 0 )
               {    basevector& rd1 = RC ? EErc[id1] : EE[id1];
                    qualvector qs1;
                    if ( Q.size() > static_cast<size_t>(id1) ) qs1 = Q[id1];
                    if ( RC ) qs1.ReverseMe();
                    basevector& rd2 = EE[id2];
                    int p1 = aa.pos1( ), p2 = aa.pos2( );
                    int nblocks = aa.Nblocks( );
                    for ( int j = 0; j < nblocks; j++ )
                    {    if ( aa.Gaps(j) > 0 )
                         {    if ( vis[rel_id2][p2] == 0 ) vis[rel_id2][p2] = 1;
                              p2 += aa.Gaps(j);    }
                         if ( aa.Gaps(j) < 0 )
                         {    if ( vis[rel_id2][p2] == 0 ) vis[rel_id2][p2] = 1;
                              p1 -= aa.Gaps(j);    }
                         for ( int x = 0; x < aa.Lengths(j); x++ )
                         {    if ( rd1[p1] == rd2[p2] &&
                                   ( qs1.size() == 0 || qs1[p1] > VISUAL_COVERAGE_QUAL ) )
                                vis[rel_id2][p2] = 2;
                              else if ( vis[rel_id2][p2] == 0 ) vis[rel_id2][p2] = 1;
                              ++p1;
                              ++p2;    }    }    }

               if ( (OUT_PROPER_MATCHES && est_over > 0)
                    || (OUT_IMPROPER_MATCHES && est_over == 0) )
               {
                    if ( BINARY_ALIGNMENTS_FILE != "" )
                    {    float s = ( QUAL1 != "" ) ? qual_score : score;
                         alignment b;
                         b.Set(aa);
                         if (RC) b.ReverseThis(EE[id1].size( ), EE[id2].size( ));
                         alignment_plus ap( CONTIG1, CONTIG2, EE[id1].size( ),
                              EE[id2].size( ), RC, b, s );
                         bin_aligns.push_back(ap);
                         if ( (int) bin_aligns.size( ) >= max_bin_aligns_in_memory )
                         {    WriteAppend( PRE + "/" + BINARY_ALIGNMENTS_FILE,
                                   bin_aligns );
                              bin_aligns.clear( );   }    }

                    else if (SUMMARY_ALIGNMENTS)
                    {    static vector<int> mgg;
                         if ( !RC ) mgg = a.MutationsGap1Gap2( EE[id1], EE[id2] );
                         else mgg = a.MutationsGap1Gap2( EErc[id1], EE[id2] );
                         out << (RC ? "rc to " : "")
                              << ids1[CONTIG1] << "[" << pos1 << "-" << Pos1 - 1
                              << "] vs " << ids2[CONTIG2] << "["
                              << pos2 << "-" << Pos2 - 1 << "]; "
                              << Pos1 - pos1 << " bases, "
                              << Sum(mgg) << " errors: "
                              << mgg[0] << " mismatches, "
                              << mgg[1] << " deletions, "
                              << mgg[2] << " insertions\n";     }

                    else if (SUMMARY_ALIGNMENTS_BRIEF)
                    {    out << (RC ? "1 " : "0 ")
                              << ids1[CONTIG1] << " " << pos1 << " " << Pos1 - 1
                              << " " << ids2[CONTIG2] << " "
                              << pos2 << " " << Pos2 - 1 << " "
                              << Pos1 - pos1 << " " << score
                              << " " << qual_score << "\n";    }

                    else
                    {
                    out << "\n***** " << (RC ? "rc to " : "")
                         << "first sequence " << CONTIG1
                         << " vs second sequence " << CONTIG2
                         << " *****\n\n";
                    out << "[" << ids1[CONTIG1] << " vs " << ids2[CONTIG2] << "]\n";
                    out << "overlap goes from " << pos1 << " to " << Pos1 - 1
                         << " on the " << (RC ? "rc to the " : "" )
                         << "first sequence\n";
                    out << "overlap goes from " << pos2 << " to " << Pos2 - 1
                         << " on the second sequence\n";
                    out << "first sequence has length " << EE[id1].size( ) << "\n";
                    out << "second sequence has length " << EE[id2].size( )
                         << "\n";
                    out << "overlap_length = " << Pos2 - pos2 << "\n";

                    if (RMR)
                    {    out << "rmr:\n";
                         const basevector &rd1 = ( !RC ? EE[id1] : EErc[id1] );
                         const basevector &rd2 = EE[id2];
                         out << Float(100) * ReciprocalMatchRate( aa, rd1, rd2 )
                              << "%\n";    }

                    if (PRINT_MISMATCH_POSITIONS)
                    {    out << "mismatches:\n";
                         const basevector &rd1 = ( !RC ? EE[id1] : EErc[id1] );
                         const basevector &rd2 = EE[id2];
                         int p1 = aa.pos1( ), p2 = aa.pos2( );
                         for ( int j = 0; j < aa.Nblocks( ); j++ )
                         {    if ( aa.Gaps(j) > 0 ) p2 += aa.Gaps(j);
                              if ( aa.Gaps(j) < 0 ) p1 -= aa.Gaps(j);
                              for ( int x = 0; x < aa.Lengths(j); x++ )
                              {    if ( rd1[p1] != rd2[p2] )
                                   {    out << "(" << ( !RC ? p1 : rd1.size( )
                                             - p1 - 1 ) << "," << p2 << ")\n";    }
                                   ++p1;
                                   ++p2;    }    }    }

                    if ( PRINT_PERFECT_INTERVALS2 != undefined )
                    {    const basevector &rd1 = ( !RC ? EE[id1] : EErc[id1] );
                         static vec<ho_interval> perfs;
                         aa.PerfectIntervals2( rd1, EE[id2], perfs );
                         out << "intervals of perfect match on second sequence"
                              << " of size >= " << PRINT_PERFECT_INTERVALS2 << ":";
                         for ( int u = 0; u < perfs.isize( ); u++ )
                         {    if ( perfs[u].Length( )
                                   >= (int) PRINT_PERFECT_INTERVALS2 )
                              {    out << " [" << perfs[u].Start( )
                                        << "," << perfs[u].Stop( ) << ")";    }    }
                         out << "\n";    }

                    out << "score = " << score << "\n";
                    if ( QUAL1 != "" ) out << "qual_score = " << qual_score << "\n";
                    if ( mismatch_rate >= 0 )
                         out << setprecision(2) << 100.0 * mismatch_rate
                              << "% mismatch rate at bases of quality " << ">= 50\n";
                    if ( est_over == 0 )
                    {    out << "overlap goes from position " << pos1 << " to "
                              << Pos1-1 << " on the first sequence\n";
                         out << "THIS IS A PECULIAR OVERLAP!\n\n";    }
                    if ( Q.empty() )
                    {    if ( !RC )
                              PrintVisualAlignment( ABBREVIATE_ALIGNMENTS, out,
                                   EE[id1], EE[id2], aa );
                         else PrintVisualAlignment( ABBREVIATE_ALIGNMENTS, out,
                                   EErc[id1], EE[id2], aa,
                                   qualvector(0), qualvector(0),
                                   0, False, 0, False, 2.0, ABBREVIATE_GOOD );    }
                    else if ( QUAL2 == "" )
                    {    if ( !RC )
                              PrintVisualAlignment( ABBREVIATE_ALIGNMENTS, out,
                                   EE[id1], EE[id2], aa, Q[id1],
                                   qualvector(0),
                                   0, False, 0, False, 2.0, ABBREVIATE_GOOD );
                         else
                         {    static qualvector Qrcid1;
                              Qrcid1.SetToReverseOf( Q[id1] );
                              PrintVisualAlignment( ABBREVIATE_ALIGNMENTS, out,
                                   EErc[id1], EE[id2], aa, Qrcid1,
                                   qualvector(0),
                                   0, False, 0, False, 2.0, ABBREVIATE_GOOD );  }  }
                    else
                    {    if ( !RC )
                              PrintVisualAlignment( ABBREVIATE_ALIGNMENTS, out,
                                   EE[id1], EE[id2], aa, Q[id1], Q[id2],
                                   0, False, 0, False, 2.0, ABBREVIATE_GOOD );
                         else
                         {    static qualvector Qrcid1;
                              Qrcid1.SetToReverseOf( Q[id1] );
                              PrintVisualAlignment( ABBREVIATE_ALIGNMENTS, out,
                                   EErc[id1], EE[id2], aa, Qrcid1, Q[id2],
                                   0, False, 0, False, 2.0, ABBREVIATE_GOOD );  }  }

                    if (SUMMARY_ALIGNMENTS_BRIEF_PLUS)
                    {    (*bout_ptr) << (RC ? "1 " : "0 ")
                              << ids1[CONTIG1] << " " << pos1 << " " << Pos1 - 1
                              << " " << ids2[CONTIG2] << " "
                              << pos2 << " " << Pos2 - 1 << " "
                              << Pos1 - pos1 << " " << score
                              << " " << qual_score << "\n";    }    }    }

               if ( !self && ( est_over > 0 || CALL_IMPROPER_PROPER ) )
               {    proper_hits[rel_id2].insert(id1);
	            if ( OUT_COVERAGE2 )
		         for ( int j = pos2; j < Pos2; j++ )
                              ++proper_second_genome_hits[rel_id2][j];    }    }

          if ( !SILENT ) Dot( cout, pass );    }

     Remove(aligns_file);

     if (COUNT_NQS)
     {    PRINT2( NQS_see, NQS_look );
          cout << "NQS(30, 25) SNP rate = " << PERCENT_RATIO( 5, NQS_see, NQS_look )
               << "\n";    }

     if ( BINARY_ALIGNMENTS_FILE != "" )
          WriteAppend( PRE + "/" + BINARY_ALIGNMENTS_FILE, bin_aligns );

     if (VISUAL_COVERAGE||VISUAL_COVERAGE_IMPROPER)
     {    for ( int i = 0; i < N1; i++ )
          {    out << "Visual coverage: - denotes a correctly hit base,\n"
                    << "                 * denotes an incorrectly hit base\n\n";
               out << "Contig " << i + START2 << " has " << vis[i].size( )
                    << " bases, covered as:\n\n";
               int uncovered = 0, covered_correctly = 0, covered_incorrectly = 0;
               for ( unsigned int j = 0; j < vis[i].size( ); j++ )
               {
                    if ( j % 80 == 0 )
                    {    unsigned int k;
                         for ( k = j; k < vis[i].size( ); k++ )
                              if ( vis[i][k] != 2 ) break;
                         if ( k - j >= 80 )
                         {    int n = (k - j) - ( (k - j) % 80 );
                              out << "\n(" << n << " matching bases)";
                              covered_correctly += n;
                              j += n - 1;
                              continue;    }    }

                    if ( j % 80 == 0 )
                    {    unsigned int k;
                         for ( k = j; k < vis[i].size( ); k++ )
                              if ( vis[i][k] != 0 ) break;
                         if ( k - j >= 80 )
                         {    int n = (k - j) - ( (k - j) % 80 );
                              out << "\n(" << n << " MISSING bases)";
                              uncovered += n;
                              j += n - 1;
                              continue;    }    }

                    if ( j > 0 && j % 80 == 0 ) out << "\n";

                    if ( vis[i][j] == 0 )
                    {    out << " ";
                         ++uncovered;    }
                    else if ( vis[i][j] == 1 )
                    {    out << "*";
                         ++covered_incorrectly;    }
                    else
                    {    out << "-";
                         ++covered_correctly;    }    }
               out << "\n\n";
               out << PERCENT_RATIO( 3, uncovered, (int) vis[i].size( ) )
                    << " of bases are uncovered\n";
               out << PERCENT_RATIO( 3, covered_correctly, (int) vis[i].size( ) )
                    << " of bases are covered correctly\n";
               out << PERCENT_RATIO( 3, covered_incorrectly, (int) vis[i].size( ) )
                    << " of bases are covered incorrectly\n\n";    }    }

     if (OUT_MISSING1)
     {    out << "First genome contigs for which no matches at all were found:\n";
          vec<int> not_found;
          for ( int i = 0; i < N0; i++ )
               if ( !something_found[i] ) not_found.push_back(i);
          PrettyPrint( out, not_found );
          out << "\n\n";    }

     if (OUT_COVERAGE2)
     {    out << "Coverage of second genome contigs by reads which properly "
                << "overlap it:\n\n";
          for ( int i = 0; i < N1; i++ )
          {    int zeroes = 0;
               for ( unsigned int j = 0; j < proper_second_genome_hits[i].size( );
                    j++ )
                    if ( proper_second_genome_hits[i][j] == 0 ) ++zeroes;
               out << "contig " << i + START2<< " has " << EE[N0 + i].size( )
                    << " bases, of which " << zeroes
                    << " are not covered.\nThe number of reads which overlap it is "
                    << proper_hits[i].size( ) << ".\n\n";
               out << "missing bases:\n\n";
               for ( unsigned int j = 0; j < proper_second_genome_hits[i].size( );
                    j++ )
                    if ( proper_second_genome_hits[i][j] == 0 )
                    {    out << "> ";
                         int count = 4;
                         for ( ; j < proper_second_genome_hits[i].size( )
                              && proper_second_genome_hits[i][j] == 0; j++ )
                         {    out << as_base( EE[N0 + i][j] );
                              ++count;
                              if ( count == 80 )
                              {    out << "\n";
                                   count = 0;    }    }
                         if ( count != 0 ) out << "\n";
                         out << "\n";
                         j--;    }
               out << "\n";
               if (OUT_READS1)
               {    out << "overlapping reads:\n";
                    vec<int> v;
                    for ( set<int>::iterator j = proper_hits[i].begin( );
                         j != proper_hits[i].end( ); j++ )
                         v.push_back( *j );
                    PrettyPrint( out, v );
                    out << "\n";    }    }    }    }
