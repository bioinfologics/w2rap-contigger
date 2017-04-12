///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Discovar messages are tagged "DISCOVAR MESSAGE".

#include <omp.h>

#include "CoreTools.h"
#include "FastIfstream.h"
#include "ParseSet.h"
#include "VecUtilities.h"
#include "lookup/SAM.h"
#include "paths/long/DiscovarTools.h"

namespace DiscovarTools{

void DiscovarUnhappy( )
{    // DISCOVAR MESSAGE
     std::cout << "Sorry, Discovar cannot proceed.\n" << std::endl;
     Scram(1);    }

void ExitAssemblyEmpty( )
{    // DISCOVAR MESSAGE
     std::cout << "\nDear user, we are sad to report that your assembly is now "
          << "empty.\nThis could be because your coverage is too low, or "
          << "that something\nelse is wrong with your data.  Assembly of very "
          << "small regions\n(e.g. ~1 kb) can also cause empty assemblies.\n"
          << "Discovar complete.  Have a nice day.\n\n";
     Scram( );    }

void ExitPathsEmpty( )
{    // DISCOVAR MESSAGE
     std::cout << "\nDear user, we are sad to report the current assembly is not supported by reads.\n"
          << "This could happen if the reads do not align well with the reference (if supplied), or "
          << "that something\nelse is wrong with your data.\n"
          << "Discovar complete.  Have a nice day.\n\n";
     Scram( );    }


void ExitNoCorrectedReads( )
{    // DISCOVAR MESSAGE
     std::cout << "\nDear user, we are sad to report that after error correction, "
          << "no reads remain.\n"
          << "Further assembly is thus impossible.\n"
          << "Discovar complete.  Have a nice day.\n\n";
     Scram( );    }



void ExitShortReads(const String& additional_info )
{
    // DISCOVAR MESSAGE
    std::cout << "\nDiscovar has found that all reads in your input are too short";
    if ( additional_info != "" )
        std::cout << ":\n" << additional_info << "\n";
    std::cout << "\nDiscovar complete.  Have a nice day.\n\n";
    Scram( );    }


void CheckDiscovarRegions( const String& REGIONS )
{    vec<String> regions;
     size_t region_size=0;
     ParseStringSet( "{"+REGIONS +"}", regions );
     for ( auto r : regions )
     {    if ( r.Contains( ":" ) )
          {    String range = r.After( ":" );
               Bool bad = False;
               if ( !range.Contains( "-" ) ) bad = True;
               if ( !bad )
               {    if ( !range.Before( "-" ).IsInt( ) 
                    || !range.After( "-" ).IsInt( ) )
                    bad = True;    }
               if ( !bad )
               {    int start = range.Before( "-" ).Int( );
                    int stop = range.After( "-" ).Int( );
                    if (  start > stop  || stop <0 || start<0){
                        bad = True;
                    }
                    else{
                        region_size += stop - start;
                    }
               }
               if ( !bad )
                    for ( auto c : range ) if ( isalpha(c) ) bad = True;
               if (bad)
               {    // DISCOVAR MESSAGE
                    std::cout << "REGIONS argument " << r << " doesn't make sense: "
                         << "there needs to be a\nstart-stop integer range after "
                         << "the colon, with start <= stop." << std::endl;
                    DiscovarUnhappy( );    }    }
     }
     CheckRegionSize(region_size);
}
               

void DiscovarRefTraceControl::CheckConsistency() const{
    bool bConsistent =   getRefSeqs().size()==getRefTags().size()
                      && ref_seqs_ext.size() == getRefSeqs().size()
                      && ref_seqs_ext_ext.size() == ref_seqs_ext.size()
                      && getRefStart().size() == ref_seqs_ext_ext.size()
                      && ref_length.size() == getRefStart().size()
                      && getRefIndex().size() == ref_length.size();
    ForceAssert(bConsistent);
    for(size_t ii=0; ii<getRefSeqs().size(); ++ii){
        ForceAssert( ref_seqs_ext[ii].size() >= getRefSeqs()[ii].size() );
        if( ref_seqs_ext[ii].size() > getRefSeqs()[ii].size() + 2*EXT_LENGTH){
            std::cout << "WARNING: extended reference sequence is too long." << std::endl;
        }
    }
};

DiscovarRefTraceControl::DiscovarRefTraceControl(const String& REF_FASTA, const String& REGIONS, const long_logging& logc
                                                , const String& sVariantOutFile_
                                                ):base_t(sVariantOutFile_)
{
    bRefTraceOn = REF_FASTA != "";
    if(bRefTraceOn) {
        if( REF_FASTA == fastaindex::filename_fix(REF_FASTA)){
            std::cout << "Please create a FASTA index for the reference genome." << std::endl;
            std::cout << "    e.g.    samtools faidx " << REF_FASTA << std::endl << std::endl;
            DiscovarUnhappy( );
        }
        CheckDiscovarRegions(REGIONS);
        vec<String> regions;
        ParseStringSet( "{" + REGIONS + "}", regions );

        fastaindex fai( REF_FASTA );


        if (logc.STATUS_LOGGING)
        {    std::cout << Date() << ": RefTrace is on. Reading " << REF_FASTA 
                  << " for REGIONS=" << REGIONS << std::endl;    }

        size_t region_size=0;

        for(const auto& entry: regions){
            String chr;
	    int start=-1,end=-1;
            if( entry.Contains(":")){
        	chr = entry.Before(":");
        	//note that this is using using base-1 coordinate as base-0 coordinate.
        	//This is a decision made by the group to achieve compliance with LongProto
                start = entry.After(":").Before("-").Int();
                end = entry.After(":").After("-").Int();
            }
            else{
        	chr = entry;
                start=0;
                end = -1;
            }
            if ( !fai.has_key(chr) ) {    // DISCOVAR MESSAGE
                std::cout << "There appears to be an incompatibility between your "
                     << "REF_FASTA and REGIONS arguments.\nSpecifically, the "
                     << "region " << entry << " refers to a reference record,\n"
                     << "which is not found in the FASTA "
                     << "file\n\n" << REF_FASTA << "." << "\n" << std::endl;
                DiscovarUnhappy( );
            }
            int ilen = fai[chr].len;

            if( end!=-1){
                //note that this is a hack for using base-1 coordinate as base-0 coordinate.
                //This is a decision made by the group to achieve compliance with LongProto
                if (ilen==end+1) {--end;};
            }
            else{
        	end = ilen;
            }
            if( end > ilen ) {
                std::cout << "There appears to be an incompatibility between your "
                     << "REF_FASTA and REGIONS arguments.\nSpecifically, the "
                     << "region " << entry << " refers to a sequence,\n"
                     << "which is longer than that specified in the FASTA "
                     << "file\n\n" << REF_FASTA << "." << "\n" << std::endl;
                DiscovarUnhappy( );
            }


            if( end-start < MIN_LENGTH){
                std::cout << "Reference sequence shorter than " << MIN_LENGTH << "-base is not allowed." << std::endl;
                DiscovarUnhappy( );

            }

            getRefTags().push_back(chr);

            if (logc.STATUS_LOGGING)
            {
            std::cout << Date() << ": about to load from indexed Fasta: " << getRefTags().back() << ":" <<
        	    start << "-" << end << std::endl;
            }
            int ext_start = std::max(0,start-EXT_LENGTH);
            int ext_end = std::min(ilen,end+EXT_LENGTH);

            fastavector fv;
            LoadRegionFromIndexedFastaFile(REF_FASTA,fai, chr, ext_start, ext_end, fv);
            ref_seqs_ext.push_back(fv);
            ref_seqs_ext_ext.push_back(start-ext_start);

            getRefSeqs().push_back(fastavector());
            getRefSeqs().back().SetToSubOf(ref_seqs_ext.back(),start-ext_start,end-start);

            region_size+=end-start;
            getRefIndex().push_back(fai[chr].idx);
            getRefStart().push_back(start);
            ref_length.push_back(end-start);
        }

        // if region not specified
        if( getRefSeqs().size()==0){
	    vec<fastavector> loc_seqs;
	    vec<String> loc_tags;
	    LoadFromFastaFile(REF_FASTA,loc_seqs,loc_tags);

            getRefStart().resize(loc_tags.size(),0);
            ref_length.resize(loc_tags.size());
            getRefIndex().resize(loc_seqs.size());
            for(size_t ii=0;ii<loc_seqs.size();++ii){
                getRefIndex()[ii]=ii;
                getRefStart()[ii]=0;
                ref_length[ii]=loc_seqs[ii].size();
                region_size+=loc_seqs[ii].size();
            }
            ref_seqs_ext = loc_seqs;
            ref_seqs_ext_ext.assign( ref_seqs_ext.size(),0);
            swap(loc_seqs,getRefSeqs());
            swap(loc_tags,getRefTags());
        }
        CheckRegionSize(region_size,"When reference tracing is on, ");
    }
    CheckConsistency();
};

void CheckRegionSize(size_t size, const String& sPrefix){
    if( size > DiscovarTools::MaxRegionSize){
        std::cout << sPrefix
                  << "REGIONS totaled more than " << DiscovarTools::MaxRegionSize << "\n"
                  << " bases, which is not officially supported at this time."<< std::endl;
        DiscovarUnhappy();
    }
}



void CheckBAMSize(size_t size, const String& REGIONS){
    if( size > DiscovarTools::MaxBAMSize && (REGIONS=="" || REGIONS=="all")){
        std::cout << "Using no REGIONS restriction when the BAM files total more than " << DiscovarTools::MaxBAMSize << "\n"
                  << "bytes, which is not officially supported at this time."<< std::endl;
        DiscovarUnhappy();
    }
}


}

void SetThreads( uint& NUM_THREADS, const Bool malloc_per_thread_check ) {
    if (malloc_per_thread_check)
    {    static bool malloc_warning = false;
         String mtt = Getenv( "MALLOC_PER_THREAD", "undefined" );
         if ( mtt != "1" && ! malloc_warning) {
	     // DISCOVAR MALLOC WARNING
	     std::cout << Date( ) << ": Warning: recommend doing "
	          << "'setenv MALLOC_PER_THREAD 1'\n"
	          << Date( ) << ": before Discovar, to improve "
	          << "computational performance." << std::endl; 
	     malloc_warning = true;
         }    }
    NUM_THREADS = configNumThreads(NUM_THREADS);
    omp_set_num_threads( NUM_THREADS );
}
