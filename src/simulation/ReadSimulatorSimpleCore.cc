///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

//
//    Simulates reads from a reference.  No quals.
//

#include "Basevector.h"
#include "CoreTools.h"
#include "FastaFileset.h"
#include "feudal/BinaryStream.h"
#include "random/RNGen.h"
#include "random/NormalRandom.h"
#include "simulation/ReadSimulatorSimpleCore.h"
// #include <vector>

void reference_generate_random(const size_t G_RANDOM,
			       BaseVecVec * bvv_p,
                               const unsigned RANDOM_SEED) 
{
  bvv_p->resize(1);
  (*bvv_p)[0].resize(G_RANDOM);
  RandomGen random(RANDOM_SEED + 1);
  for (size_t i = 0; i < G_RANDOM; i++) 
    (*bvv_p)[0].set(i, random.unsignedN(4));
}

void reference_get(BaseVecVec * bvv_p,
                   const String & FASTB_REF,
                   const String & FASTA_REF,
                   const int REF_ID,
                   const int REF_START,
                   const int REF_STOP,
                   const size_t G_RANDOM,
                   const String & HEAD_RANDOM_REF_OUT,
                   const unsigned RANDOM_SEED) 
{                         
  if (FASTB_REF != "" || FASTA_REF != "") {
 
   if (FASTB_REF != "") {
      cout << Tag() << "Reading reference from '" << FASTB_REF << "'." << std::endl;
      bvv_p->ReadAll(FASTB_REF);
    }
   else {
      cout << Tag() << "Reading reference from '" << FASTA_REF << "'." << std::endl;
      vecString readNames;
      FastFetchReads(*bvv_p, &readNames, FASTA_REF);
    }
    
    if (REF_ID >= 0) {
      *bvv_p = BaseVecVec(1, (*bvv_p)[REF_ID]);
      
      if (REF_START >= 0 && REF_STOP >= 0) {
        (*bvv_p)[0].SetToSubOf((*bvv_p)[0], REF_START, REF_STOP - REF_START);
      }
    }
    
  }
  
  else {
    cout << Tag() << "Generating random reference of size " 
         << G_RANDOM << "." << std::endl;
    reference_generate_random(G_RANDOM, bvv_p, RANDOM_SEED);
    
    if (HEAD_RANDOM_REF_OUT != "") {
      cout << Tag() << "Writing random genome to '" << HEAD_RANDOM_REF_OUT << ".fastb'." << std::endl;
      bvv_p->WriteAll(HEAD_RANDOM_REF_OUT + ".fastb");
    }
  }
}   


static void make_empty_read( const size_t read_size, const size_t ibv_ref, const size_t ib_ref, bool rc,
	BaseVecVec *bvv_p,  vec<size_t> * ibvs_ref_p, vec<size_t> * ibs_ref_p, vec<Bool> * rc_p )
{
    bvv_p->push_back(BaseVec());	// TODO: investigate if there was a reason not to combine these two statements
    bvv_p->back().resize(read_size);

    ibvs_ref_p->push_back(ibv_ref);
    ibs_ref_p->push_back(ib_ref);

    rc_p->push_back( rc );
}

void reads_size_and_position_compute(const BaseVec & bv_ref, 
				     const bool is_bv_circular,
                                     const size_t ibv_ref,
				     const float coverage, 
				     RandomLength & len_rnd_gen,
                                     const unsigned RANDOM_SEED,
                                     const bool FW_ONLY,
				     const double prob_tile,
				     BaseVecVec * bvv_p, 
				     vec<size_t> * ibvs_ref_p, 
				     vec<size_t> * ibs_ref_p,
				     vec<Bool> * rc_p,
				     RandomLength pair_len_rnd_gen,
				     const size_t pad_bases,   // here we can pad out the read for an "extra" number of bases to support future deletions.
				     const bool trimmable,     // should we allow trimming of the fragment if the non-circular contig runs out?
				     const bool paired
				     )
{
  const size_t nb_ref = bv_ref.size();
  const size_t nb_reads_max = nb_ref * coverage;

  // ---- Generating random lengths
  //      Resizing base vectors now is much faster than 
  //      building a read base by base with push_back()

  const bool tiled = (prob_tile > 0.0);

  RandomGen random(RANDOM_SEED + 2); 
 
  size_t nb_reads = 0;

  size_t jb_ref = 0;   
  size_t ib0_ref = 0;   // the starting index in the reference base vec

  while (nb_reads < nb_reads_max) {
    
    if (nb_reads == 0 ||     // 1st time 
        jb_ref >= nb_ref) {  // reached the end
      jb_ref = 0;
      ib0_ref = (is_bv_circular ? random.float01() * nb_ref : 0);
    }


    size_t nb_read = len_rnd_gen.value();      // get a random read length
    size_t tile_skip = nb_read;			// by default we skip a full read if we're tiling

    // prob_tile == 1.0 => tiled reads
    // prob_tile == 0.5 => only half of tiled reads
    if (!tiled || random.float01() < prob_tile) { // include read 

      size_t ib_ref = 0;

      if (tiled) {
	ib_ref = (ib0_ref + jb_ref) % nb_ref;

	if (jb_ref + nb_read  > nb_ref) { // trim read length to fit?
	 // needs trimming
	  if ( trimmable )
	    nb_read = nb_ref - jb_ref;
	  else
	    nb_read = 0;	// if we can't trim, but need to, we'll skip it
	}
      }
      else { // not tiled 
        ib_ref = (is_bv_circular ? random.unsignedN(nb_ref) : random.unsignedN(nb_ref - nb_read));
        //cout << "ib_ref= " << ib_ref << std::endl;
      }

      if ( nb_read > 0 ) {
	if ( !paired ) {
	  make_empty_read( nb_read+pad_bases, ibv_ref, ib_ref, FW_ONLY ? false : random.unsignedN(2), bvv_p, ibvs_ref_p, ibs_ref_p, rc_p );

	  nb_reads += nb_read;		// for coverage purposes, how many bases did we cover -- don't include padding

	} else { // paired case
	  size_t len1 = pair_len_rnd_gen.value();
	  size_t len2 = pair_len_rnd_gen.value();

//	  std::cout << "len1+pad_bases=" << len1+pad_bases << ", len2+pad_bases=" << len2+pad_bases << std::endl;
	  if ( len1 <= nb_read && len2 <= nb_read ) {
	    if ( FW_ONLY || random.unsignedN(2) == 1U ) {
	      make_empty_read( len1+pad_bases, ibv_ref, ib_ref, false, bvv_p, ibvs_ref_p, ibs_ref_p, rc_p );
	      make_empty_read( len2+pad_bases, ibv_ref, ib_ref + nb_read - len2 - pad_bases , true,  bvv_p, ibvs_ref_p, ibs_ref_p, rc_p );
	    } else {
	      make_empty_read( len2+pad_bases, ibv_ref, ib_ref + nb_read - len2 - pad_bases , true,  bvv_p, ibvs_ref_p, ibs_ref_p, rc_p );
	      make_empty_read( len1+pad_bases, ibv_ref, ib_ref, false, bvv_p, ibvs_ref_p, ibs_ref_p, rc_p );
	    }
	    nb_reads += len1 + len2;	// for coverage purposes, how many bases did we cover -- don't include padding
	    tile_skip = len1;		// if we're tiling, skip the left read only(?)
	  }

	} // paired or unpaired

      } // if nb_read > 0
    } // if tiled read accepted

    jb_ref += tile_skip;
  }
}

void reads_sample_parallel(const BaseVecVec & bvv_ref, 
			   const vec<bool> & /* is_bv_circular */,
			   const vec<size_t> & ibvs_ref,
			   const vec<size_t> & ibs_ref,
                           vec<Bool>& rc,
			   BaseVecVec * bvv_p,
                           const unsigned RANDOM_SEED,
                           const size_t i_thread, 
			   const size_t n_threads,
			   const size_t trim_size		// generally not used, but *could* be used to trim bases that were padded to support deletions, but for which no
							        // deletions will ultimately be made.  Mostly here for completeness.
			   )
{
  RandomGen random(RANDOM_SEED + 3 + i_thread); // a generator for this thread, no need to lock. faster this way!

//  const size_t nbv_ref = bvv_ref.size();
  size_t nbv_tot = bvv_p->size();
  size_t ibv0, ibv1, nbv;

  ibv0 =  i_thread     * nbv_tot / n_threads;
  ibv1 = (i_thread + 1)* nbv_tot / n_threads;

  nbv = ibv1 - ibv0;

  for (size_t jbv = 0; jbv < nbv; jbv++) {
    // if (i_thread == 0) dots_pct(jbv, nbv);

    const size_t ibv = ibv0 + jbv;

    const BaseVec & bv_ref = bvv_ref[ibvs_ref[ibv]];
    const size_t nb_ref = bv_ref.size();
    
    BaseVec & bv = (*bvv_p)[ibv];
    const size_t nb = bv.size() - trim_size;
    if ( trim_size )
      bv.resize( nb);

    const size_t ib0 = ibs_ref[ibv];
    
    if ( !rc[ibv] ) {
      for ( size_t ib = 0; ib < nb; ib++ )
	bv.set(ib, bv_ref[(ib0+ib) % nb_ref]);
    } else {
      for (size_t ib = 0; ib < nb; ib++)
        bv.set(nb - ib - 1, 3u - bv_ref[(ib0 + ib) % nb_ref]);
    }

  }
}

void reads_perturb_parallel(BaseVecVec * bvv_p,
                            ReadErrorVecVec* pVREV,
			    const float ERR_DEL, 
			    const float ERR_INS, 
			    const float ERR_SUB,
                            const unsigned RANDOM_SEED,
			    const size_t i_thread, 
			    const size_t n_threads,
			    const size_t trim_size	// trim_size indicates how much padding was added in reads_size_and_position_compute to support deletions.  All reads will be trimmed.
			    )
{
  RandomGen random(RANDOM_SEED + 4 + i_thread); // a generator for this thread, no need to lock. faster this way!

  BaseVec bv_orig;
  ReadErrorVec errs;

  const size_t nbv_tot = bvv_p->size();
  const size_t ibv0 =  i_thread     * nbv_tot / n_threads;
  const size_t ibv1 = (i_thread + 1)* nbv_tot / n_threads;
  
  const size_t nbv = ibv1 - ibv0;


  for (size_t ibv = 0; ibv < nbv; ibv++) {
    BaseVec & bv_err = (*bvv_p)[ibv0 + ibv];
    bv_orig = bv_err;

    const size_t nb_orig = bv_orig.size();
    if ( trim_size > 0 )			// the final read will be shorter by trim_size bases.
      bv_err.resize(bv_orig.size() - trim_size);
    const size_t nb_err = bv_err.size();

    int tailroom = trim_size;
    size_t ib = 0;
    for (size_t jb = 0; ib < nb_err && jb < nb_orig; ib++) {
      float rand = random.float01();
      if (rand < ERR_INS) {  // insertion
	bv_err.set(ib, random.unsignedN(4));
	if ( pVREV )
	    errs.push_back(ReadError(ib,ReadError::INSERTION,bv_err[ib],ReadError::GAP_CODE));
      }
      else if ((rand -= ERR_INS) < ERR_SUB) { // substitution
	bv_err.set(ib, (bv_orig[jb] + 1 + random.unsignedN(3)) & 3u);
        if ( pVREV )
            errs.push_back(ReadError(ib,ReadError::SUBSTITUTION,bv_err[ib],bv_orig[jb]));
	jb++;
      }
      else if ((rand -= ERR_SUB) < ERR_DEL) { // deletion
	if (ib > 0  && ( trim_size == 0 || tailroom > 0) ) {
          if ( pVREV )
              errs.push_back(ReadError(ib,ReadError::DELETION,ReadError::GAP_CODE,bv_orig[jb]));
	  ib--;
	  jb++;
	  tailroom--;
	}
      }
      else {  // no errors
	bv_err.set(ib, bv_orig[jb]);
	jb++;
      }
    }
    if (ib < nb_err)
      bv_err.resize(ib); // reducing size is hopefully thread safe.

    // if (i_thread == 0) dots_pct(ibv, nbv);
    if ( pVREV )
    {
        (*pVREV)[ibv0+ibv] = errs;
        errs.clear();
    }
  }

}
