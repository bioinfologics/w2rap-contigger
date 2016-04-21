///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

//
//    Simulates reads from a reference.  No quals.

#ifndef READ_SIMULATOR_SIMPLE_CORE_H
#define READ_SIMULATOR_SIMPLE_CORE_H

#include "Basevector.h"
#include "CoreTools.h"
#include "FastaFileset.h"
#include "ReadError.h"
#include "random/RNGen.h"
#include "random/NormalRandom.h"
#include <iostream>
// #include <vector>

static inline String Tag(String S = "NS") { return Date() + " (" + S + "): "; } 

class RandomGen : public RNGen
{
public:
  RandomGen(const unsigned seed = 2672) : RNGen(seed) {}
  float float01() { return float(next()) / float(RNGen::RNGEN_RAND_MAX); }
  unsigned unsignedN(const size_t n) { return next() * n / RNGen::RNGEN_RAND_MAX; }
};


class RandomLength
{
  NormalRandom _normal;
  const unsigned _mu;
  const unsigned _sig;
  const unsigned _min;

public:
  RandomLength(const unsigned LEN, 
	       const unsigned LEN_SIG,
	       const unsigned LEN_MIN,
               const unsigned RANDOM_SEED) : 
    _normal(LEN, LEN_SIG),
    _mu(LEN),
    _sig(LEN_SIG),
    _min(LEN_MIN)
  {
    if (_sig == 0) ForceAssertLe(_min, LEN);
    _normal.seed(RANDOM_SEED);
  }

  unsigned value()
  {
    if (_sig == 0) return _mu;
    float len = 0;
    while (len < _min) len = _normal.value();
    return len;
  }

  unsigned mu() const { return _mu; };
  unsigned sig() const { return _sig; };
  unsigned min() const { return _min; };
};

void reference_generate_random(const size_t G_RANDOM,
			       BaseVecVec * bvv_p,
                               const unsigned RANDOM_SEED);

void reference_get(BaseVecVec * bvv_p,
                   const String & FASTB_REF,
                   const String & FASTA_REF,
                   const int REF_ID,
                   const int REF_START,
                   const int REF_STOP,
                   const size_t G_RANDOM,
                   const String & HEAD_RANDOM_REF_OUT,
                   const unsigned RANDOM_SEED);

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
				     RandomLength pair_len_rnd_gen = RandomLength(0,0,0,0U),
				     const size_t pad_bases = 0,	// increase the fragment by this length to support future deletions
				     const bool trimmable = true,	// should we allow trimming of the fragment if the non-circular contig runs out?
				     const bool paired = false);

void reads_sample_parallel(const BaseVecVec & bvv_ref, 
			   const vec<bool> & is_bv_circular,
			   const vec<size_t> & ibvs_ref,
			   const vec<size_t> & ibs_ref,
                           vec<Bool>& rc,
			   BaseVecVec * bvv_p,
                           const unsigned RANDOM_SEED,
                           const size_t i_thread, 
			   const size_t n_threads,
			   const size_t trim_size = 0);

void reads_perturb_parallel(BaseVecVec * bvv_p,
                            ReadErrorVecVec* pVREV,
			    const float ERR_DEL, 
			    const float ERR_INS, 
			    const float ERR_SUB,
                            const unsigned RANDOM_SEED,
			    const size_t i_thread, 
			    const size_t n_threads,
			    const size_t trim_size = 0);


#if 0		// not implemented?
void reads_sample_and_perturb_parallel( const BaseVecVec & bvv_ref,
    const vec<bool> & is_bv_circular,
    const vec<size_t> & ibvs_ref,
    const vec<size_t> & ibs_ref,
    vec<Bool>& rc,
    BaseVecVec * bvv_p,
    const bool FW_ONLY,
    ReadErrorVecVec* pVREV,
    const float ERR_DEL,
    const float ERR_INS,
    const float ERR_SUB,
    const unsigned RANDOM_SEED,
    const size_t i_thread,
    const size_t n_threads,
    bool paired = false );
#endif


#endif
