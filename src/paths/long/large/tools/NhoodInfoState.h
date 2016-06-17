///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef NHOOD_INFO_STATE_H
#define NHOOD_INFO_STATE_H

#include "CoreTools.h"
#include "kmers/KmerRecord.h"
#include "paths/HyperBasevector.h"
#include "paths/long/large/Lines.h"

class NhoodInfoEngine;

class nhood_info_state {

  public:

    String SEEDS;
    int DEPTH;
    Bool EXT;
    Bool GREEN;
    Bool BNG;
    String PURPLE;
    Bool NEATO;
    Bool REL;
    Bool LINENO;
    Bool SHOW_INV;
    Bool SHOW_CN;
    Bool SHOW_ALIGN;
    Bool COUNT;
    double FONTSIZE;
    double SCALE;
    Bool ALTREF;
    Bool COV2;
    double ASPECT;
    Bool SVG;
    vec<int> used;

    void Initialize( );
    Bool SetState( const String& line, std::ostream& out, const Bool SERVER,
                   const HyperBasevectorX& hb, const NhoodInfoEngine& engine );
    void SetUsed( const vec<int>& u ) {
        used = u;
    }

};

class NhoodInfoEngine {

  private:
    static const int L = 20;

    String dir;
    HyperBasevectorX hb;
    vec<int> inv;
    vec< triple<kmer<L>,int,int> > kmers_plus;
    vec<vec<vec<vec<int>>>> lines;
    vec<int> tol, npairs, llens;
    vec<String> genome_names, genome_names_alt;
    vec< std::pair<int,ho_interval> > ambint, ambint_alt;
    vec< vec< std::pair<int,int> > > hits, hits_alt;
    vec<covcount> cov; // for backward compatibility
    vec<vec<covcount>> covs;
    vec<String> subsam_names;
    vec<vec<int>> count;

    class nhood_info_state state;

  public:
    void Initialize(const String& DIR_IN, const bool SEQ_LOOKUP, const bool EXT,
                    const bool COV2);

    Bool HasCounts( ) const {
        return count.nonempty( );
    }
    const vec<vec<int>>& Count( ) const {
        return count;
    }
    Bool HasLines( ) const {
        return lines.nonempty( );
    }
    String Dir( ) const {
        return dir;
    }

    void SetState(const String& SEEDS,
                  const int DEPTH,
                  const bool EXT,
                  const bool COUNT,
                  const bool GREEN,
                  const bool BNG,
                  const String PURPLE,
                  const bool NEATO,
                  const bool REL,
                  const bool LINENO,
                  const bool SHOW_INV,
                  const bool SHOW_CN,
                  const bool SHOW_ALIGN,
                  const double FONTSIZE,
                  const double SCALE,
                  const bool ALTREF,
                  const bool COV2,
                  const double ASPECT,
                  const bool SVG );

    void RunAsClient(const Bool INTERACTIVE, const int DEPTH,
                     const String& OUT, const bool PNG, const bool PDF,
                     const Bool SVG, const int RANDOM_SEED,
                     const String& SEEDS_MINUS, const Bool LABEL_CONTIGS,
                     const String DOTEXTRA = "" );

    void RunAsServer(const String& SERVER_DIR, const Bool LABEL_CONTIGS );
};

#endif
