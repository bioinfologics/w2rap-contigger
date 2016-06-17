///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef NHOOD_INFO_STUFF
#define NHOOD_INFO_STUFF

#include "CoreTools.h"
#include "paths/HyperBasevector.h"
#include "paths/long/large/Lines.h"
#include "paths/long/large/tools/NhoodInfoState.h"

void TestDot( );

void CreateEdgeLabels( const HyperBasevectorX& hb, const vec<int>& inv,
                       const vec<vec<vec<vec<int>>>>& lines, const vec<int>& tol,
                       const vec<int>& llens, const vec<covcount>& cov, const vec<vec<covcount>>& covs,
                       vec< vec< std::pair<int,int> > >& hits, const vec<String>& subsam_names,
                       const vec<vec<int>>& count, const vec<String>& genome_names,
                       const vec<int>& used, const nhood_info_state& state, vec<String>& edge_names );

void ColorEdges( const HyperBasevectorX& hb, const vec<vec<int>>& count,
                 const vec<int>& seeds, const vec<Bool>& invisible,
                 const nhood_info_state& state, vec<String>& edge_color,
                 vec<int>& pen_widths );

Bool DefineSeeds( const HyperBasevectorX& hb, const vec<int>& inv,
                  const vec< triple<kmer<20>,int,int> >& kmers_plus,
                  const vec<vec<vec<vec<int>>>>& lines,
                  const vec<int>& tol, const vec<String>& genome_names,
                  const vec< std::pair<int,ho_interval> >& ambint, Bool& ambflag,
                  const vec< vec< std::pair<int,int> > >& hits, const nhood_info_state& state,
                  const int RANDOM_SEED, const String& SEEDS_MINUS, vec<int>& seeds,
                  const int max_seeds, std::ostream& tout );

void MakeDot( const HyperBasevectorX& hb, const vec<int>& inv,
              const vec<vec<vec<vec<int>>>>& lines, const vec<int>& tol,
              const vec<int>& llens, const vec<covcount>& cov, const vec<vec<covcount>>& covs,
              vec< vec< std::pair<int,int> > >& hits, const vec<String>& genome_names,
              const vec<int>& used, const nhood_info_state& state, const vec<int>& seeds,
              const vec<Bool>& invisible, const vec<String>& subsam_names,
              const vec<vec<int>>& count, const Bool LABEL_CONTIGS, const String& dot_file );

void MakeFasta( const HyperBasevectorX& hb, const vec<int>& used,
                const String& fasta_file );

#endif
