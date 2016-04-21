///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS
#include "CoreTools.h"
#include "VecUtilities.h"
#include "ParallelVecUtilities.h"
#include "lookup/LookAlign.h" 
#include "efasta/EfastaTools.h"
#include "paths/AssemblyCleanupTools.h"
#include "pairwise_aligners/SmithWatBandedA.h"

struct scompare : public std::binary_function<superb,superb,bool>
{
    bool operator()( const superb& s1, const superb& s2 ) const
    { return s1.FullLength() > s2.FullLength(); }
};


Assembly::Assembly( const String in_superb_file, const String in_contigs_file, const String in_scaff_graph_file ){
  // Loading scaffolds
  cout << Date( ) << ": loading superb file" << std::endl;
  ReadSuperbs( in_superb_file, scaffolds_ );
  
  // reading contig information

  if ( in_contigs_file.Contains(".efasta") ) {
    LoadEfastaIntoStrings(in_contigs_file, efastas_);
    fastgs_.resize( efastas_.size() );
    for ( size_t i = 0; i < efastas_.size(); i++ )
      fastgs_[i] = recfastg( ToString(i), basefastg( efastas_[i] ) );
  }else if ( in_contigs_file.Contains(".fastg") ) {
    LoadFastg(in_contigs_file, fastgs_);
    efastas_.resize( fastgs_.size() );
    efasta tmp;
    for ( size_t i = 0; i < fastgs_.size(); i++ )
    {
      tmp.clear();
      fastgs_[i].AsEfasta( tmp , fastg_meta::MODE_2);
      efastas_[i] = tmp;
    }
  }

  tigMap_.resize( efastas_.size() );
  for ( size_t tid = 0; tid != efastas_.size(); tid++ )
    tigMap_[tid] = ToString(tid);

  scaffMap_.resize( scaffolds_.size() );
  for ( int sid = 0; sid < scaffolds_.isize(); sid++ )
    scaffMap_[sid] = ToString(sid);

  if ( in_scaff_graph_file.nonempty() )
    BinaryReader::readFile( in_scaff_graph_file, &SG_ );
  else{
    SG_.Initialize( scaffolds_.isize() );
  }
}


Assembly::Assembly( const vec<superb>& scaffoldsIn, const VecEFasta& efastasIn, const vec<String>* scaffMap, const vec<String>* tigMap, const digraphE<sepdev>* SGin ){

  scaffolds_ = scaffoldsIn;
  efastas_   = efastasIn;
  fastgs_.resize( efastas_.size() );
  for ( size_t i = 0; i < efastas_.size(); i++ )
    fastgs_[i] = recfastg( ToString(i), basefastg( efastas_[i] ) );

  if ( tigMap == 0 ){
    tigMap_.resize( efastas_.size() );
    for ( size_t tid = 0; tid != efastas_.size(); tid++ )
      tigMap_[tid] = ToString(tid);
  }else tigMap_ = *tigMap;

  if ( scaffMap == 0 ){
    scaffMap_.resize( scaffolds_.size() );
    for ( int sid = 0; sid < scaffolds_.isize(); sid++ )
      scaffMap_[sid] = ToString(sid);
  }else scaffMap_ = *scaffMap;

  if ( SGin != 0 ){
    SG_ = *SGin;
  }else{
    SG_.Initialize( scaffolds_.isize() );
  }
}

Assembly::Assembly( const vec<superb>& scaffoldsIn, const vec<recfastg>& fastgsIn, const vec<String>* scaffMap, const vec<String>* tigMap, const digraphE<sepdev>* SGin ){

  scaffolds_ = scaffoldsIn;
  fastgs_    = fastgsIn;
  efastas_.resize( fastgs_.size() );
  efasta tmp;
  for ( size_t i = 0; i < fastgs_.size(); i++ )
  {
    tmp.clear();
    fastgs_[i].AsEfasta( tmp , fastg_meta::MODE_2);
    efastas_[i] = tmp;
  }
  if ( tigMap == 0 ){
    tigMap_.resize( efastas_.size() );
    for ( size_t tid = 0; tid != efastas_.size(); tid++ )
      tigMap_[tid] = ToString(tid);
  }else tigMap_ = *tigMap;

  if ( scaffMap == 0 ){
    scaffMap_.resize( scaffolds_.size() );
    for ( int sid = 0; sid < scaffolds_.isize(); sid++ )
      scaffMap_[sid] = ToString(sid);
  }else scaffMap_ = *scaffMap;

  if ( SGin != 0 ){
    SG_ = *SGin;
  }else{
    SG_.Initialize( scaffolds_.isize() );
  }
}



size_t Assembly::scaffoldsTotLen() const{
  size_t len = 0;
  for ( size_t is = 0; is < scaffolds_.size(); is++ )
    len += scaffolds_[is].FullLength();
  return len;
}

size_t Assembly::scaffoldsRedLen() const{
  size_t len = 0;
  for ( size_t is = 0; is < scaffolds_.size(); is++ )
    len += scaffolds_[is].ReducedLength();
  return len;
}

size_t Assembly::scaffoldsNtigs() const{
  size_t ntigs = 0;
  for ( size_t is = 0; is < scaffolds_.size(); is++ )
    ntigs += scaffolds_[is].Ntigs();
  return ntigs;
}

// check integrity of scafolds and contigs data: contig size in superb == contig size in efasta,
//  each contig used once and only once
void Assembly::check_integrity() const{
  
  cout << Date() << ": checking integrity" << std::endl;
  ForceAssertEq( efastas_.size(), fastgs_.size() );
  vec<int> tigs_used( efastas_.size(), 0);
  for ( size_t i = 0; i < efastas_.size( ); i++ ){    
    vec<String> s(1);
    s[0] = efastas_[i];
    ValidateEfastaRecord(s);
    ForceAssert( fastgs_[i].IsGapless() );
    ForceAssertEq( fastgs_[i].Length1(), efastas_[i].Length1() );
    ForceAssertEq( fastgs_[i].MinLength(fastg_meta::MODE_2), efastas_[i].MinLength() );
    ForceAssertEq( fastgs_[i].MaxLength(fastg_meta::MODE_2), efastas_[i].MaxLength() );
    efasta fe;
    fastgs_[i].AsEfasta( fe ,fastg_meta::MODE_2);
    basevector b1, b2;
    fe.FlattenTo( b1 );
    efastas_[i].FlattenTo( b2 );
    ForceAssertEq( b1, b2 );
    ForceAssertEq( fe.Length1(), efastas_[i].Length1() );
    ForceAssertEq( fe.MinLength(), efastas_[i].MinLength() );
    ForceAssertEq( fe.MaxLength(), efastas_[i].MaxLength() );
  }
  for ( size_t si = 0; si < scaffolds_.size(); si++ ){
    const superb & s = scaffolds_[si];
    ForceAssertGt( s.Ntigs(), 0 );
    for ( int tpos = 0; tpos < s.Ntigs(); tpos++ ){
      size_t tid = s.Tig(tpos);
      ForceAssertLt( tid, efastas_.size() );
      if ( efastas_[tid].Length1() != s.Len(tpos) ){
	DPRINT5( si, tpos, tid, s.Len(tpos), efastas_[tid].Length1() );
	ForceAssertEq( efastas_[tid].Length1(), s.Len(tpos) );
      }
      tigs_used[tid]++;
    }
  }
  vec<size_t> unused_tigs, overused_tigs;
  for ( size_t tid = 0; tid < tigs_used.size(); tid++ ){
    if ( tigs_used[tid] == 0 )
      unused_tigs.push_back( tid );
    else if ( tigs_used[tid] > 1 )
      overused_tigs.push_back(tid);
    
  }
  
  if ( unused_tigs.size() > 0 || overused_tigs.size() > 0 ){
    
    if ( unused_tigs.size() > 0 ){
      int max_un_size = efastas_.at( unused_tigs[0] ).Length1();
      for ( size_t i = 0; i < unused_tigs.size(); i++ )
	if (  efastas_.at( unused_tigs[i] ).Length1() > max_un_size )
	  max_un_size = efastas_.at( unused_tigs[i] ).Length1();

      cout << Date() << ": maximum size of unused contig is : " << max_un_size << std::endl;
    }

    DPRINT2( unused_tigs.size(), overused_tigs.size() );
    ForceAssert( unused_tigs.size() == 0 && overused_tigs.size() == 0 );
  }
  
  ForceAssertEq( scaffolds_.isize(), SG_.N() );

  return;
}


void Assembly::remove_small_scaffolds(const int min_scaffold_size) {
    cout << Date() << ": removing small scaffolds: " << std::endl;
    vec<int> verts_to_remove; 
    for ( int si = 0; si < scaffolds_.isize(); si++ )
	if ( scaffolds_.at(si).ReducedLength() < min_scaffold_size )
	    verts_to_remove.push_back(si);

    remove_scaffolds(verts_to_remove);

}

void Assembly::remove_scaffolds( const vec<Bool>& to_remove ) {
  cout << Date() << ": initial number of scaffolds = " << scaffolds_.size() << std::endl;
  ForceAssertEq( scaffolds_.size(), to_remove.size() );
  ForceAssertEq( efastas_.size(), fastgs_.size() );
  vec<int> verts_to_remove; 
  for ( size_t i = 0; i < to_remove.size(); i++ )
    if ( to_remove[i] ) {
      verts_to_remove.push_back(i);
      SG_.DeleteEdgesAtVertex( i );
    }
  
  EraseIf( scaffolds_, to_remove );
  EraseIf( scaffMap_, to_remove );
  SG_.RemoveEdgelessVertices( verts_to_remove );
  ForceAssertEq( scaffolds_.isize(), SG_.N() );
  ForceAssertEq( efastas_.size(), fastgs_.size() );
  cout << Date() << ": final number of scaffolds = " << scaffolds_.size() << std::endl;
}

void Assembly::remove_scaffolds( const vec<int>& s_to_remove ) {
  if ( s_to_remove.size() == 0 ) return;

  ForceAssertGe( Min(s_to_remove), 0 );
  ForceAssertLt( Max(s_to_remove), scaffolds_.isize() );

  vec<Bool> to_remove( scaffolds_.size(), False );
  for ( size_t i = 0; i < s_to_remove.size(); i++ )
    to_remove[ s_to_remove[i] ] = True;

  remove_scaffolds( to_remove );
}

void Assembly::remove_contigs( const vec<Bool>& to_remove )
{
  ForceAssertEq( efastas_.size(), fastgs_.size() );
  vec<int> offsets( efastas_.size(), 0 );

  int offset = 0;
  for ( size_t tid = 0; tid < efastas_.size(); tid++ ){
    if ( to_remove[tid] ){
      offsets[tid] = -1;
      offset++;
    }
    else{ offsets[tid] = offset; }
  }
    
  ForceAssertEq( offsets.size(), efastas_.size() );
  for ( int tid = 0; tid < offsets.isize(); tid++ ){
    if ( offsets[tid] > 0 ){
      efastas_[ tid - offsets[tid] ] = efastas_[tid];
      fastgs_[ tid - offsets[tid] ] = fastgs_[tid];
      tigMap_[tid - offsets[tid] ] = tigMap_[tid];
    }      
  }
  efastas_.resize( efastas_.size() - offset );
  fastgs_.resize( fastgs_.size() - offset );
  tigMap_.resize( tigMap_.size() - offset );
  ForceAssertEq( efastas_.size(), tigMap_.size() );
  
  for ( size_t si = 0; si < scaffolds_.size(); si++ ){
    for ( int tpos = 0; tpos < scaffolds_[si].Ntigs(); tpos++ ){
      int tid = scaffolds_[si].Tig(tpos);
      if ( offsets[tid] >= 0 ){
	int newtid = tid - offsets[tid];
	ForceAssertGe( newtid, 0 );
	scaffolds_[si].SetTig( tpos, newtid );
      }
      else{
	scaffolds_[si].RemoveTigByPos( tpos );
	tpos--;
      }
    }
  }
  ForceAssertEq( efastas_.size(), fastgs_.size() );
}

   
void Assembly::remove_small_contigs( const int min_size_solo, const int min_size_in ){
  cout << Date() << ": removing small contigs and renumbering" << std::endl;
  vec<int> offsets( efastas_.size(), 0 );
  vec<Bool> tigsToRemove( efastas_.size(), False);
  for ( size_t si = 0; si < scaffolds_.size(); si++ ){
    if ( scaffolds_[si].Ntigs() == 1 ){
      if ( scaffolds_[si].Len(0) < min_size_solo )
	tigsToRemove[ scaffolds_[si].Tig(0) ] = True;
    }
    else{
      for ( int tpos = 0; tpos < scaffolds_[si].Ntigs(); tpos++ ){
	if ( scaffolds_[si].Len(tpos) < min_size_in )
	  tigsToRemove[ scaffolds_[si].Tig(tpos) ] = True;
      }
    }
  }
  remove_contigs(tigsToRemove);
}

void Assembly::remove_unused_contigs(){
  cout << Date() << ": removing unused contigs and renumbering" << std::endl;
  ForceAssertEq( efastas_.size(), fastgs_.size() );
  //ForceAssertEq( scaffolds_.isize(), SG_.N() );
  vec<int> tigs_used( efastas_.size(), 0);
  for ( size_t si = 0; si < scaffolds_.size(); si++ ){
    const superb & s = scaffolds_[si];
    ForceAssertGt( s.Ntigs(), 0 );
    for ( int tpos = 0; tpos < s.Ntigs(); tpos++ ){
      size_t tid = s.Tig(tpos);
      ForceAssertLt( tid, efastas_.size() );
      if ( efastas_[tid].Length1() != s.Len(tpos) ){
	DPRINT5( si, tpos, tid, s.Len(tpos), efastas_[tid].Length1() );
	ForceAssertEq( efastas_[tid].Length1(), s.Len(tpos) );
      }
      tigs_used[tid]++;
    }
  }

  size_t unusedCt = 0;
  vec<int> offsets( efastas_.size(), 0 );
  int offset = 0;
  for ( size_t tid = 0; tid < efastas_.size(); tid++ ){
    if ( ! tigs_used[tid] ){
      offsets[tid] = -1;
      offset++;
      unusedCt++;
    }
    else{ offsets[tid] = offset; }
  }
  cout << Date( ) << ": found " << unusedCt << " unused contigs, removing" << std::endl;
  ForceAssertEq( offsets.size(), efastas_.size() );
  for ( int tid = 0; tid < offsets.isize(); tid++ ){
    if ( offsets[tid] > 0 ){
      efastas_[ tid - offsets[tid] ] = efastas_[tid];
      fastgs_[ tid - offsets[tid] ] = fastgs_[tid];
      tigMap_[tid - offsets[tid] ] = tigMap_[tid];
    }      
  }
  efastas_.resize( efastas_.size() - offset );
  fastgs_.resize( fastgs_.size() - offset );
  tigMap_.resize( tigMap_.size() - offset );
  ForceAssertEq( efastas_.size(), tigMap_.size() );
  ForceAssertEq( efastas_.size(), fastgs_.size() );

  cout << Date() << ": updating scaffolds tig ids" << std::endl;
  for ( size_t si = 0; si < scaffolds_.size(); si++ ){
    for ( int tpos = 0; tpos < scaffolds_[si].Ntigs(); tpos++ ){
      int tid = scaffolds_[si].Tig(tpos);
      if ( offsets[tid] >= 0 ){
	int newtid = tid - offsets[tid];
	ForceAssertGe( newtid, 0 );
	scaffolds_[si].SetTig( tpos, newtid );
      }
      else{
	scaffolds_[si].RemoveTigByPos( tpos );
	tpos--;
      }
    }
  }
  cout << Date() << ": done with removing unused contigs" << std::endl;
}


void Assembly::reorder(){
  // Sorting scaffolds according to size and renumbering contigs according
  //   to sequential appearance in scaffolds
  cout << Date() << ": sorting scaffolds" << std::endl;
  vec<int> neworder( scaffolds_.size(), vec<int>::IDENTITY );
  SortSync( scaffolds_, neworder, scompare() );
  PermuteVec( scaffMap_, neworder );
  SG_.ReorderVertices( neworder );
  renumber();
}



// renumber all the contigs sequentially according to the scaffold
void Assembly::renumber(){
  cout << Date() << ": renumbering contigs for ordered scaffolds" << std::endl;
  vec<String> otigMap;  
  vec<superb> oscaffolds = scaffolds_;
  VecEFasta oefastas;
  vec<recfastg> ofastgs;
  int cTid = -1;
  for ( size_t is = 0; is < scaffolds_.size(); is++ ){
    for ( int tpos = 0; tpos < scaffolds_[is].Ntigs(); tpos++ ){
      cTid++;
      int oTid = scaffolds_[is].Tig(tpos);
      oscaffolds.at(is).SetTig( tpos, cTid );
      oefastas.push_back( efastas_.at(oTid) );
      ofastgs.push_back( fastgs_.at(oTid) );
      otigMap.push_back( tigMap_.at(oTid) );
    }
  }
  efastas_.resize(0);
  fastgs_.resize(0);
  tigMap_.resize(0);
  size_t newScaffoldsTotLen = 0, newScaffoldsRedLen = 0;
  for ( size_t is = 0; is < scaffolds_.size(); is++ ){
    newScaffoldsTotLen += scaffolds_[is].FullLength();
    newScaffoldsRedLen += scaffolds_[is].ReducedLength();
  }
  scaffolds_ = oscaffolds; 
  oscaffolds.resize(0);
  efastas_ = oefastas;
  fastgs_ = ofastgs;
  oefastas.resize(0);
  ofastgs.resize(0);
  tigMap_ = otigMap;
}

struct sData{
  
  sData() : AmbEventCount_(0), MaxAmbCount_(0), RedLen_(0), MinRedLen_(0), MaxRedLen_(0), 
	    FullLen_(0), MinFullLen_(0), MaxFullLen_(0) {};

  int AmbEventCount_, MaxAmbCount_;
  int RedLen_, MinRedLen_, MaxRedLen_;
  int FullLen_, MinFullLen_, MaxFullLen_;
  basevector seq_, minseq_, maxseq_;
  basevector seq_rc_, minseq_rc_, maxseq_rc_;
};

void Assembly::dedup( const Bool Exact ) {
  // Remove duplicate scaffolds.
  cout << Date() << ": removing duplicate scaffolds: " << std::endl;
  cout << Date() << ": initial number of scaffolds = " << scaffolds_.size() << std::endl;

  // cononicalize efastas
  VecEFasta cefastas( efastas_.size() );
  VecEFasta cefastas_rc( efastas_.size() );
  for ( size_t ci = 0; ci < efastas_.size(); ci++ ){
    efastas_[ci].Canonicalize( cefastas[ci] );
    ValidateEfastaRecord( cefastas[ci] );
    efasta eseq_rc = efastas_[ci];
    eseq_rc.ReverseComplement();
    eseq_rc.Canonicalize( cefastas_rc[ci] );
    ValidateEfastaRecord( cefastas_rc[ci] );
  }

  // compute scaffold data

  vec< vec< vec<basevector> > > s2altTigs( scaffolds_.size() );
  vec<sData> s2data( scaffolds_.size() );

  for ( int si = 0; si < scaffolds_.isize(); si++ ){
    s2altTigs[si].resize( scaffolds_[si].Ntigs() );
    for ( int i = 0; i < scaffolds_[si].Ntigs(); i++ ){
      int ti = scaffolds_[si].Tig(i);
      efasta& tigi = cefastas[ ti ];
      s2data[si].AmbEventCount_ += tigi.AmbEventCount();
      s2data[si].MaxAmbCount_   += tigi.AmbCount();
      s2data[si].MinRedLen_     += tigi.MinLength();
      s2data[si].MaxRedLen_     += tigi.MaxLength();
      s2data[si].MinFullLen_    += tigi.MinLength();
      s2data[si].MaxFullLen_    += tigi.MaxLength();
      s2data[si].RedLen_        += tigi.Length1();
      s2data[si].FullLen_       += tigi.Length1();
      if ( i < scaffolds_[si].Ntigs() -1 ){
	s2data[si].MinFullLen_ += scaffolds_[si].Gap(i) - scaffolds_[si].Dev(i);
	s2data[si].MaxFullLen_ += scaffolds_[si].Gap(i) + scaffolds_[si].Dev(i);
	s2data[si].FullLen_    += scaffolds_[si].Gap(i);
      }

      fastavector fseq;
      tigi.FlattenMinTo( fseq );
      s2altTigs[si][i].push_back( fseq.ToBasevector() );
      if ( tigi.Ambiguities() > 0 ){
	tigi.FlattenTo( fseq );
	s2altTigs[si][i].push_back( fseq.ToBasevector() );
	tigi.FlattenMaxTo( fseq );
	s2altTigs[si][i].push_back( fseq.ToBasevector() );
      }else{
	s2altTigs[si][i].push_back( fseq.ToBasevector() );
	s2altTigs[si][i].push_back( fseq.ToBasevector() );
      }
    }
    DPRINT3( si, s2data[si].AmbEventCount_, s2data[si].MaxAmbCount_ );
  }
  
  for ( int si = 0; si < scaffolds_.isize(); si++ ) {
    int Ni = scaffolds_[si].Ntigs();
    basevector seq, minseq, maxseq;
    for ( int ci = 0; ci < Ni; ci++ ){
      minseq = Cat( minseq, s2altTigs[si][ci][0] );
      seq    = Cat( seq, s2altTigs[si][ci][1] );
      maxseq = Cat( maxseq, s2altTigs[si][ci][2] );
      if ( ci < Ni -1 )
	for ( int pi = 0; pi < scaffolds_[si].Gap(ci); pi++ ){
	  seq.push_back(0); minseq.push_back(0); maxseq.push_back(0);
	}
    }
    s2data[si].seq_ = seq; s2data[si].minseq_ = minseq; s2data[si].maxseq_ = maxseq; 
    basevector seq_rc, minseq_rc, maxseq_rc;
    for ( int i = 0; i < Ni; i++ ){
      int ci = Ni -1 -i;
      basevector ci_minseq_rc = s2altTigs[si][ci][0]; ci_minseq_rc.ReverseComplement();
      basevector ci_seq_rc = s2altTigs[si][ci][1];    ci_seq_rc.ReverseComplement();
      basevector ci_maxseq_rc = s2altTigs[si][ci][2]; ci_maxseq_rc.ReverseComplement();
      minseq_rc = Cat( minseq_rc, ci_minseq_rc );
      seq_rc    = Cat( seq_rc, ci_seq_rc );
      maxseq_rc = Cat( maxseq_rc, ci_maxseq_rc );
      if ( ci > 0 )
	for ( int pi = 0; pi < scaffolds_[si].Gap(ci-1); pi++ ){
	  seq_rc.push_back(0); minseq_rc.push_back(0); maxseq_rc.push_back(0);
	}
    }
    ForceAssertEq( minseq.size(), minseq_rc.size() );
    ForceAssertEq( seq.size(), seq_rc.size() );
    ForceAssertEq( maxseq.size(), maxseq_rc.size() );
    s2data[si].seq_rc_ = seq_rc; s2data[si].minseq_rc_ = minseq_rc; s2data[si].maxseq_rc_ = maxseq_rc; 
  }


  int removed_count = 0;
  int removed_size  = 0;
  vec <Bool> to_remove ( scaffolds_.size(), False );

  for ( int si = 0; si < scaffolds_.isize(); si++ ) {
    if ( to_remove[si] ) continue;
    int Ni = scaffolds_[si].Ntigs();
    for ( int sj = si + 1; sj < scaffolds_.isize(); sj++ ) {
      if ( to_remove[sj] ) continue;
      int Nj = scaffolds_[sj].Ntigs();
      ForceAssertEq( Nj, s2altTigs[sj].isize() );
      
      // roughly check if scaffolds could be duplicates of each other 
      if ( Ni != Nj ) 
	continue;
      if ( s2data[si].FullLen_ != s2data[sj].FullLen_ && 
	   ( s2data[si].MaxFullLen_ < s2data[sj].MinFullLen_ || 
	     s2data[sj].MaxFullLen_ < s2data[si].MinFullLen_ ) )
	continue;
      if ( s2data[si].RedLen_ != s2data[sj].RedLen_ && 
	   ( s2data[si].MaxRedLen_ < s2data[sj].MinRedLen_ || 
	     s2data[sj].MaxRedLen_ < s2data[si].MinRedLen_ ) )
	continue;
      
      cout << "\n\n\n";
      cout << Date() << ":  ---------     comparing scaffolds:" << std::endl;
      DPRINT2( si, sj );
      DPRINT2( s2data[si].AmbEventCount_, s2data[sj].AmbEventCount_ );
      DPRINT2( s2data[si].MaxAmbCount_, s2data[sj].MaxAmbCount_ );
      DPRINT2( s2data[si].seq_.size(), s2data[sj].seq_.size());

      // more detailed check      
      vec<int> areRcs( Ni, 0 );
      vec<int> areEqual( Ni, 0 );
      for ( int ci = 0; ci < Ni; ci++ ){
	
	if ( ci > 0 &&  areEqual[ci-1] != 1 &&  areRcs[ci-1] != 1 ) continue;
	
	int ti = scaffolds_[si].Tig(ci);
	vec< vec<String> > iblocks;
	cefastas[ti].GetBlocks( iblocks );

	if ( ci == 0 || areEqual[ci-1] == 1 ){
	  // checking duplication
	  int cje = ci;
	  int tje = scaffolds_[sj].Tig(cje);
	  if ( cefastas[ti] == cefastas[tje] ){
	    areEqual[ci] = 1;
	  }else if ( cefastas[ti].AmbEventCount() == cefastas[tje].AmbEventCount() ){
	    vec< vec<String> > eblocks;
	    cefastas[tje].GetBlocks( eblocks );
	    Bool pathEqual = True;
	    for ( int ib = 0; ib < iblocks.isize(); ib++ ){
	      if ( iblocks[ib].size() == 1 && eblocks[ib].size() == 1 ){
		if ( iblocks[ib][0] != eblocks[ib][0] ){
		  pathEqual =False;
		  break;
		}
	      }else{
		Bool foundEqualRoute = False;
		for ( int ir1 = 0; ir1 < iblocks[ib].isize(); ir1++ )
		  for ( int ir2 = 0; ir2 < eblocks[ib].isize(); ir2++ )
		    if ( iblocks[ib][ir1] == eblocks[ib][ir2] ){
		      foundEqualRoute = True;
		      break;
		    }
		if ( ! foundEqualRoute ){
		  pathEqual = False;
		  break;
		}
	      }
	    }
	    if ( pathEqual ) areEqual[ci] = 1;
	  }
	}
	
	if ( ci == 0 || areRcs[ci-1] == 1 ){
	  // checking reverse duplication	  
	  int cjr = Nj -ci -1;
	  int tjr = scaffolds_[sj].Tig(cjr);
	  efasta eseq_rc = cefastas_rc[tjr];
	  if ( cefastas[ti] == eseq_rc ){
	    areRcs[ci] = 1;
	  }else if ( cefastas[ti].AmbEventCount() == eseq_rc.AmbEventCount() ){
	    vec< vec<String> > rblocks;
	    eseq_rc.GetBlocks( rblocks );
	    Bool pathEqual = True;
	    for ( int ib = 0; ib < iblocks.isize(); ib++ ){
	      if ( iblocks[ib].size() == 1 && rblocks[ib].size() == 1 ){
		if ( iblocks[ib][0] != rblocks[ib][0] ){
		  pathEqual =False;
		  break;
		}
	      }else{
		Bool foundEqualRoute = False;
		for ( int ir1 = 0; ir1 < iblocks[ib].isize(); ir1++ )
		  for ( int ir2 = 0; ir2 < rblocks[ib].isize(); ir2++ )
		    if ( iblocks[ib][ir1] == rblocks[ib][ir2] ){
		      foundEqualRoute = True;
		      break;
		    }
		if ( ! foundEqualRoute ){
		  pathEqual = False;
		  break;
		}
	      }
	    }
	    if ( pathEqual ) areRcs[ci] = 1;
	  }
	}
      }    
      
      if ( Sum( areRcs ) > 0 || Sum( areEqual ) > 0 )
	DPRINT4( areRcs.size(), Sum( areRcs ), areEqual.size(), Sum( areEqual ) );

      if ( Sum( areRcs ) == areRcs.isize() || Sum( areEqual ) == areEqual.isize() ){

	String type = Sum( areRcs ) == areRcs.isize() ? "reverse complement" : "duplicate";
	
	cout << Date() << ": scaffold " << sj << "(l=" << scaffolds_[si].FullLength() << ") is " << type << " of " << si
	     << " (length " << scaffolds_[si].FullLength() << ")" << std::endl;
	
	to_remove[sj] = True;
	removed_count++;
	removed_size += scaffolds_[sj].ReducedLength();
	
      }else if ( ! Exact) { 
	cout << "\n";
	cout << Date() << ": ------------- RUNNING BANDED SMITH-WATERMAN -------------------\n\n";
	int amb_diff_count = Max( s2data[si].maxseq_.isize() - s2data[si].minseq_.isize(), 
				  s2data[sj].maxseq_.isize() - s2data[sj].minseq_.isize() );
	int amb_event_count = Max( s2data[si].AmbEventCount_, s2data[sj].AmbEventCount_ );
	int amb_max_count   = Max( s2data[si].MaxAmbCount_, s2data[sj].MaxAmbCount_ );
	DPRINT3( amb_diff_count, amb_event_count, amb_max_count );
	// no exact match, check alignment
	for ( int iter = 0; iter < 3; iter++ ){
	  if ( to_remove[sj] ) continue;
	  basevector si_seq;  
	  if ( iter == 0 )      si_seq = s2data[si].seq_;
	  else if ( iter == 1 ) si_seq = s2data[si].minseq_;
	  else if ( iter == 2 ) si_seq = s2data[si].maxseq_;
	  basevector sj_seq;
	  if ( iter == 0 )      sj_seq = s2data[sj].seq_;
	  else if ( iter == 1 ) sj_seq = s2data[sj].minseq_;
	  else if ( iter == 2 ) sj_seq = s2data[sj].maxseq_;

	  // flank sequences with the same ends for the case where efasta amibugity is at the edge
	  Bool FlankSequences = False;
	  if ( efastas_[ scaffolds_[si].Tig(0) ].front() == '{' || efastas_[ scaffolds_[sj].Tig(0) ].front() == '{' ||
	       efastas_[ scaffolds_[si].Tig(0) ].back() == '}' || efastas_[ scaffolds_[sj].Tig(0) ].back() == '}' )
	    FlankSequences = True;

	  basevector flank;
	  flank.SetFromString("AAAAAAAAAA");
	  if ( FlankSequences ){
	    cout << Date()  << " Flanking" << std::endl;
	    si_seq = Cat( flank, si_seq, flank );
	    sj_seq = Cat( flank, sj_seq, flank );
	  }

	  align aF;
	  int min_overlap = Min( si_seq.size(), sj_seq.size() );
	  DPRINT( min_overlap );
	  int errorsF = 0;
	  int lenDiff = abs( si_seq.isize() - sj_seq.isize() );
	  int bandwidth = 2 * lenDiff; 
	  int maxErr = amb_max_count + amb_event_count + lenDiff + 1 + round( 0.0001 * (double)si_seq.isize() );
	  DPRINT( maxErr );
	  bandwidth = Max( bandwidth, 10 );
	  DPRINT( bandwidth );
	  cout << Date() << ": aligning forward" << std::endl;
	  float scoreF = 
	    SmithWatBandedA2<unsigned short>(si_seq, sj_seq, 0, bandwidth, aF, errorsF );  
	  look_align laF;
	  laF.ResetFromAlign( aF, si_seq, sj_seq );
	  DPRINT5( scoreF, laF.a.extent1(), laF.a.extent2(), errorsF, maxErr ); 
	  if ( ( laF.a.extent1() >= min_overlap || laF.a.extent2() >= min_overlap ) && 
	       errorsF <= maxErr && laF.Fw1() ){
	    cout << Date() << ": scaffold " << sj << "(l=" << scaffolds_[si].FullLength() << ") is a copy of " << si
		 << " (length " << scaffolds_[si].FullLength() << ")" << std::endl;
	    
	    to_remove[sj] = True;
	    removed_count++;
	    removed_size += scaffolds_[sj].ReducedLength();
	  }else{
	    basevector sj_seq_rc;
	    if ( iter == 0 )      sj_seq_rc = s2data[sj].seq_rc_;
	    else if ( iter == 1 ) sj_seq_rc = s2data[sj].minseq_rc_;
	    else if ( iter == 2 ) sj_seq_rc = s2data[sj].maxseq_rc_;
	    
	    if ( FlankSequences )
	      sj_seq_rc = Cat( flank, si_seq, flank );

	    align aR;
	    int errorsR = 0;
	    cout << Date() << ": aligning reverse" << std::endl;
	    int scoreR = 
	      SmithWatBandedA2<unsigned short>(si_seq, sj_seq_rc, 0, bandwidth, aR, errorsR );  
	    look_align laR;
	    laR.ResetFromAlign( aR, si_seq, sj_seq_rc );
	    DPRINT5( scoreR, laR.a.extent1(), laR.a.extent2(), errorsR, maxErr ); 
	    if ( ( laR.a.extent1() >= min_overlap || laR.a.extent2() >= min_overlap ) && 
		 errorsR <= maxErr && laR.Fw1() ){
	      cout << Date() << ": scaffold " << sj << "(l=" << scaffolds_[si].FullLength() << ") is reverse-complement of " << si
		   << " (length " << scaffolds_[si].FullLength() << ")" << std::endl;
	      
	      to_remove[sj] = True;
	      removed_count++;
	      removed_size += scaffolds_[sj].ReducedLength();
	    }
	  }
	  cout << Date() << ":  vvvvvvvvvvvvvv  DONE WIHT ALIGNMENT vvvvvvvvvvvvvvv\n\n\n" << std::endl;
	}
      }
    }
  }

  EraseIf( scaffolds_, to_remove );
  EraseIf( scaffMap_, to_remove );
  EraseIf( s2data, to_remove );
  EraseIf( s2altTigs, to_remove );

  vec<int> verts_to_remove;  
  for ( int i = 0; i < to_remove.isize(); i++ )
    if ( to_remove[i] ){
      SG_.DeleteEdgesAtVertex( i );
      verts_to_remove.push_back( i );
    }

  SG_.RemoveEdgelessVertices( verts_to_remove );
  

  cout << Date() << ": removed " << removed_count << " duplicate scaffolds"
       << " (" << removed_size << " bases)" << std::endl;
  cout << Date() << ": final number of scaffolds = " << scaffolds_.size() << std::endl;
  remove_unused_contigs();
}


void Assembly::dedup2() {
  // Remove scaffolds that are possibly duplicate
  // The criteria are:
  // 1. two scaffolds have same number of contigs (n_contig)
  // 2. n_contig >= 2
  // 3. Each contig and gap much match 
  //    - Gaps are matched whan [ gap_size +/- 3 * std ] overlap
  //    - Contigs are matched when 
  //      - perfect efasta match for contig size < 50,000
  //      - 1/100 mismatch kmer rate for contig size >= 50,000 (!!!!! this is only meant to be temporary fix to the problem)
  // 4. Both rc and fw duplicates are checked

  int VERBOSITY = 0;
  const int EfastaMatchSize = 50 * 1000; 
  const int MaxDev  = 3;
  const double MaxMismatchRate = 0.01;
  cout << Date() << ": " << "Remove possible duplicate scaffolds" << std::endl;
  cout << Date() << ": " << "initial number of scaffolds = " << scaffolds_.size() << std::endl;
  
  int removed_count = 0;
  int removed_size = 0;
  vec <Bool> to_remove ( scaffolds_.size(), False );
  for ( int si = 0; si < scaffolds_.isize(); si++ ) {
    if ( to_remove[si] ) continue;
    int Ni = scaffolds_[si].Ntigs();
    //if ( Ni < 2 ) continue;
    for ( int sj = si + 1; sj < scaffolds_.isize(); sj++ ) {
      if ( to_remove[sj] ) continue;
      int Nj = scaffolds_[sj].Ntigs();
      if ( Nj != Ni ) continue;
      //VERBOSITY = ( si == 1 && sj == 6 ? 1: 0 );
      // check scaffolds duplicate in two passes: fw in pass=0, rc in pass=1
      for( int pass = 0; pass < 2; pass++ ) {
	if ( VERBOSITY >= 1 )
	cout << Date() << ": " << "pass= " << pass << std::endl;
	if ( VERBOSITY >= 1 )
	  cout << Date() << ": " << "check gap " << std::endl ;
	// compare the gap of two scaffolds
	bool gap_match = true;
	for ( int igap = 0; igap < Ni-1; ++igap ) {
	  int gap1 = scaffolds_[si].Gap(igap);
	  int dev1 = scaffolds_[si].Dev(igap);
	  int gap2 = scaffolds_[sj].Gap( pass == 0 ? igap : Ni -2 - igap);
	  int dev2 = scaffolds_[sj].Dev( pass == 0 ? igap : Ni -2 - igap);
	  if ( IntervalOverlap( gap1 - MaxDev * dev1, gap1 + MaxDev * dev1 + 1,
		gap2 - MaxDev * dev2, gap2 + MaxDev * dev2 + 1) == 0 ) {
	    gap_match = false;
	    break;
	  }
	}
	if ( !gap_match) continue;
	if ( VERBOSITY >= 1 )
	  cout << Date() << ": " << "check contigs " << std::endl ;

	// check if the contigs matchs
	bool tig_match = true;
	if ( VERBOSITY >= 1 )
	  cout<< Date()  << " ncontigs= " << Ni << std::endl;
	for ( int itig = 0; itig < Ni; ++itig) {
	  if ( VERBOSITY >= 1 )
	    cout << Date() << ": check contig " << itig << ": "
	    << " vs " << scaffolds_[sj].Tig( pass == 0 ? itig :  Nj -1 -itig ) << std::endl;
	  efasta &tigi = efastas_[ scaffolds_[si].Tig(itig) ];
	  efasta &tigj = efastas_[ scaffolds_[sj].Tig( pass == 0 ? itig :  Nj -1 -itig ) ] ;

	  if ( VERBOSITY >= 1 )
	    cout << Date() << ": tigi.size()= " << tigi.size() << std::endl;
	  if ( VERBOSITY >= 1 )
	    cout << Date() << ": tigj.size()= " << tigj.size() << std::endl;
	  
	  // check if contig sizes match
	  int tig_size = ( tigi.size() + tigj.size()  ) /2;
	  int MaxMismatch =  tig_size * MaxMismatchRate ;
	  if ( abs( (int)tigi.size() - (int)tigj.size() ) > MaxMismatch ) { tig_match = false; break; }

	  // require perfect efasta match if contig size less than EfastaMatchSize
	  if ( tig_size < EfastaMatchSize ) {
	    if ( VERBOSITY >= 1 )
	      cout << Date() << ":  check efasta " << std::endl;
	    vec<basevector> v_ibases, v_jbases;
	    tigi.ExpandTo(v_ibases);
	    tigj.ExpandTo(v_jbases);
	    if ( pass == 1 ) 
	      for ( size_t k = 0; k < v_jbases.size(); ++k ) { v_jbases[k].ReverseComplement(); }
	    bool foundEqual = False;
	    for ( size_t vi = 0; vi < v_ibases.size() && ! foundEqual; vi++ ){
	      for ( size_t vj = 0; vj < v_jbases.size(); vj++ ){
		if ( v_jbases[vj].size() != v_ibases[vi].size() )
		  continue;
		basevector jbases = v_jbases[vj];
		if ( jbases == v_ibases[vi] ){
		  foundEqual = True;
		  break;
		}
	      }
	    }
	    if ( ! foundEqual ) {
	      tig_match = false;
	      break;
	    }
	  } 
	  // larger contig size. do kmer matching
	  else {
	    if ( VERBOSITY >= 1 )
	      cout << Date() << ":  check kmers " << std::endl;
	    basevector base1, base2;
	    tigi.FlattenTo( base1 );
	    tigj.FlattenTo( base2 );
	    if ( pass == 1 ) base2.ReverseComplement();
	    const int K = 24;
	    ForceAssertGt( base1.isize( ), K );
	    vec< basevector > kmers1( base1.isize( ) - K + 1);
            #pragma omp parallel for
	    for ( int jz = 0; jz <= base1.isize( ) - K; jz += 1000 ) 
	      for ( int j = jz; j <= Min( base1.isize( ) - K, jz + 1000 ); j++ ) 
		kmers1[j].SetToSubOf( base1, j, K ); 
	      ParallelUniqueSort(kmers1);    
	    ForceAssertGt( base2.isize( ), K );
	    vec< basevector > kmers2( base2.isize( ) - K + 1);
            #pragma omp parallel for
	    for ( int jz = 0; jz <= base2.isize( ) - K; jz += 1000 ) 
	      for ( int j = jz; j <= Min( base2.isize( ) - K, jz + 1000 ); j++ ) 
	        kmers2[j].SetToSubOf( base2, j, K ); 
	    ParallelUniqueSort(kmers2);    

	    // compare how many kmers are identical for the two sorted list
	    int nkmer1 = kmers1.size(), nkmer2 = kmers2.size();
	    int nkmer = (nkmer1 + nkmer2)/2;
	    int count = 0;
	    for( size_t i = 0, j = 0; i < kmers1.size() && j < kmers2.size(); ) {
	      if ( kmers1[i] > kmers2[j] ) j++;
	      else if ( kmers1[i] < kmers2[j] )  i++;
	      else count++, i++, j++;
	    }
	    if ( VERBOSITY >= 1 ) {
	      cout << Date() << ": nkmer= " << nkmer << std::endl;
	      cout << Date() << ": duplicate= " << count << std::endl;
	    }
	    if ( abs(nkmer - count) > int( nkmer * MaxMismatchRate) ) {
	      tig_match = false;
	      break;
	    }
	  }
	} // for itig
	if ( !tig_match ) continue;
	// ---------------------------------------------------------------------------
	// Now we concluded that the two scaffolds are duplicate. Remove the later one
	// --------------------------------------------------------------------------
	{
	  String type = pass == 0 ? "fw duplicate" : "rc duplicate";
	  cout << Date() << ": scaffold " << sj << " is " << type << " of " << si
	    << " (length " << scaffolds_[si].FullLength() << ")" << std::endl;
	  to_remove[sj] = True;
	  removed_count++;
	  removed_size += scaffolds_[sj].ReducedLength();
	  break; // do not go second pass 
	}
      } // end pass 2
    } // end for sj
  } // end for si
  
  // now remove the duplicate scaffolds 
  {
    EraseIf( scaffolds_, to_remove );
    EraseIf( scaffMap_, to_remove );
    vec<int> verts_to_remove;  
    for ( int i = 0; i < to_remove.isize(); i++ )
      if ( to_remove[i] ){
	SG_.DeleteEdgesAtVertex( to_remove[i] );
	verts_to_remove.push_back( i );
      }
 
    SG_.RemoveEdgelessVertices( verts_to_remove );
  }

  cout << Date() << ": removed " << removed_count << " duplicate scaffolds"
       << " (" << removed_size << " bases)" << std::endl;
  cout << Date() << ": final number of scaffolds = " << scaffolds_.size() << std::endl;
  remove_unused_contigs();
}

void Assembly::dedup_exact() {
  // Remove duplicate scaffolds...currently only handles case of singleton contigs
  // which are duplicates fw or rc. --bruce 8 Jun 2011
  cout << Date() << ": removing duplicate scaffolds: " << std::endl;
  cout << Date() << ": initial number of scaffolds = " << scaffolds_.size() << std::endl;
  int removed_count = 0;
  int removed_size = 0;
  vec<int> remaining( scaffolds_.size(), vec<int>::IDENTITY );
  vec<int> verts_to_remove;
  for ( int si = 0; si < scaffolds_.isize(); si++ ) {
    if (scaffolds_[si].Ntigs() != 1) continue;
    for ( int sj = si + 1; sj < scaffolds_.isize(); sj++ ) {
      if (scaffolds_[sj].Ntigs() != 1) continue;
      efasta &tigi = efastas_[scaffolds_[si].Tig(0)];
      efasta &tigj = efastas_[scaffolds_[sj].Tig(0)];
      if (tigi.size() != tigj.size()) continue;
      basevector ibases, jbases;
      tigi.FlattenTo(ibases);
      tigj.FlattenTo(jbases);
      bool rc = False;
      if (ibases != jbases) {
	jbases.ReverseComplement();
	rc = True;
	if (ibases != jbases) continue;
      }
      /*
      cout << "scaffold " << sj << " duplicate of " << si
	   << " (length " << scaffolds_[si].FullLength()
	   << (rc ? " rc" : "") << ")" << std::endl;
      */
      scaffolds_.erase( scaffolds_.begin() + sj );
      scaffMap_.erase( scaffMap_.begin() + sj );
      verts_to_remove.push_back( remaining[sj] );
      SG_.DeleteEdgesAtVertex( remaining[sj] );
      remaining.erase( remaining.begin() + sj );
      sj--;
      removed_count++;
      removed_size += ibases.size();
    }
  }

  SG_.RemoveEdgelessVertices( verts_to_remove );

  cout << Date() << ": removed " << removed_count << " duplicate scaffolds"
       << " (" << removed_size << " bases)" << std::endl;
  cout << Date() << ": final number of scaffolds = " << scaffolds_.size() << std::endl;
  remove_unused_contigs();
}


void Assembly::set_min_gaps( const int min_gap ){
  cout << Date() << ": resetting gaps < " << min_gap << std::endl;
  for ( size_t is = 0; is < scaffolds_.size(); is++ )
    for ( int tpos = 0; tpos < scaffolds_[is].Ntigs() -1; tpos++ )
      if ( scaffolds_[is].Gap(tpos) < min_gap ) {
	scaffolds_[is].SetGap( tpos, min_gap );
	scaffolds_[is].SetDev( tpos, 0 );
      }
}


void Assembly::WriteFastg( const String head_out ) const {
//   fastg_meta FGM;
   Ofstream(out_g, head_out + ".assembly.fastg");
//   out_g << FGM.GetFileHeader( head_out ) << "\n";
   out_g << fastg_meta::GetFileHeader( head_out ) << "\n";
   for (size_t is = 0; is < scaffolds_.size(); is++) {
     const superb& S = scaffolds_[is];
     headfastg hfg( ToString(is), SG_.From(is) );
     basefastg bfg;
     for ( int it = 0; it < S.Ntigs(); it++) {
       int tigId = S.Tig( it );
       bfg += basefastg( efastas_[tigId] );
       if ( it < S.Ntigs() -1 )
	 bfg += basefastg( S.Gap(it), S.Dev(it) );
     } 
     recfastg rfg( hfg, bfg );
     rfg.Print( out_g );
   }
//   out_g << FGM.GetFileFooter() << "\n";
   out_g << fastg_meta::GetFileFooter() << "\n";
}

void Assembly::Write( const String head_out ) const {

  // writing output
  cout << Date() << ": writing output files" << std::endl;
  WriteSuperbs( head_out + ".superb", scaffolds_ );
  WriteSummary( head_out + ".summary", scaffolds_ );


  
  Ofstream( efout, head_out + ".contigs.efasta" );
  for ( size_t id = 0; id < efastas_.size(); id++ )
    efastas_[id].Print(efout, "contig_" + ToString(id) );
}

void Assembly::WriteExtra( const String head_out ) const{

  vec<fastavector> fastas(efastas_.size());
  for ( size_t id = 0; id < efastas_.size(); id++ )
    efastas_[id].FlattenTo( fastas[id] );
  Ofstream( fout, head_out + ".contigs.fasta" );
  for ( size_t id = 0; id < fastas.size(); id++ )
    fastas[id].Print(fout, "contig_" + ToString(id) );

  {
    vecfastavector vec_tigs;
    for ( size_t i = 0; i < fastas.size( ); i++ )
      vec_tigs.push_back_reserve( fastas[i] );
    vec_tigs.WriteAll( head_out + ".contigs.vecfasta" ); 
  }

  vecbasevector bases( efastas_.size() );
  for ( size_t id = 0; id < efastas_.size(); id++ )
    efastas_[id].FlattenTo( bases[id] );
  
  bases.WriteAll( head_out + ".contigs.fastb" );
  
  Ofstream( cmout, head_out + ".contigs.mapping" );
  for ( size_t id = 0; id < tigMap_.size(); id++ )
    cmout << ToString(id) + " from " + ToString( tigMap_[id] ) << "\n";

  Ofstream( smout, head_out + ".superb.mapping" );
  for ( size_t is = 0; is < scaffMap_.size(); is++ )
    smout << ToString(is) + " from " + ToString( scaffMap_[is] ) << "\n";

  WriteScaffoldedEFasta( head_out + ".assembly.efasta", efastas_, scaffolds_ );
  WriteScaffoldedFasta( head_out + ".assembly.fasta", fastas, scaffolds_ );

  BinaryWriter::writeFile( head_out + ".scaffold_graph.", SG_ );

}
