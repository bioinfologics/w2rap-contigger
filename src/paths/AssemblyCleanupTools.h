///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty ora guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef ASSEMBLY_CLEANUP_TOOLS_H
#define ASSEMBLY_CLEANUP_TOOLS_H

#include "CoreTools.h"
#include "efasta/EfastaTools.h"
#include "VecUtilities.h"
#include "paths/Sepdev.h"
#include "graph/Digraph.h"
#include "fastg/FastgTools.h"

/**
 * Struct: Assembly
 *
 * This is a class for manimpulating a superb and associate contigs file).
 */

class Assembly{
  
 private:
  
  vec<superb> scaffolds_;
  VecEFasta efastas_;
  vec<recfastg> fastgs_;
  digraphE<sepdev> SG_;
  

  vec<String> scaffMap_;
  vec<String> tigMap_;
  
 public:
  // -------------- Constructors

  Assembly( const vec<superb>& scaffolds, const VecEFasta& efastas, const vec<String>* scaffMap = 0, const vec<String>* tigMap = 0, const digraphE<sepdev>* SG = 0 );
  Assembly( const String scaffoldsFile, const String contigsFile, const String scaffoldGraphFile = "" );
  Assembly( const vec<superb>& scaffolds, const vec<recfastg>& fastgs, const vec<String>* scaffMap = 0, const vec<String>* tigMap = 0, const digraphE<sepdev>* SG = 0 );

  // --- useful functions
  void check_integrity() const;
  size_t scaffoldsTotLen() const;
  size_t scaffoldsRedLen() const;
  size_t scaffoldsNtigs() const;

  void remove_small_scaffolds( const int MIN_SCAFFOLD_SIZE );
  void remove_scaffolds( const vec<Bool>& );
  void remove_scaffolds( const vec<int>& );
  void remove_contigs( const vec<Bool>& to_remove );
  void remove_small_contigs( const int MIN_CONTIG_SIZE_SOLO, const int MIN_CONTIG_SIZE_IN );
  void remove_unused_contigs();
  void dedup(const Bool exact = True );
  void dedup2();
  void dedup_exact();
  void reorder();
  // renumber all the contigs sequentially according to the scaffold 
  void renumber();
  void set_min_gaps( const int min_gap );

  const vec<superb>& scaffolds() const { return scaffolds_; }
  const VecEFasta& efastas() const { return efastas_; }
  const digraphE<sepdev>& SG() const { return SG_; }
  const vec<String>& scaffMap() const { return scaffMap_; }
  const vec<String>& tigMap()  const { return tigMap_; }

  void Write( const String head_out ) const;
  void WriteExtra( const String head_out ) const;
  void WriteFastg( const String head_out ) const;
  void WriteAll( const String head_out ) const {
    Write( head_out );
    WriteExtra( head_out );
    WriteFastg( head_out );
  }

};

#endif
