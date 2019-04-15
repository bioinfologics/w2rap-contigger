///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef EQUIV_H
#define EQUIV_H

#include "Vec.h"
#include "system/TraceVal.h"


/**
   Class: equiv_rel

   Represents an equivalence relation on 1,...,n.

   Example of use: connected components of a graph -- see digraph::ComponentRelation().
*/
template<class INT>
class equiv_rel_template
{
  
public:
  
  equiv_rel_template( ) { };
  explicit equiv_rel_template(INT n);
  void Initialize(INT n);
  INT Size( ) const;
  
  void Next( INT& ) const;
  
  // MethodDecl : Equiv
  // Test whether two elements are in the same equivalence class.
  Bool Equiv(INT, INT) const;
  
  // MethodDecl: Join
  // Join the equivalence classes of the given elements.
  // Return False if the two elements are already in the same
  // equivalence class, True otherwise.
  Bool Join(INT, INT);
  
  // Method: ClassId
  // Return the representative of the given element's equivalence class.
  INT ClassId(INT i) const { return y_[i]; }
  
  // Method: Representative
  // Test whether the element is the representative of its equivalence class.
  Bool Representative( INT i ) const { return y_[i] == i; }
  
  // MethodDecl: OrbitCount
  // Return the number of equivalence classes.
  INT OrbitCount( ) const;
  
  // MethodDecl: Orbit
  // Find all elements equivalent to the given one.
  void Orbit( INT, vec<INT>& ) const;    // compute an orbit
  
  // MethodDecl: OrbitSize
  // Return the size of the given element's equivalence class.
  INT OrbitSize( INT ) const;            // compute an orbit's size 
  
  
  void Singletons( vec<INT>& ) const;    // list orbits of size 1
  bool Singletons() const;               // Do any exist?
  
  // OrbitReps and OrbitRepsAlt return orbit representatives, but in general
  // they return different representatives.
  //
  // In particular, OrbitRepsAlt returns a vector whose entries
  // are the possible return values of ClassId, in sorted order,
  // so eg you can BinPosition( orbit_reps_alt, ClassId(i) ).
  //
  // It's unfortunate that OrbitReps doesn't do this, but lots
  // of things already use these methods, and it seems like too
  // much work to figure out whether any of them depend on any
  // undocumented properties of OrbitReps.
  void OrbitReps( vec<INT>& reps ) const;
  void OrbitRepsAlt( vec<INT>& reps ) const;
  
private:
  
  vec<INT> x_, y_;
  
};



// TYPEDEFS
// These go hand in hand with the template instantiations at the end of Equiv.cc
typedef equiv_rel_template<int> equiv_rel;


#endif
