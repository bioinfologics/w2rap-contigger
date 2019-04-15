///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Equiv.h"
#include "Vec.h"

// ===============================================================================

template<class INT>
equiv_rel_template<INT>::equiv_rel_template(INT n)
{    x_.resize(n);
     y_.resize(n);
     for ( INT i = 0; i < n; i++ )
          x_[i] = y_[i] = i;    }

template<class INT>
void equiv_rel_template<INT>::Initialize(INT n)
{    x_.resize(n);
     y_.resize(n);
     for ( INT i = 0; i < n; i++ )
          x_[i] = y_[i] = i;    }

template<class INT>
INT equiv_rel_template<INT>::Size() const
{    return x_.size();    }

template<class INT>
Bool equiv_rel_template<INT>::Equiv( INT i, INT j ) const
{    return y_[i] == y_[j];    }

template<class INT>
void equiv_rel_template<INT>::Next( INT& n ) const
{    n = x_[n];    }

template<class INT>
Bool equiv_rel_template<INT>::Join(INT a, INT b)
{    if ( y_[a] == y_[b] ) return false; // already equivalent
     std::swap( x_[a], x_[b]); 
     for( INT n = x_[a]; y_[n] != y_[a]; n = x_[n] )
          y_[n] = y_[a];
     return true;
}

template<class INT>
void equiv_rel_template<INT>::Orbit( INT a, vec<INT>& o ) const
{    o.resize(0);
     o.push_back(a);
     INT b = a;
     while(1)
     {    b = x_[b];
          if ( b == a ) break;
          o.push_back(b);    }    }



template<class INT>
INT equiv_rel_template<INT>::OrbitSize( INT a ) const
{    INT ans=1;
     INT b = a;
     while(1)
     {    b = x_[b];
          if ( b == a ) return ans;
          ans++;    }    }

template<class INT>
INT equiv_rel_template<INT>::OrbitCount( ) const
{    INT count = 0;
     for ( INT i = 0; i < (INT) x_.size( ); i++ )
          if ( x_[i] == y_[i] ) ++count;
     return count;    }

template<class INT>
void equiv_rel_template<INT>::OrbitReps( vec<INT>& reps ) const
{    reps.clear( );
     for ( INT i = 0; i < (INT) x_.size( ); i++ )
          if ( x_[i] == y_[i] ) reps.push_back(i);    }


template<class INT>
void equiv_rel_template<INT>::OrbitRepsAlt( vec<INT>& reps ) const
{    reps.clear( );
     for ( INT i = 0; i < (INT) x_.size( ); i++ )
          if ( i == y_[i] ) reps.push_back(i);    }



template<class INT>
void equiv_rel_template<INT>::Singletons( vec<INT>& sings ) const {
  sings.clear();
  for (INT i=0; i < (INT)x_.size(); i++)
    if( x_[i]==i )
      sings.push_back(i);
}

template<class INT>
bool equiv_rel_template<INT>::Singletons() const {
  for (INT i=0; i < (INT)x_.size(); i++)
    if( x_[i]==i )
      return true;
  return false;
}




// TEMPLATE INSTANTIATIONS
// These go hand in hand with the typedefs at the end of Equiv.h
template class equiv_rel_template<int>;
