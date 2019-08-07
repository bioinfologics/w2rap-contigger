// Set.h defines some utilities often used with sets, e.g. Member()
// and i/o operators.

#ifndef SET_H
#define SET_H

#include "Vec.h"
#include <iostream>
#include <set>

template <class T, class C=std::less<T>>
using StdSet = std::set<T,C>;

template <class T, class C=std::less<T>>
using StdMultiset = std::multiset<T,C>;

template<class T,class C,class A> bool Member( std::set<T,C,A>& the_set, const T& value )
{    return the_set.find(value) != the_set.end( );    }

template<class T> std::ostream& operator<<(std::ostream& out, const std::set<T>& the_set)
{    out << the_set.size( ) << "\n";
     typename std::set<T>::const_iterator the_set_iter = the_set.begin(); 
     for ( ; the_set_iter != the_set.end(); ++the_set_iter )
       out << *the_set_iter << "\n";
     return out;    }

template<class T> std::istream& operator>>(std::istream& in, std::set<T>& s)
{    int n;
     in >> n;
     char c;
     in.get(c);
 
     T temp;
     for ( int idx = 0; idx < n; ++idx )
     {
       in >> temp;
       s.insert( temp );
     }
     return in;    }


/// Get all values in a set
template<class Value, typename Cmp> vec<Value> values(const std::set<Value,Cmp> &  m) {
  vec<Value> ret(m.size());
  int i = 0;
  for (typename std::set<Value,Cmp>::const_iterator iter = m.begin();
       iter != m.end();
       ++i, ++iter) {
    ret[i] = *iter;
  }
  return ret;
}

#endif
