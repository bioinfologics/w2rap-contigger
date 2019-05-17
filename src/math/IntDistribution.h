///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
// author: Filipe Ribeiro      03/2011
// 
//
//


#ifndef _MATH__INT_DISTRIBUTION_H
#define _MATH__INT_DISTRIBUTION_H

#include "feudal/BinaryStream.h"
#include "feudal/CharString.h"
#include "Vec.h"
#include <deque>


#define __INT_FUNCTION_BINARY_VERSION__ 2


template<class T>
class IntFunction
{
    bool _neg_inf;
    bool _pos_inf;
    int _x0;
    std::deque<T> _v;

public:
    IntFunction(const int x0 = 0,
                const int x1 = 0,
                const T v = 0)
            : _neg_inf(false),
              _pos_inf(false),
              _x0(x0),
              _v((x1 >= x0) ? x1 - x0 + 1 : 0, v)
    {
        ForceAssertGe(x1, x0);
    }

    ~IntFunction()
    {
        // std::cout << "~int_func(): " << x_min() << ", " << x_max()
        //     << " _neg.sz(): " << _negative.size()
        //     << " _pos_sz(): " << _positive.size() << std::endl;
    }

    int size() const { return _v.size(); }

    int x_min() const { return _x0; }
    int x_max() const { return _x0 + size() - 1; }
    T f_x_min() const { return _v.front(); }
    T f_x_max() const { return _v.back(); }

    int x_f_min() const;
    int x_f_max() const;
    T f_min() const { return (*this)[x_f_min()]; }
    T f_max() const { return (*this)[x_f_max()]; }

    T sum_in(const int x0, const int x1) const
    {
        T sum = 0;
        for (int x = x0; x < x1; x++) sum += (*this)[x];
        return sum;
    }
    T sum_below(const int x1) const { return sum_in(x_min(), x1); }
    T sum_above(const int x0) const { return sum_in(x0, x_max()); }
    T sum() const { return sum_in(x_min(), x_max()); }


    void expand_neg_infinity() { _neg_inf = true; }
    void expand_pos_infinity() { _pos_inf = true; }
    void expand_infinity() { _neg_inf = _pos_inf = true; }

    T operator [] (const int x) const
    {
        if (x < x_min()) return (_neg_inf) ? f_x_min() : 0;
        if (x > x_max()) return (_pos_inf) ? f_x_max() : 0;
        return _v[x - _x0];
    }

    T & operator [] (const int x)
    {
        if (x > x_max()) {
            const T f1 = (_pos_inf) ? f_x_max() : 0;
            _v.insert(_v.end(), x - x_max(), f1);
        }
        if (x < x_min()) {
            const T f0 = (_neg_inf) ? f_x_min() : 0;
            const int n = x_min() - x;
            for (int i = 0; i != n; i++) _v.push_front(f0);
            //_v.insert(_v.rend(), x_min() - x, f0);
            _x0 = x;
        }
        return _v[x - _x0];
    }

    IntFunction & operator += (const IntFunction & f)
    {
        for (int x = f.x_min(); x <= f.x_max(); x++) (*this)[x] += f[x];
        return *this;
    }

    virtual void to_text_file(const String & fn) const;


public:
    // ---- SELF_SERIALIZABLE method
    void writeBinary(BinaryWriter& writer) const;

    // ---- SELF_SERIALIZABLE method
    void readBinary(BinaryReader& reader);

    // ---- SELF_SERIALIZABLE method
    static size_t externalSizeof() { return 0; }
};


SELF_SERIALIZABLE(IntFunction<int64_t>);
SELF_SERIALIZABLE(IntFunction<int32_t>);
SELF_SERIALIZABLE(IntFunction<int16_t>);
SELF_SERIALIZABLE(IntFunction<int8_t>);
SELF_SERIALIZABLE(IntFunction<uint64_t>);
SELF_SERIALIZABLE(IntFunction<uint32_t>);
SELF_SERIALIZABLE(IntFunction<uint16_t>);
SELF_SERIALIZABLE(IntFunction<uint8_t>);
SELF_SERIALIZABLE(IntFunction<float>);
SELF_SERIALIZABLE(IntFunction<double>);




template<class T>
int IntFunction<T>::x_f_max() const
{
    int x_f_max = x_min();
    T f_max = (*this)[x_f_max];

    const int x1 = x_max();
    for (int x = x_f_max + 1; x <= x1; x++) {
        const T f = (*this)[x];
        if (f > f_max) {
            x_f_max = x;
            f_max = f;
        }
    }
    return x_f_max;
}


template<class T>
int IntFunction<T>::x_f_min() const
{
    int x_f_min = x_min();
    T f_min = (*this)[x_f_min];

    const int x1 = x_max();
    for (int x = x_f_min + 1; x <= x1; x++) {
        const T f = (*this)[x];
        if (f < f_min) {
            x_f_min = x;
            f_min = f;
        }
    }
    return x_f_min;
}



template<class T>
void IntFunction<T>::to_text_file(const String & fn) const
{
    const int x0 = x_min();
    const int x1 = x_max();

    std::ofstream os;
    os.open(fn.c_str());
    os << "# x_min = " << x0 << std::endl;
    os << "# x_max = " << x1 << std::endl;
    os << "# 1:x  2:f_x" << std::endl;
    os << std::fixed;

    for (int x = x0; x <= x1; x++) {
        os << std::setw(10) << x << " "
           << std::setw(16) << (*this)[x]
           << std::endl;
    }
    os.close();
}




// ---- SELF_SERIALIZABLE method
template<class T>
void IntFunction<T>::writeBinary(BinaryWriter& writer) const
{
    const int version = __INT_FUNCTION_BINARY_VERSION__;
    writer.write(version);

    writer.write(_neg_inf);
    writer.write(_pos_inf);
    writer.write(_x0);
    const size_t n = _v.size();
    writer.write(n);
    writer.writeItr(_v.begin(), _v.end());
}


// ---- SELF_SERIALIZABLE method
template<class T>
void IntFunction<T>::readBinary(BinaryReader& reader)
{
    int version;
    reader.read(&version);
    if (version != __INT_FUNCTION_BINARY_VERSION__) {
        std::cout << Date() << "**** binary data on disk is from a different code version." << std::endl;
        ForceAssertEq(version, __INT_FUNCTION_BINARY_VERSION__);
    }

    reader.read(&_neg_inf);
    reader.read(&_pos_inf);
    reader.read(&_x0);
    size_t n;
    reader.read(&n);
    _v.resize(n, 0);
    reader.readItr(_v.begin(), _v.end());
}









template<class T>
class IntFunctionPrimitive
{
    const IntFunction<T> & _f;
    IntFunction<T>         _f_sum;

public:
    IntFunctionPrimitive(const IntFunction<T> & f)
            : _f(f),
              _f_sum(f.x_min(), f.x_max())
    {
        const int x0 = _f.x_min();
        const int x1 = _f.x_max();

        _f_sum[x0] = _f[x0];
        for (int x = x0 + 1; x <= x1; x++)
            _f_sum[x] = _f_sum[x - 1] + _f[x];
    }

    T f_sum(const int a, const int b) const;

};




template<class T>
T IntFunctionPrimitive<T>::f_sum(const int a, const int b) const
{
    if (a == b) return _f[a];

    const int x0 = _f.x_min();
    const int x1 = _f.x_max();
    const T f0 = _f.f_x_min();
    const T f1 = _f.f_x_max();


    const int aa = (a < b) ? a - 1 : b;
    const int bb = (a < b) ? b     : a - 1;

    const int n1a = (aa > x1) ? aa - x1 : 0;
    const int n1b = (bb > x1) ? bb - x1 : 0;

    const int n0a = (aa < x0) ? aa - x0 : 0;
    const int n0b = (bb < x0) ? bb - x0 : 0;

    const T sum_a = (aa >= x1) ? _f_sum.f_x_max() : (aa <= x0) ? _f_sum.f_x_min() : _f_sum[aa];
    const T sum_b = (bb >= x1) ? _f_sum.f_x_max() : (bb <= x0) ? _f_sum.f_x_min() : _f_sum[bb];


    return (sum_b - sum_a) + f1 * (n1b - n1a) + f0 * (n0b - n0a);
}

class IntFrequencies : public IntFunction<size_t>
{
public:
    size_t freq(const int x) const { return (*this)[x]; }

    void to_text_file(const String & head) const
    {
        const String fn = head + ".freq";

        std::ofstream os;
        os.open(fn.c_str());
        os << "# 1:x  2:freq(x) 3:cum(x) 4:freq_norm(x) 5:cum_norm(x)" << std::endl;

        const int x0 = x_min();
        const int x1 = x_max();

        os << std::fixed;
        size_t total = 0;
        size_t cum = 0;
        for (int x = x0; x <= x1; x++) cum += freq(x);
        double cum_double = cum;

        cum = 0;
        for (int x = x0; x != x1; x++) {
            const size_t fx = freq(x);
            cum += fx;
            os << std::setw(10) << x << " "
               << std::setw(16) << fx << " "
               << std::setw(16) << cum << " "
               << std::setw(16) << std::setprecision(12) << double(fx) / cum_double << " "
               << std::setw(16) << std::setprecision(12) << double(cum) / cum_double << std::endl;
        }

        os.close();
    }
};



SELF_SERIALIZABLE(IntFrequencies);





class IntDistribution 
{
public:
  typedef double real_t;

private: 
  IntFunction<real_t> _prob;
  IntFunction<real_t> _prob_le;
  IntFunction<real_t> _prob_le_sum;

  void _normalize();
public:
  
  IntDistribution()
    : _prob(),
      _prob_le(),
      _prob_le_sum()
  {}
  
  IntDistribution(const IntFunction<real_t> & func) 
    : _prob(func), 
      _prob_le(func.x_min(), func.x_max(), 0.0),
      _prob_le_sum(func.x_min(), func.x_max(), 0.0)
  {
    _normalize();
  }


  ~IntDistribution()
  {
    //cout << "~int_dist(): " << _prob.x_min() << ", " << _prob.x_max() 
    //     << " le: " << _prob_le.x_min() << ", " << _prob_le.x_max() 
  }

  operator bool() const { return (_prob.size() > 1); }
  
  static IntDistribution gaussian(const int mean, 
                                  const int sigma,
                                  const real_t n_sigma = 4.0);
  static IntDistribution uniform(const int a,
                                 const int b);
  
  template <class INT_t>
  void from_hits(const vec<INT_t> & hits)
  {
    IntFunction<real_t> f;
    for (size_t i = 0; i != hits.size(); i++)
      f[hits[i]]++;
    _prob = f;
    _normalize();
  }


  int x_min() const { return _prob.x_min(); }
  int x_max() const { return _prob.x_max(); }


  real_t prob(const int x) const 
  { 
    if (x < x_min() || x > x_max()) return 0.0;
    return _prob[x]; 
  }
  
  real_t prob_le(const int x) const 
  {
    if (x < x_min()) return 0.0;
    if (x > x_max()) return 1.0;
    return _prob_le[x]; 
  }

  real_t prob_lt     (const int x)                const { return prob_le(x - 1); }
  real_t prob_ge     (const int x)                const { return 1.0 - prob_lt(x); }
  real_t prob_gt     (const int x)                const { return 1.0 - prob_le(x); }
  real_t prob_in     (const int x0, const int x1) const { return prob_le(x1) - prob_lt(x0); }
  real_t prob_not_in (const int x0, const int x1) const { return 1.0 - prob_in(x0, x1); }
  real_t prob_max    ()                           const { return prob(x_prob_max()); }
  
  // x_prob_le(random01) is a simple way of getting a random sample from the distribution. 
  int    x_prob_le   (const real_t F)             const { return quantile(F); } 
  int    x_prob_max  ()                           const { return _prob.x_f_max(); } 

  real_t mean        ()                                     const;
  real_t moment      (const real_t c, const unsigned order) const;
  real_t variance    ()                                     const { return moment(mean(), 2); }
  int    quantile    (const real_t Fq)                      const;
  int    median      ()                                     const { return quantile(0.5); }
  int    mode        ()                                     const { return x_prob_max(); }

  real_t prob_in_sum(const int x0, const int x1, const int n) const;
  


  IntFunction<real_t> probs() const { return _prob; }

  void split(IntDistribution * p_dist_left,
             IntDistribution * p_dist_right,
             const int x_split = 0) const;

  IntDistribution reverse() const;
  

  void to_text_file(const String & head) const;


public:
  // ---- SELF_SERIALIZABLE method
  void writeBinary( BinaryWriter& writer ) const
  { _prob.writeBinary(writer); }

  // ---- SELF_SERIALIZABLE method
  void readBinary( BinaryReader& reader )
  { _prob.readBinary(reader); _normalize(); }

  // ---- SELF_SERIALIZABLE method
  static size_t externalSizeof() { return 0; }

};



SELF_SERIALIZABLE(IntDistribution);









// Compute the distribution of the sum of 2 random variables
//   Z = X + Y    =>   p(Z) = p(X) * p(Y)
IntDistribution distribution_of_sum(const IntDistribution & x_dist,
                                    const IntDistribution & y_dist);

// Compute the distribution of the sum of 2 random variables
//   Z = X - Y    =>   p(Z) = p(X) * p(Y)
IntDistribution distribution_of_difference(const IntDistribution & x_dist,
                                           const IntDistribution & y_dist);










class IntLogDistribution : public IntFunction<double>
{
public:
  typedef double real_t;

public:
  IntLogDistribution(const int x0 = 0, const int x1 = 0, const real_t val = 0) 
    : IntFunction<real_t>(x0, x1, val) 
  {
    this->expand_infinity();
  }

  IntLogDistribution(const IntDistribution & dist)
    : IntFunction<real_t>(dist.x_min(), dist.x_max())
  {
    const int x0 = dist.x_min();
    const int x1 = dist.x_max();
    
    real_t p_min = -1; // undefined
    for (int x = x0; x <= x1; x++) {
      const real_t p = dist.prob(x);
      ForceAssertGe(p, 0);
      if (p > 0.0 && (p_min < 0 || p < p_min))
        p_min = p;
    }

    IntLogDistribution & me = *this;
    for (int x = x0; x <= x1; x++) {
      const real_t p = dist.prob(x);
      me[x] = (p > 0.0) ? log(p) : log(p_min);
    }        
    me.expand_infinity();
  }


  IntLogDistribution(const IntFunction<real_t> & func)
    : IntFunction<real_t>(func.x_min(), func.x_max())
  {
    const int x0 = func.x_min();
    const int x1 = func.x_max();
    
    real_t f_min = -1; // undefined
    for (int x = x0; x <= x1; x++) {
      const real_t f = func[x];
      ForceAssertGe(f, 0);
      if (f > 0.0 && (f_min < 0 || f < f_min))
        f_min = f;
    }

    IntLogDistribution & me = *this;
    for (int x = x0; x <= x1; x++) {
      const real_t f = func[x];
      me[x] = (f > 0.0) ? log(f) : log(f_min);
    }        
    me.expand_infinity();
  }


  ~IntLogDistribution()
  {
    // std::cout << "~int_log_dist(): " << x_min() << ", " << x_max() << std::endl;
  }


  IntDistribution int_distribution() const 
  {
    const int x0 = x_min();
    const int x1 = x_max();
    IntFunction<real_t> func(x0, x1);
    const IntLogDistribution & me = *this;

    real_t f_max = me.f_max(); 

    for (int x = x0; x <= x1; x++)
      func[x] = exp(me[x] - f_max);   // = 1 at me[x] = f_max

    // here I rely on the IntDistribution constructor that accepts an IntFunction
    return func;
  }


  void to_text_file(const String & fn) const
  {
    const real_t fmax = f_max();
    
    const int x0 = x_min();
    const int x1 = x_max();

    std::ofstream os;
    os.open(fn.c_str());
    os << "# x_min = " << x0 << std::endl;
    os << "# x_max = " << x1 << std::endl;
    os << "# 1:x  2:f_x" << std::endl;
    os << std::fixed;
    
    for (int x = x0; x <= x1; x++) {
      os << std::setw(10) << x << " "
	 << std::setw(16) << ((*this)[x] - fmax)
	 << std::endl;
    }
    os.close();
  }




};





#endif
