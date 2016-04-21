///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
//
//  Author: Neil Weisenfeld - Oct 15, 2012 
//


#ifndef REFERENCEITERATOR_H_
#define REFERENCEITERATOR_H_

class ReferenceIterator {
public:
  ReferenceIterator(const BaseVec& ref,			// reference contig to walk
		    size_t start_pos,			// start position
		    bool rc = false,			// rc (returns complement and walks backwards)
		    bool circular = false,		// is contig circular?
		    double err_del=0.0,			// rate of deletions to generate
		    double err_ins=0.0,			// ...insertions
		    double err_sub=0.0,			// ...substitutions
		    RandomGen* err_rand=0)		// seed for RNG for errors
  : _ref(ref),
    _ref_size(ref.size()),
    _pos(start_pos),
    _start_pos(start_pos),
    _rc(rc),
    _inc(rc ? -1 : 1),
    _err_del(err_del),
    _err_ins(err_ins),
    _err_sub(err_sub),
    _err_rand(err_rand),
    _base(0U),
    _serial(0),
    _done(false),
    _our_err_rand(false) {

    // set last position based on direction and circularity
    if ( circular ) {
      _last_pos = (_pos + _ref_size - _inc ) % _ref_size ;	// go once around the horn; unsigned so add _ref_size first
    } else if ( _rc ) {
      _last_pos = 0;
    } else {
      _last_pos = _ref_size - 1;
    }

    if ( !_err_rand ) {
      _err_rand = new RandomGen();
      _our_err_rand = true;
    }

    make_base();
  };

  ~ReferenceIterator() {
    if (_our_err_rand)
      delete _err_rand;
  }



  unsigned char operator*() {
//    ForceAssertEq(_done, false);
    return out_base(_base);
  }

  ReferenceIterator& operator++() {		// prefix increment
    if ( !_done ) {
      make_base();
    }
    return *this;
  }


  bool done() const {
    return _done;			// we've exhausted the bases
  }

  bool last() const {
    return ( _pos == _last_pos );	// we're at the last base
  }

  const ReadErrorVec& getErrors() { return _errs;};


private:
  const BaseVec& _ref;
  const size_t _ref_size;


  size_t _pos;
  size_t _start_pos;


  const bool _rc;
  const int _inc;

  const double _err_del;
  const double _err_ins;
  const double _err_sub;

  RandomGen* _err_rand;

  unsigned char _base;

  size_t _serial;

  bool _done;

  size_t _last_pos;
  ReadErrorVec _errs;

  bool _our_err_rand;


  unsigned char out_base(unsigned char base) {
    if ( base > 4 )
      return base;
    else if ( _rc )
      return 3-base;
    else
      return base;
  }

  void increment(void) {
    ForceAssertEq(_done, false);
    if ( _pos == _last_pos )
      _done = true;
    else {
      _pos = ( _pos + _ref_size + _inc ) % _ref_size;
      _serial++;
    }
  }


  void make_base() {
    ForceAssertEq(_done, false);

    if ( _err_del + _err_ins + _err_sub == 0.0 ) {		// error free bases

      _base = _ref[_pos];
      increment();

    } else {

      float rand = _err_rand->float01();

      if ( rand < _err_ins ) { 					// insertion - generate base w/o increment
	_base = _err_rand->unsignedN(4);
	_errs.push_back( ReadError( _serial, ReadError::INSERTION, out_base(_base), ReadError::GAP_CODE ));

      } else if ( (rand -= _err_ins) < _err_sub ) { 		// substitution - change base and increment
	_base = (_ref[_pos] + 1 + _err_rand->unsignedN(3)) & 3u;
	_errs.push_back( ReadError( _serial, ReadError::SUBSTITUTION, out_base(_base), _ref[_pos] ));
	increment();

      } else if ( (rand -= _err_sub) < _err_del && !last()) { 	// deletion - increment, emit base, increment again -- really should loop again instead.
	_errs.push_back( ReadError( _serial, ReadError::DELETION, ReadError::GAP_CODE, out_base(_ref[_pos]) ) );
	increment();
	_base = _ref[_pos];
	increment();

      } else { // no error
	_base = _ref[_pos];
	increment();

      }

    }

  } // make_base()


};


#endif /* REFERENCEITERATOR_H_ */
