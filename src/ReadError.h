///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file ReadError.h
 * \author tsharpe
 * \date Oct 18, 2012
 *
 * \brief
 */
#ifndef READERROR_H_
#define READERROR_H_

#include "dna/Bases.h"
#include "feudal/MasterVec.h"
#include "feudal/SerfVec.h"
#include "dvString.h"
#include <ostream>

class ReadError
{
public:
    ReadError() : mLocation(0), mType(0), mRdBase(4), mRefBase(4) {}

    static unsigned char const GAP_CODE = 4;

    enum ErrType { SUBSTITUTION, INSERTION, DELETION };
    ReadError( unsigned location, ErrType type,
                unsigned char rdBase, unsigned char refBase )
    : mLocation(location), mType(type), mRdBase(rdBase), mRefBase(refBase) {}

    // compiler-supplied copying and destructor are OK

    unsigned getLocation() const { return mLocation; }
    ErrType getType() const { return ErrType(mType); }
    unsigned char getReadBase() const { return mRdBase; }
    unsigned char getRefBase() const { return mRefBase; }

    // modifier
    void setLocation(unsigned pos) { mLocation = pos; };

    // Add comparator so that the ReadError vector can be sorted and compared
    friend bool operator<( const ReadError& e1, const ReadError& e2) {
        if( e1.mLocation != e2.mLocation ) return e1.mLocation < e2.mLocation; 
        if( e1.mType != e2.mType ) return e1.mType < e2.mType;
        if ( e1.mType == DELETION ) return false;
        return e1.mRdBase < e2.mRdBase;
    }

    friend std::ostream& operator<<( std::ostream& os, ReadError const& err )
    { switch ( err.getType() )
      {
      case SUBSTITUTION:
          os << "S(" << Base::val2Char(err.getReadBase())
              << '/' << Base::val2Char(err.getRefBase());
          break;
      case INSERTION:
          os << "I(" << Base::val2Char(err.getReadBase());
          break;
      case DELETION:
          os << "D(" << Base::val2Char(err.getRefBase());
          break;
      }
      return os << ")@" << err.getLocation(); }

private:
    unsigned mLocation; // offset into the read
    unsigned char mType;
    unsigned char mRdBase; // the incorrect base present in the read
    unsigned char mRefBase; // the reference base that should've been there
    // above two member can contain 0=A, 1=C, 2=G, 3=T or 4=missing
    // mRdBase=4 for a deletion
    // mRefBase=4 for an insertion
};
TRIVIALLY_SERIALIZABLE(ReadError);

typedef SerfVec<ReadError> ReadErrorVec;
typedef MasterVec<ReadErrorVec> ReadErrorVecVec;
//extern template class SmallVec<ReadError,MempoolAllocator<ReadError> >;
//extern template class OuterVec<ReadErrorVec>;

#endif /* READERROR_H_ */
