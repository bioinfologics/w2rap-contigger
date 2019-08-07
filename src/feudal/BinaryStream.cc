/*
 * \file BinaryStream.cc
 * \author tsharpe
 * \date Aug 13, 2009
 *
 * \brief
 */
#include "feudal/BinaryStream.h"
#include "system/System.h"
#include "BinaryStream.h"

void BinaryReader::readLoop( char* buf, size_t len )
{
    size_t remain = 0;
    while ( len )
    {
        remain = fillBuf(BUF_SIZ);
        if ( !remain )
            FatalErr("BinaryReader attempted to read past the end of file "
                      << mFR.getFilename());

        if ( remain > len )
            remain = len;
        memcpy(buf, mpBuf, remain);
        buf += remain;
        len -= remain;
    }
    mpBuf += remain;
}

void BinaryReader::testToken()
{
    MagicToken tok;
    if ( !read(&tok).isValid() )
        FatalErr("Reading binary file " << mFR.getFilename()
                  << " failed: Initial token is invalid.");
}
/*
template<class Itr>
void BinaryWriter::writeItr(Itr begin, Itr const& end ) { while (begin != end ) { write(*begin); ++begin; } }

template<class T>
void BinaryWriter::write(T const* begin, T const* end ) { writeArray(begin, end, typename Serializability<T>::type()); }

template<class Itr>
void BinaryReader::readItr(Itr begin, Itr const& end ) { while (begin != end ) { read(&*begin); ++begin; } }

template<class T>
void BinaryReader::read(T* begin, T* end ) { readArray(begin, end, typename Serializability<T>::type()); }
 */