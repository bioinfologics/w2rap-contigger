/*
 * \file Algorithms.h
 * \author tsharpe
 * \date Mar 18, 2010
 *
 * \brief
 */
#ifndef FEUDAL_ALGORITHMS_H_
#define FEUDAL_ALGORITHMS_H_

#include <cstddef>

/// Numbers of kmers in some double-vector
/// Itr has to point to something that has a size() method.
template <class Itr>
inline size_t kmerCount( Itr itr, Itr end, unsigned K )
{
    size_t result = 0;
    while ( itr != end )
    {
        if ( itr->size() >= K )
            result += itr->size() - K + 1;
        ++itr;
    }
    return result;
}

#endif /* FEUDAL_ALGORITHMS_H_ */
