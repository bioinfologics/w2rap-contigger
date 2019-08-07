/*
 * Hash.h
 *
 *  Created on: Mar 19, 2014
 *      Author: tsharpe
 */

#ifndef MATH_HASH_H_
#define MATH_HASH_H_

#include <cstddef>

// FNV-1a algorithm (64-bit)
// Itr's value-type should be byte-like.
// The purpose of providing you an overridable initial value is to allow you
// to calculate a full FNV-1a piecewise:  Use the default initial value for your
// first chunk, and the previously returned hash value as the initial value for
// subsequent chunks.
template <class Itr>
uint64_t FNV1a( Itr itr, Itr const& end, uint64_t val = 14695981039346656037ul )
{
    while ( itr != end )
    {
        val = 1099511628211ul * (val ^ *itr);
        ++itr;
    }
    return val;
}


#endif /* MATH_HASH_H_ */
