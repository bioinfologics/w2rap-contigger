/*
 * \file IsSizeT.h
 * \author tsharpe
 * \date Oct 19, 2009
 *
 * \brief
 */
#ifndef ISSIZET_H_
#define ISSIZET_H_

struct YesSizeT {};
struct NoSizeT {};

template <class T>
struct IsSizeT : public NoSizeT {};

template<>
struct IsSizeT<int> : public YesSizeT {};

template<>
struct IsSizeT<unsigned int> : public YesSizeT {};

template<>
struct IsSizeT<long> : public YesSizeT {};

template<>
struct IsSizeT<unsigned long> : public YesSizeT {};

#endif /* ISSIZET_H_ */
