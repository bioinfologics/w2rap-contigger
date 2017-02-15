///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file SpinLockedData.h
 * \author tsharpe
 * \date Nov 16, 2011
 * \new spinlock by bj on Feb 11, 2016
 * \brief
 */
#ifndef SYSTEM_SPINLOCKEDDATA_H_
#define SYSTEM_SPINLOCKEDDATA_H_

#include "system/Assert.h"
#include <atomic>
#include <iostream>


/// A spin-lock.
class SpinLockedData
{
    std::atomic<bool> m_flag;

public:
    SpinLockedData() : m_flag(false) {};
    inline void lock()     noexcept {   while(std::atomic_exchange_explicit(&m_flag, true, std::memory_order_acquire)); }
    inline void unlock()   noexcept {   std::atomic_exchange_explicit(&m_flag, false, std::memory_order_release); }
};

/// Something that operates a spin-lock, and never forgets to unlock it.
class SpinLocker
{
public:
    SpinLocker( SpinLockedData& lock ) : mLock(lock) { mLock.lock(); }
    ~SpinLocker() { mLock.unlock(); }

private:
    SpinLocker( SpinLocker const& ); // unimplemented -- no copying
    SpinLocker& operator=( SpinLocker const& ); // unimplemented -- no copying

    SpinLockedData& mLock;
};

#endif /* SYSTEM_SPINLOCKEDDATA_H_ */
