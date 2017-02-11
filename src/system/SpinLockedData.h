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
    std::atomic_flag m_flag = ATOMIC_FLAG_INIT;

public:
    void lock()     noexcept {   while(m_flag.test_and_set(std::memory_order_acquire)); }
    void unlock()   noexcept {         m_flag.clear(std::memory_order_release);         }
    bool try_lock() noexcept { return !m_flag.test_and_set(std::memory_order_acquire);  }
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
