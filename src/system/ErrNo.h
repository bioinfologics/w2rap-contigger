/*
 * \file ErrNo.h
 * \author tsharpe
 * \date Dec 13, 2011
 *
 * \brief Handle the errno variable safely.
 */
#ifndef SYSTEM_ERRNO_H_
#define SYSTEM_ERRNO_H_

#include <ostream>
#include <string>

class ErrNo
{
public:
    /// initialized with current value of errno global
    ErrNo();

    /// initialized with whatever value you want
    explicit ErrNo( int err ) : mErrNo(err) {}

    // compiler-supplied copying and destructor are OK

    /// the value of the stored errno
    int val() const { return mErrNo; }

    /// thread-safe method to produce the text associated with the errno
    std::string text() const;

    friend std::ostream& operator<<( std::ostream& os, ErrNo const& errNo )
    { return os << errNo.text(); }

private:
    int mErrNo;
};

#endif /* SYSTEM_ERRNO_H_ */
