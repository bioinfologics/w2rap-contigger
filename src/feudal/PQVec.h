/*
 * PQVec.h, new implementation by bj
 *
 *  Created on: Jan 30, 2017
 *      Author: bjclavijo
 */

#ifndef PQVEC_H_
#define PQVEC_H_

#include "Qualvector.h"
#include <algorithm>
#include <cstddef>
#include <cstring>
#include <iterator>
#include <memory>
#include <numeric>

// helper class to do the buffer encoding and decoding// a compressed QualVec.
// the only way you can make one is with a QualVec (or by copying).
// the only way you can use one is to turn it into a QualVec.

class PQVec
{
public:
    explicit PQVec( QualVec const& qv ) {
        //Encode from qv
        mSize=0;
        if (0==qv.size()) return;
        uint8_t buff[qv.size()];
        uint8_t last_qual=qv[0];
        uint8_t lqcount=1;
        for (uint16_t i=1;i<qv.size();++i){
            if (qv[i]!=last_qual or lqcount==4){
                buff[mSize]=last_qual+64*(lqcount-1);
                ++mSize;
                last_qual=qv[i];
                lqcount=1;
            }
            else ++lqcount;
        }
        buff[mSize]=last_qual+64*(lqcount-1);
        ++mSize;
        mData=(uint8_t *) malloc(mSize);
        memcpy(mData,buff,mSize);
    }

    ~PQVec(){
        if (mSize) free(mData);
    }

    void unpack( QualVec &qv ) const
    {
        //Decode into qv
        //first pass: get the size
        uint16_t esize=0;
        for (auto i=0;i<mSize;i++) esize+=mData[i]/64+1;
        qv.resize(esize);
        uint16_t epos=0;
        for (auto i=0;i<mSize;i++){
            uint8_t q=mData[i]%64;
            for (auto p=mData[i]/64+1;p>0;--p){
                qv[epos]=q;
                ++epos;
            }
        }
    }

    inline void unpack( QualVec * qv ) const {
        unpack(*qv);
    }

    //operator QualVec() const { QualVec qv; unpack(&qv); return qv; }

    // all the rest of this crap is boilerplate


private:
    uint8_t * mData;
    uint16_t mSize;
};


//using PQVec = PQVecA<>;
using VecPQVec = std::vector<PQVec>;



#endif /* PQVEC_H_ */
