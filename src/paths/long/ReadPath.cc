//
// Created by Bernardo Clavijo (TGAC) on 22/10/2016.
//
#include "ReadPath.h"

void WriteReadPathVec(const ReadPathVec &rpv, const char * filename){
    std::ofstream f(filename, std::ios::out | std::ios::trunc | std::ios::binary);
    uint64_t pathcount=rpv.size();
    f.write((const char *) &pathcount, sizeof(pathcount));
    uint16_t ps;
    int mOffset;
    for (auto const &rp:rpv){
        mOffset=rp.getOffset();
        ps=rp.size();
        f.write((const char *) &mOffset, sizeof(mOffset));
        f.write((const char *) &ps, sizeof(ps));
        f.write((const char *) rp.data(),ps*sizeof(int));
    }
    f.close();
}


void LoadReadPathVec(ReadPathVec &rpv, const char * filename){
    std::ifstream f(filename, std::ios::in | std::ios::binary);
    uint64_t pathcount;
    f.read((const char *) &pathcount, sizeof(pathcount));
    rpv.resize(pathcount);
    uint16_t ps;
    int mOffset;
    for (auto const &rp:rpv){
        f.read((const char *) &mOffset, sizeof(mOffset));
        rp.setOffset(mOffset);
        f.read((const char *) &ps, sizeof(ps));
        rp.resize(ps);
        f.read((const char *) rp.data(),ps*sizeof(int));
    }
    f.close();
}

