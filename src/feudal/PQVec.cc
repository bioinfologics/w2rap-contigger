/*
 * PQVec.cc
 *
 *  Created on: Jan 30, 2016
 *      Author: tsharpe
 */

#include "feudal/PQVec.h"

void save_quals(VecPQVec const &vpqv, std::string filename){
    std::ofstream batch_file(filename,std::ios::out | std::ios::trunc | std::ios::binary);
    uint64_t vsize=vpqv.size();
    batch_file.write((const char *)&vsize, sizeof(vsize));
    for (auto &pq:vpqv){
        batch_file.write((const char *)&(pq.mSize), sizeof(pq.mSize));
        batch_file.write((const char *) pq.mData, pq.mSize);
    }
    batch_file.close();
}
void load_quals(VecPQVec &pqv, std::string filename){
    std::ifstream batch_file(filename,std::ios::in | std::ios::binary);
    uint64_t vsize;
    batch_file.read((char *)&vsize, sizeof(vsize));
    pqv.clear();
    pqv.reserve(vsize);
    while (vsize--){
        PQVec pq;
        batch_file.read(( char *)&(pq.mSize), sizeof(pq.mSize));
        pq.mData=(uint8_t *) malloc(pq.mSize);
        batch_file.read(( char *) pq.mData, pq.mSize);
        pqv.push_back(std::move(pq));
    }
    batch_file.close();

}
