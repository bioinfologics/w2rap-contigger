//
// Created by Gonzalo Garcia (TGAC) on 01/11/2016.
//

#ifndef W2RAP_CONTIGGER_TESTCODE_2_H
#define W2RAP_CONTIGGER_TESTCODE_2_H

#include <iostream>
#include <gcc5/c++/cstdint>
#include <vector>
#include <map>
#include "paths/HyperBasevector.h"

class EdgeDictionary{
public:
    EdgeDictionary(HyperBasevector &hbv);
    ~EdgeDictionary();

private:
    HyperBasevector* pHbv;
    std::map<std::uint64_t, std::vector<std::uint32_t>> kmap;

};


#endif //W2RAP_CONTIGGER_TESTCODE_2_H
