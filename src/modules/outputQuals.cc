#include <iostream>
#include <feudal/PQVec.h>

int main(int argc, char **argv) {
    VecPQVec quals;
    load_quals(quals, "./pe_data.cqual");

    for (uint64_t i = 0; i < quals.size(); i++) {
        QualVec pqv;
        quals[i].unpack(&pqv);
        std::cout << "Qual " << i+1 << std::endl;
        for (const auto q : pqv) {
            std::cout << q;
        }
        std::cout << std::endl;
    }

    return 0;
}