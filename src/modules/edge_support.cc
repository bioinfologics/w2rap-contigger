#include <iostream>
#include <Vec.h>
#include <paths/HyperBasevector.h>
#include <paths/long/ReadPath.h>
#include <util/OutputLog.h>

vec<int> get_supports(const HyperBasevector &hb, vec<int> &inv, ReadPathVec &paths) {
    vec<int> support(hb.EdgeObjectCount(), 0);
    for (const auto & path : paths) {
        for (int64_t j = 0; j < (int64_t) path.size(); j++) {
            int e = path[j];
            if (j >= 1) support[e]++;
            if (inv[e] >= 0 && j < (int64_t) path.size() - 1)
                support[inv[e]]++;
        }
    }
    return support;
}

int main(int argc, char **argv) {
    std::cout << "Hello world" << std::endl;
    std::vector<int> edges;
    HyperBasevector hbv;
    vec<int> hbvinv;
    ReadPathVec paths;

    BinaryReader::readFile(argv[1], &hbv);
    //Create inversion
    OutputLog(4) <<"Creating graph involution..." << std::endl;
    hbvinv.clear();
    hbv.Involution(hbvinv);
    //load paths
    OutputLog(2) <<"Loading paths..." << std::endl;
    LoadReadPathVec(paths,argv[2]);

    auto supports=get_supports(hbv,hbvinv,paths);
    for (int i=3; i < argc; i++) {
        std::cout << "edge" << argv[i] << ": " << supports[std::stoi(argv[i])] << std::endl;
    }

    return 0;
}