#include <iostream>
#include <Vec.h>
#include <paths/HyperBasevector.h>
#include <paths/long/ReadPath.h>
#include <util/OutputLog.h>
#include <unordered_map>

vec<int> get_supports(const HyperBasevector &hb, vec<int> &inv, ReadPathVec &paths, std::unordered_map<int,std::vector<int>> &edges) {
    vec<int> support(hb.EdgeObjectCount(), 0);
    for (const auto & path : paths) {
        for (int64_t j = 0; j < (int64_t) path.size(); j++) {
            int e = path[j];
            auto edge_lookup(edges.find(e));
            if (edge_lookup != edges.cend()) {
                edge_lookup->second.push_back(j);
            }
            if (j >= 1) support[e]++;
            if (inv[e] >= 0 && j < (int64_t) path.size() - 1)
                support[inv[e]]++;
        }
    }
    return support;
}

int main(int argc, char **argv) {
    std::cout << "Hello world" << std::endl;
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
    std::unordered_map<int, std::vector<int>> edges;
    for (int i=3; i < argc; i++) {
        edges[std::stoi(argv[i])];
    }
    auto supports=get_supports(hbv,hbvinv,paths,edges);

    for (const auto &edge_paths : edges) {
        std::cout << "edge"<<edge_paths.first<<":"<<std::endl;
        for (const auto &path : edge_paths.second) {
            std::cout << path << " ";
        }
        std::cout << std::endl;
    }

    return 0;
}