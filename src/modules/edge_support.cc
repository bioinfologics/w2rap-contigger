#include <iostream>
#include <Vec.h>
#include <paths/HyperBasevector.h>
#include <paths/long/ReadPath.h>
#include <util/OutputLog.h>
#include <unordered_map>

vec<int> get_supports(const HyperBasevector &hb, vec<int> &inv, ReadPathVec &paths, std::unordered_map<int,std::vector<int64_t>> &edges) {
    vec<int> support(hb.EdgeObjectCount(), 0);
    int64_t p(0);
    for (const auto & path : paths) {
        for (int64_t j = 0; j < (int64_t) path.size(); j++) {
            int e = path[j];
            if (j >= 1) {
                support[e]++;
                auto edge_lookup(edges.find(e));
                if (edge_lookup != edges.cend()) {
                    edge_lookup->second.push_back(p);
                }
            }
            if (inv[e] >= 0 && j < (int64_t) path.size() - 1) {
                auto edge_lookup(edges.find(inv[e]));
                if (edge_lookup != edges.cend()) {
                    edge_lookup->second.push_back(p);
                }
                support[inv[e]]++;
            }
        }
        p++;
    }
    return support;
}

int main(int argc, char **argv) {
    std::cout << "Hello world" << std::endl;

    // Load graph
    HyperBasevector hbv;
    BinaryReader::readFile(argv[1], &hbv);

    //Create inversion
    OutputLog(4) <<"Creating graph involution..." << std::endl;
    vec<int> hbvinv;
    hbv.Involution(hbvinv);

    //load paths
    OutputLog(2) <<"Loading paths..." << std::endl;
    ReadPathVec paths;
    LoadReadPathVec(paths,argv[2]);

    // Load reads
    vecbvec bases;
    bases.ReadAll(argv[3]);

    std::unordered_map<int, std::vector<int64_t>> edges;
    for (int i=3; i < argc; i++) {
        edges[std::stoi(argv[i])];
    }
    auto supports=get_supports(hbv,hbvinv,paths,edges);

    for (const auto &edge_paths : edges) {
        std::cout << "edge"<<edge_paths.first<<": " << supports[edge_paths.first] << " paths support this edge"<<std::endl;
        for (const auto &path : edge_paths.second) {
            std::cout << path << bases[path].ToString() << "\n";
        }
        std::cout << std::endl;
    }

    return 0;
}