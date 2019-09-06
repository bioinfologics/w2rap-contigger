#include <iostream>
#include <Vec.h>
#include <paths/HyperBasevector.h>
#include <paths/long/ReadPath.h>
#include <util/OutputLog.h>
#include <unordered_map>

vec<int> get_supports(const HyperBasevector &hb, vec<int> &inv, ReadPathVec &paths, std::unordered_map<int,std::vector<int64_t>> &edge_support, std::unordered_map<int,std::vector<int64_t>> &edge_paths) {
    vec<int> support(hb.EdgeObjectCount(), 0);
    int64_t p(0);
#pragma omp parallel for
    for (int64_t p = 0; p < paths.size(); p++) {
        for (int64_t j = 0; j < (int64_t) paths[p].size(); j++) {
            int e = paths[p][j];
            auto edgep_lookup(edge_paths.find(e));
            auto edgeInvp_lookup(edge_paths.find(inv[e]));
            if (edgep_lookup != edge_paths.cend()) {
                edgep_lookup->second.push_back(p);
            }
            if (edgeInvp_lookup != edge_paths.cend()) {
                edgeInvp_lookup->second.push_back(p);
            }
            if (j >= 1) {
                support[e]++;
                auto edge_lookup(edge_support.find(e));
                if (edge_lookup != edge_support.cend()) {
#pragma omp critical
                    edge_lookup->second.push_back(p);
                }
            }
            if (inv[e] >= 0 && j < (int64_t) paths[p].size() - 1) {
                auto edge_lookup(edge_support.find(inv[e]));
                if (edge_lookup != edge_support.cend()) {
#pragma omp critical
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

    std::unordered_map<int, std::vector<int64_t>> edge_support;
    std::unordered_map<int, std::vector<int64_t>> edges_paths;
    for (int i=4; i < argc; i++) {
        edge_support[std::stoi(argv[i])];
        edges_paths[std::stoi(argv[i])];
    }
    auto supports=get_supports(hbv,hbvinv,paths,edge_support,edges_paths);

    for (const auto &edge_paths : edge_support) {
        std::cout << "edge"<<edge_paths.first<<": " << supports[edge_paths.first] << " paths support this edge"<<std::endl;
        for (const auto &path : edge_paths.second) {
            std::cout << path <<" " << bases[path].ToString() << "\n";
        }
        std::cout << std::endl;
    }
    for (const auto &edge_paths : edges_paths) {
        std::cout << "edge"<<edge_paths.first<<": " << supports[edge_paths.first] << " paths contain this edge"<<std::endl;
        for (const auto &path : edge_paths.second) {
            std::cout << path <<" " << bases[path].ToString() << "\n";
        }
        std::cout << std::endl;
    }

    return 0;
}