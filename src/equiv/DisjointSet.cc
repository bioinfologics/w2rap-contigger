//
// Created by Luis Yanes (EI) on 2019-03-01.
//
#include "DisjointSet.hpp"


void DisjointSet::normalise_sets() {
    for (int i = 0; i < parent.size(); i++){
        normalize_node(i);
    }
}

void DisjointSet::OrbitRepsAlt(std::vector<int> &reps) {
    reps.clear();
    for (int i = 0; i < parent.size(); i++){
        if (parent[i] == i) reps.push_back(i);
    }
}

void DisjointSet::Orbit(int a, std::vector<int> &orbit) {
    orbit.resize(0);
    for (int i = 0; i < parent.size(); i++) {
        if (parent[i]==a) {
            orbit.push_back(i);
        }
    }
}

std::size_t DisjointSet::count_sets() {
    std::size_t count = 0;
    for (int i = 0; i < parent.size(); i++)
        if (parent[i] == i)
            ++count;
    return count;
}

void DisjointSet::union_set(int x, int y) {link(find(x), find(y));}

void DisjointSet::link(int x, int y) {link_sets(x, y);}

void DisjointSet::link_sets(int i, int j) {
    i = find(i);
    j = find(j);
    if (i == j) return;
    if (rank[i] > rank[j])
        parent[j] = i;
    else {
        parent[i] = j;
        if (rank[i] == rank[j])
            ++rank[j];
    }
}

int DisjointSet::find(int x) {
    auto old = x;
    auto ancestor = parent[old];
    while (ancestor != x) {
        x = ancestor;
        ancestor = parent[x];
    }
    x = parent[old];
    while (ancestor != x) {
        parent[old] = ancestor;
        old = x;
        x = parent[old];
    }
    return ancestor;
}

void DisjointSet::normalize_node(int i) {
    if (i > parent[i] || parent[parent[i]] != parent[i])
        parent[i] = parent[parent[i]];
    else {
        parent[parent[i]] = i;
        parent[i] = i;
    }
}
