//
// Created by Luis Yanes (EI) on 2019-02-28.
//

#ifndef W2RAP_CONTIGGER_DISJOINTSET_HPP
#define W2RAP_CONTIGGER_DISJOINTSET_HPP

#include <vector>
#include <Vec.h>

/**
 * This class aims to solve the performance problems of the equiv_rel_template
 * While the idea is good and it works for small classes and a small totoal number of elements, the joins are terribly expensive
 * when classes become large.
 * This solves that by creating a heap structure. And the price is only payed once, when flatten and link is called
 *
 * @tparam INT
 */
template<class INT>
class equiv_rel_template_bj {
public:
    equiv_rel_template_bj(INT n){
        up.resize(n);
        next_in_class.resize(n);
        for (INT x=0;x<n;++x) {
            up[x]=x;
            next_in_class[x]=x;
        }
    }

    INT class_of(INT x){
        INT c;
        for (c=x;c!=up[c];c=up[c]);
        return c;
    }
    void join(INT x, INT y){
        if (x==y) return;
        //since we're going to compute the class of x and y we may as well flatten them
        up[x]=class_of(x);
        up[y]=class_of(y);
        if (up[x]<up[y]) up[up[y]]=up[x];
        else up[up[x]]=up[y];
    }
    void flatten_and_link(){
        for(INT x=0;x<up.size();++x) next_in_class[x]=x;

        for(INT x=up.size()-1;x>=0;--x){
            up[x]=class_of(x);
            std::swap(next_in_class[up[x]],next_in_class[x]);
        }
    }

    void OrbitRepsAlt( vec<INT>& reps ) const {
        reps.clear( );
        for ( INT i = 0; i < (INT) up.size( ); i++ ) if ( i == up[i] ) reps.push_back(i);
    }

    void Orbit( INT a, vec<INT>& o ) const {
        o.resize(0);
        o.push_back(a);
        INT b = a;
        while(1) {
            b = next_in_class[b];
            if ( b == a ) break;
            o.push_back(b);
        }
    }

private:
    std::vector<INT> up, next_in_class;
};

class DisjointSet {

public:
    std::vector<int> parent;
    std::vector<int> rank;
    std::vector<int> size;

    explicit DisjointSet(int N) : parent(N), rank(N,0), size(N, 1) {
        for (int n = 0; n < N; n++) {
            parent[n] = n;
        }
    }

    int find(int x);

    void link_sets(int i, int j);

    // Link - union the two sets represented by x and y
    void link(int x, int y);

    // Union-Set - union the two sets containing x and y
    void union_set(int x, int y);

    std::size_t count_sets();

    void normalise_sets();

    void normalize_node(int i);


    void OrbitRepsAlt(std::vector<int> &reps);

    void Orbit(int a, std::vector<int> &orbit);
};

#endif //W2RAP_CONTIGGER_DISJOINTSET_HPP
