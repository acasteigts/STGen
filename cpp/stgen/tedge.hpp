#ifndef _STGEN_TEDGE
#define _STGEN_TEDGE
#include <bitset>
#include <vector>
#include <array>
#include <numeric>
#include <algorithm>
#include <math.h>
#include <iostream>
#include <random>


namespace stgen{
using namespace std;

///////////////
///* Edge *////

template<int N>
struct Edge {
    int8_t u;
    int8_t v;
};
template<int N>
struct TEdge {
    int8_t u;
    int8_t v;
    int8_t t;
};

template<int N>
static constexpr int e2i(int u, int v){
    int min = std::min(u,v);
    int max = std::max(u,v);
    return ((min*N - (min*(min+1)) / 2) + max-min-1);
}

template<int N>
static constexpr Edge<N> i2e(int i){
    int k = 1;
    while (i >= (N - k)){
        i -= N - k;
        k += 1;
    }
    return Edge<N>({(int8_t)(k - 1), (int8_t)(k + i)});
}

template<int N>
static constexpr TEdge<N> i2te(int i, int t=0){
    int k = 1;
    while (i >= (N - k)){
        i -= N - k;
        k += 1;
    }
    return TEdge<N>({(int8_t)(k - 1), (int8_t)(k + i), (int8_t) t});
}

template<int N, int M=N*(N-1)/2>
const array<bitset<M>,M> ADJACENT_EDGES = ([](){
    array<bitset<M>,M> adjacent_edges;
    if (N==2){
        adjacent_edges[0].set(0);
    } else {
        for (int i=0; i<M; i++){
            bitset<M> bits;
            Edge e = i2e<N>(i);
            for (int j=0; j<M; j++){
                if (j!=i){
                    Edge f = i2e<N>(j);
                    if ((e.u == f.u) || (e.u == f.v) || (e.v == f.u) || (e.v == f.v)){
                        bits.set(j);
                    }
                }
            }
            adjacent_edges[i]=bits;
        }
    }
    return adjacent_edges;
})();


template<int N, int M=N*(N-1)/2>
bool is_matching(const bitset<N*(N-1)/2> edge_ind){
    int nodes[N]={0};
    for (int i=edge_ind._Find_first(); i<M; i=edge_ind._Find_next(i)){
        if (edge_ind.test(i)) {
            auto e = i2e<N>(i);
            if (nodes[e.u] == 1 || nodes[e.v] == 1) {
                return false;
            } else {
                nodes[e.u] = 1;
                nodes[e.v] = 1;
            }
        }
    }
    return true;
}

template<int N, int M=N*(N-1)/2>
void build_matchings(vector<bitset<M>>& matchings, int numbits, uint64_t acc){
    if (numbits == 0) {
        if (is_matching<N>(acc)){
            matchings.push_back(acc);
        }
        return;
    }
    for (int i=0; i<M; i++){
        if (acc < (1ull << i)) {
            build_matchings<N,M>(matchings, numbits - 1, acc | (1ull << i));
        }
    }
}

template<int N, int M=N*(N-1)/2>
const vector<bitset<M>> MATCHINGS = ([](){
    vector<bitset<M>> matchings;
    for (int numbits = 1; numbits <= N/2; numbits++) {
        uint64_t acc = 0;
        build_matchings<N,M>(matchings, numbits, acc);
    }
    vector<bitset<M>> filtered; // Delete useless matchings
    for (int numbits = 1; numbits <= N/2; numbits++){
        bool first_edge_hit = false;
        for (auto m : matchings){
            if (m.count() == numbits){
                if (m.test(0)){
                    if (! first_edge_hit){
                        bitset<M> bs(m);
                        filtered.push_back(bs);
                        first_edge_hit = true;
                    }
                } else {
                    bitset<M> bs(m);
                    filtered.push_back(bs);
                }
            }
        }
    }
    return filtered;
})();

template<int N, int M=N*(N-1)/2>
const array<int,M+1> DELIM = ([](){
    array<int,M+1> delimiters;
    int insert = 0;
    int numbits = 1;
    delimiters[insert++] = 0;
    for (int i=0; i<MATCHINGS<N>.size(); i++){
        bitset<M> bs = MATCHINGS<N>[i];
        if (bs.count() > numbits){
            delimiters[insert++] = i;
            numbits = bs.count();
        }
    }
    while (insert < M+1){
        delimiters[insert++] = MATCHINGS<N>.size();
    }
    return delimiters;
})();

template<int N, int M=N*(N-1)/2>
static bool is_stable_under_automorphism(const bitset<M> matching, const array<int8_t, M> aut){
    for (int i=matching._Find_first(); i<M; i=matching._Find_next(i)){
        if (! matching.test(aut[i])){
            return false;
        }
    }
    return true;
}

template<int N, int M=N*(N-1)/2>
static bool are_equivalent_matchings(const bitset<M> m1, const bitset<M> m2, const array<int8_t,M> aut){
    for (int i=m1._Find_first(); i<M; i=m1._Find_next(i)){
        if (! m2.test(aut[i])){
            return false;
        }
    }
    return true;
}

}
#endif
