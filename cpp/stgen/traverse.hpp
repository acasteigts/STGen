#ifndef _STGEN_TRAVERSE
#define _STGEN_TRAVERSE
#include <bitset>
#include <vector>
#include <array>
#include <numeric>
#include <execution>
#include <algorithm>
#include <math.h>
#include <iostream>
#include <random>
#include "tgraph_sym.hpp"

namespace stgen{
using namespace std;


//////////////////////////////////////////////////////////
//* User-defined functions for processing each graph *////

template<int N>
bool default_process_as_tgraph(const TGraph<N>& g, int64_t& nb){
    nb++;
    return true;
};
template<int N>
auto process_as_tgraph = default_process_as_tgraph<N>;
template<int N>
auto usr_process_as_tgraph = default_process_as_tgraph<N>;
template<int N>

bool default_process_as_symmetric(const SymmetricTGraph<N>& g, int64_t& nb){
    nb++;
    return true;
};
template<int N>
auto process_as_symmetric = default_process_as_symmetric<N>;
template<int N>
auto usr_process_as_symmetric = default_process_as_symmetric<N>;



//////////////////////////////////////////////
//* Depth-first generation of the graphs *////


template<int N, int M=(N*(N-1)/2)>
bitset<M> candidate_edges(const TGraph<N>& g, int& max_new){
    bitset<M> cand_edges;
    for (int i=g.max_edges._Find_first(); i<M; i=g.max_edges._Find_next(i)){
        cand_edges |= ADJACENT_EDGES<N>[i];
    }
    cand_edges &= ~g.edges;
    if (cand_edges.any()){
        int max_edges_count = g.max_edges.count();
        max_new = max_edges_count * 2;
        if (max_new > cand_edges.count()){
            max_new = cand_edges.count();
        }
        if (max_new > N/2){
            max_new = N/2;
        }
    }
    return cand_edges;
}


template<int N, int M>
int64_t visit_successors(const TGraph<N>& g, const bitset<N*(N-1)/2>, int);

template<int N, int M>
int64_t visit_successors(const SymmetricTGraph<N>& g, const bitset<N*(N-1)/2>, int, bool);



// Visit (TGraph)
template<int N, int M=(N*(N-1)/2)>
int64_t visit_graph(const TGraph<N>& g){
    int64_t nb = 0;
    if (process_as_tgraph<N>(g, nb)){
        if (!g.edges.all()){
            int max_new;
            bitset<M> cand_edges = candidate_edges(g, max_new);
            if (max_new > 0){
                if (max_new ==1){
                    const TGraph h = TGraph<N>{g, cand_edges};
                    nb += visit_graph<N>(h);
                } else {
                    nb += visit_successors(g, cand_edges, max_new);
                }
            }
        }
    }
    return nb;
}

// Visit (SymmetricTGraph)
template<int N, int M=(N*(N-1)/2)>
int64_t visit_graph(const SymmetricTGraph<N>& g, bool as_symmetric = false){
    int64_t nb = 0;
    if (! as_symmetric && g.auts.empty()){
        nb += visit_graph<N>(g.tg);
    } else {
        bool keep;
        if (as_symmetric){
            keep = process_as_symmetric<N>(g, nb);
        } else {
            keep = process_as_tgraph<N>(g.tg, nb);
        }
        if (keep){
            if (!g.tg.edges.all()){
                int max_new;
                bitset<M> cand_edges = candidate_edges(g.tg, max_new);
                if (max_new > 0){
                    nb += visit_successors(g, cand_edges, max_new, as_symmetric);
                }
            }
        }
    }
    return nb;
}




// Visit Successors (TGraph)
template<int N, int M=(N*(N-1))/2>
int64_t visit_successors_rec(const TGraph<N>& g, const bitset<N*(N-1)/2>& cand_edges, int nb_new, const uint64_t acc){
    int64_t nb = 0;
    if (nb_new == 0) {
        if (is_matching<N>(acc)){
            nb += visit_graph<N>(TGraph<N>{g, acc});
        }
    } else {
        for (int i=cand_edges._Find_first(); i<M; i=cand_edges._Find_next(i)){
            if (acc < (1ull << i)) {
                nb += visit_successors_rec<N>(g, cand_edges, nb_new - 1, acc | (1ull << i));
            }
        }
    }
    return nb;
}
template<int N, int M=(N*(N-1))/2>
int64_t visit_successors(const TGraph<N>& g, const bitset<N*(N-1)/2> cand_edges, int max_new){
    int64_t nb = 0;
    if (cand_edges.count() < N/2){ // Empiric
        for (int nb_new = 1; nb_new <= max_new; nb_new++) {
            nb += visit_successors_rec<N>(g, cand_edges, nb_new, 0L);
        }
    } else {
        for (int msize = 1; msize <= max_new; msize++){
            // bool found = false;
            for(auto it = MATCHINGS<N>.begin()+DELIM<N>[msize-1]; it != MATCHINGS<N>.begin()+DELIM<N>[msize]; it++){
                bitset<M> matching = *it;
                if ((cand_edges | matching) == cand_edges){
                    nb += visit_graph(TGraph<N>{g, matching});
                    // found = true;
                }
            }
            // if (!found){break;}
        }
    }
    return nb;
}


// Visit Successors (SymmetricTGraph)
template<int N, int M=(N*(N-1))/2>
int64_t visit_successors(const SymmetricTGraph<N>& g, const bitset<N*(N-1)/2> cand_edges, int max_new, bool as_symmetric){
    int64_t nb = 0;
    vector<bitset<M>> matchings;
    matchings.reserve(2*M); // Empiric estimation
    for (int msize = 1; msize <= max_new; msize++){
        matchings.clear();
        bool found = false;
        for(auto it = MATCHINGS<N>.begin()+DELIM<N>[msize-1]; it != MATCHINGS<N>.begin()+DELIM<N>[msize]; it++){
            bitset<M> matching = *it;
            if ((cand_edges | matching) == cand_edges){
                if (!g.contains_equivalent_matching(matchings, matching)) {
                    matchings.push_back(matching);
                    found = true;
                }
            }
        }
        if (found){
            for (auto m : matchings){
                nb += visit_graph(SymmetricTGraph<N>{g, m}, as_symmetric);
            }
        } else {break;}
    }
    return nb;
}


}
#endif