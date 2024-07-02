#ifndef _STGEN_TGRAPH
#define _STGEN_TGRAPH
#include <bitset>
#include <array>
#include "tedge.hpp"

namespace stgen{


///////////////////
////* TGraph */////

template<int N, int M=(N*(N-1)/2)>
struct TGraph {
    std::bitset<M> edges;
    std::bitset<M> max_edges;
    TEdge<N> tedges[M];
    int lifetime;
    int nb_edges;

    TGraph(){
        max_edges.set();
        lifetime = 0;
        nb_edges = 0;
    }
    TGraph(const TGraph& g, const std::bitset<M> matching){
        edges = g.edges | matching;
        max_edges = matching;
        nb_edges = g.nb_edges;
        std::copy(g.tedges, g.tedges+nb_edges, tedges);
        lifetime = g.lifetime + 1;
        for (int i=matching._Find_first(); i<M; i=matching._Find_next(i)){
            TEdge<N> e = i2te<N>(i, lifetime);
            tedges[nb_edges++] = e;
        }         
    }

    bool is_tc() const{
        bitset<N> predecessors[N];
        for (int i=0; i<N; i++) {
            predecessors[i].set(i);
        }
        for (int i=0; i<nb_edges; i++){
            TEdge<N> e = tedges[i];
            predecessors[e.u] |= predecessors[e.v];
            predecessors[e.v] = predecessors[e.u];
        }
        for (int i=0; i<N; i++){
            if (! predecessors[i].all()){
                return false;
            }
        }
        return true;
    }

};


}
#endif