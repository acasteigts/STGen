#ifndef _STGEN_TGRAPH
#define _STGEN_TGRAPH
#include <bitset>
#include "tedge.hpp"

/* This file contains the minimal version of the TGraph class.
You can replace the tgraph.hpp file with this one if you want
to try it. */    

namespace stgen{


///////////////////
////* TGraph */////

template<int N, int M=(N*(N-1)/2)>
struct TGraph {
    bitset<M> edges;
    bitset<M> max_edges;

    TGraph(){
        max_edges.set();
    }
    TGraph(const TGraph& g, const std::bitset<M> matching){
        edges = g.edges | matching;
        max_edges = matching;
    }

};


}
#endif