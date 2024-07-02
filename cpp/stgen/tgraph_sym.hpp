#ifndef _STGEN_TGRAPH_SYM
#define _STGEN_TGRAPH_SYM
#include <bitset>
#include "tgraph.hpp"


namespace stgen{
using namespace std;


/////////////////////
/* SymmetricTGraph */


template<int N, int M=(N*(N-1))/2>
struct SymmetricTGraph {
    TGraph<N> tg;
    vector<array<int8_t, M>> auts;

    SymmetricTGraph(){
        tg = TGraph<N>();
        auts.reserve(tgamma(N+1)); // N!
        array<int, N> p;
        iota(begin(p), end(p), 0);
        next_permutation(p.begin(), p.end());
        do {
            array<int8_t, M> aut;
            for (int i=0; i<M; i++){
                Edge e = i2e<N>(i);
                aut[i] = e2i<N>(p[e.u],p[e.v]);
            }
            auts.push_back(aut);
        } while (next_permutation(p.begin(), p.end()));
    }

    SymmetricTGraph(const SymmetricTGraph& g, const bitset<M> new_edges){
        tg = TGraph<N>{g.tg, new_edges};
        auts.reserve(auts.size());
        for (auto aut : g.auts){
            if (is_stable_under_automorphism<N,M>(new_edges, aut)){
                auts.push_back(aut);   
            }
        }
    }

    bool is_rigid() const{
        return auts.size()==0;
    }
    
public:
    bool contains_equivalent_matching(const vector<bitset<M>>& matchings, const bitset<M> m) const{
        for (auto m2 : matchings) { 
            if (m.count() == m2.count()) {
                for (auto aut : auts) {
                    if (are_equivalent_matchings<N,M>(m, m2, aut)) {
                        return true;
                    }
                }
            }
        }
        return false;
    }

    vector<bitset<M>> successors_as_matchings(const bitset<M> cand_edges, int matching_size) const{
        vector<bitset<M>> matchings;
        matchings.reserve(N); // Empiric
        for(auto it = MATCHINGS<N>.begin()+DELIM<N>[matching_size-1]; it != MATCHINGS<N>.begin()+DELIM<N>[matching_size]; it++){
            bitset<M> m = *it;
            if ((cand_edges | m) == cand_edges){
                if (!contains_equivalent_matching(matchings, m)) {
                    matchings.push_back(m);
                }
            }
        }
        return matchings;
    }
};


}
#endif