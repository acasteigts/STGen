#ifndef _STGEN_ITER
#define _STGEN_ITER
#include <bitset>
#include <vector>
#include <array>
#include <numeric>
#include <execution>
#include <algorithm>
#include <math.h>
#include <iostream>
#include <random>
#include "traverse.hpp"

namespace stgen{
using namespace std;



/////////////////////////
// Sequential iterator 

template<int N>
class TGraphIterator{
    bool as_symmetric;
public:
    TGraphIterator(bool (*process)(const TGraph<N>&, int64_t&)){
        process_as_tgraph<N> = process;
        as_symmetric = false;
    }

    TGraphIterator(bool (*process)(const SymmetricTGraph<N>&, int64_t&)){
        process_as_symmetric<N> = process;
        as_symmetric = true;
    }

    int64_t execute(){
        auto g = SymmetricTGraph<N>();
        return visit_graph(g, as_symmetric);
    }
};



///////////////////////
// Parallel iterator

template<int N>
class TGraphParIterator{
    static vector<TGraph<N>> boundary;

    static bool process_sym(const SymmetricTGraph<N>& g, int64_t& nb){
        if (g.auts.empty()){
            boundary.push_back(g.tg);
            return false;
        } else {
            return usr_process_as_tgraph<N>(g.tg, nb);
        }
    }

public:
    TGraphParIterator(bool (*process)(const TGraph<N>&, int64_t&)){
        usr_process_as_tgraph<N> = process;
    }

    int64_t execute(){
        boundary.clear();
        process_as_symmetric<N> = process_sym;
        auto g = SymmetricTGraph<N>();
        int64_t nb = visit_graph<N>(g, true);
        process_as_tgraph<N> = usr_process_as_tgraph<N>;
        nb += transform_reduce(execution::par_unseq, boundary.begin(), boundary.end(), 0LL, plus{}, [](TGraph<N> g) {
            return visit_graph<N>(g);
        });
        return nb;
    }
};
template<int N>
vector<TGraph<N>> TGraphParIterator<N>::boundary;






////////////////////////////
// Parallel batch iterator

template<int N>
class TGraphBatchIterator{
    static const int THRESHOLD_1 = N;
    static const int THRESHOLD_2 = N+4;
    static vector<SymmetricTGraph<N>> boundary_1;
    static vector<TGraph<N>> boundary_2;
    int nb_graphs_in_top = 0; // user-defined count in the top of the tree
    int nb_batches = 0; // number of batches for subsequent execution

private:
    static bool process_top(const SymmetricTGraph<N>& g, int64_t& nb){
        if (g.tg.edges.count() < THRESHOLD_1){
            return usr_process_as_tgraph<N>(g.tg, nb);
        } else {
            boundary_1.push_back(g);
            return false;
        }
    };
    static bool process_batch(const SymmetricTGraph<N>& g, int64_t& nb){
        if ((! g.auts.empty()) || (g.tg.edges.count() < THRESHOLD_2)){
            return usr_process_as_tgraph<N>(g.tg, nb);
        } else {
            boundary_2.push_back(g.tg);
            return false;
        }
    };

public:
    TGraphBatchIterator(bool (*process)(const TGraph<N>&, int64_t&)){
        usr_process_as_tgraph<N> = process;
        process_as_symmetric<N> = process_top;
        boundary_1.clear();
        auto g = SymmetricTGraph<N>();
        nb_graphs_in_top = visit_graph(g, true);
        nb_batches = boundary_1.size(); 
    }

    int number_of_batches() const{
        return nb_batches;
    }
    int64_t execute(int num_batch){
        boundary_2.clear();
        process_as_symmetric<N> = process_batch;
        auto h = boundary_1[num_batch];
        int64_t nb = visit_graph(h, true);
        process_as_tgraph<N> = usr_process_as_tgraph<N>;
        nb += transform_reduce(execution::par_unseq, boundary_2.begin(), boundary_2.end(), 0LL, plus{}, [](TGraph<N> g) {
            return visit_graph(g);
        });
        if (num_batch==0){
            nb += nb_graphs_in_top;
        }
        return nb;
    }
    int64_t execute_all(){
        cout << nb_batches << " batches to be executed..." << endl;
        int64_t total = 0;
        process_as_symmetric<N> = process_batch;
        for (int num_batch = 0; num_batch < nb_batches; num_batch++){
            int64_t nb = execute(num_batch);
            cout << "batch " << num_batch << "/" << nb_batches-1 << ": " << nb << endl;
            total += nb;
        }
        total += nb_graphs_in_top;
        cout << "total: " << total << " (including " << nb_graphs_in_top << " in the top before the batches)" << endl;
        return total;
    }
};
template<int N>
vector<SymmetricTGraph<N>> TGraphBatchIterator<N>::boundary_1;
template<int N>
vector<TGraph<N>> TGraphBatchIterator<N>::boundary_2;

}
#endif