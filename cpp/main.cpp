#include "stgen/stgen.h"

/* This program is a basic example of program using STGen 
to count all happy graphs on 6 vertices */

const int N = 6;

using namespace stgen;

bool count_all(const TGraph<N>& g, int64_t& nb){
    nb++;
    return true;
};


int main(){
    TGraphIterator<N> iter(count_all);
    int64_t nb = iter.execute();
    std::cout << nb << std::endl;
}