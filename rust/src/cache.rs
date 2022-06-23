use crate::{M, N};
use crate::helper::{are_adjacent, gen_edges, i2e};
use crate::tedges::TEdge;

const CACHE_SIZE: usize = match N{4 => 9, 5 => 25, 6 => 75, 7 => 231, 8 => 763, _ => 75};

pub struct Cache {
    pub(crate) indicators: [u32; CACHE_SIZE],
    pub(crate) delimiters: [usize; M+1],
    pub(crate) edges: [TEdge; M],
    pub(crate) e2i: [[usize; N]; N],
    pub(crate) adjacent_bits: [u32; M],
}


#[test]
pub fn test() {
    println!("{:?}", compute_cache().indicators);
}

pub fn compute_cache() -> Cache {
    let mut nb = 0;
    let mut indicators: [u32; CACHE_SIZE] = [0 as u32; CACHE_SIZE];
    let mut delimiters: [usize; M+1] = [CACHE_SIZE; M+1];

    let mut delim = 0;
    for nb_new_edges in 0..(N/2)+1 {
        for ind in 1..2_u32.pow(M as u32) {
            if ind.count_ones() == nb_new_edges as u32 && is_independent(ind) {
                indicators[nb] = ind;
                nb += 1;
            }
        }
        delimiters[delim] = nb;
        delim += 1;
    }
    let mut e2i: [[usize; N]; N] = [[0; N]; N];
    let mut i = 0;
    for u in 0..N{
        for v in (u+1)..N{
            e2i[u][v] = i;
            e2i[v][u] = i;
            i += 1;
        }
    }
    let edges = gen_edges(N as u8);
    let mut adjacent_bits:[u32; M] = [0; M];
    for i in 0..M{
        let mut adj = 0;
        for j in 0..M{
            if are_adjacent(&edges[i], &edges[j]){
                adj |= 1 << j;
            }
        }
        adjacent_bits[i] = adj;
    }
    Cache{indicators, delimiters, edges, e2i, adjacent_bits}
}

pub fn get_indicators(nb_cand_edges: usize, cache: &Cache) -> &[u32] {
    return &cache.indicators[0..cache.delimiters[nb_cand_edges]];
}

pub fn is_independent(edge_ind: u32) -> bool{
    let mut nodes: [u8; N] = [0; N];
    for i in 0..M{
        if edge_ind & 1 << i != 0 {
            let e = i2e(N as u8, i as u8);
            if nodes[e.0 as usize] == 1 || nodes[e.1 as usize] == 1 {
                return false;
            } else {
                nodes[e.0 as usize] = 1;
                nodes[e.1 as usize] = 1;
            }
        }
    }
    true
}

#[test]
pub fn clementine(){
    for i in 0..10{
        println!("{}", i);
    }
}