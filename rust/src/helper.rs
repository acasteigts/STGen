use crate::{M, N};
use crate::tedges::TEdge;

pub fn gen_edges(n: u8) -> [TEdge; M] {
    let mut edges: [TEdge; M] = [TEdge(0, 0, 0); M];
    let mut i = 0;
    for u in 0..n {
        for v in (u + 1)..n {
            edges[i].0 = u;
            edges[i].1 = v;
            i += 1;
        }
    }
    edges
}

pub fn e2i(n: u8, u: u8, v: u8) -> usize {
    // julia: (u - 1) * (n - (u / 2)) + (v - u) // TODO simplify rust using this
    if u < v {
        ((u * n - (u * (u + 1)) / 2) + v - u - 1) as usize
    }else{
        ((v * n - (v * (v + 1)) / 2) + u - v - 1) as usize
    }
}
pub fn i2e(n: u8, mut i: u8) -> TEdge {
    let mut k = 1;
    while i >= (n - k){
        i -= n - k;
        k += 1;
    }
    TEdge(k - 1, k + i, 0)
}

pub fn bits_to_indices(edge_bits: u32) -> Vec<usize> {
    let mut edge_inds: Vec<usize> = vec![];
    for i in 0..M {
        if edge_bits & 1 << i != 0 {
            edge_inds.push(i);
        }
    }
    edge_inds
}
pub fn indices_to_bits(edge_inds: &[usize]) -> u32 {
    let mut edge_bits: u32 = 0;
    for i in edge_inds {
        edge_bits |= 1 << *i;
    }
    edge_bits
}

pub fn powerset<T>(s: &[T]) -> Vec<Vec<T>> where T: Clone {
    (0..2usize.pow(s.len() as u32)).map(|i| {
        s.iter().enumerate().filter(|&(t, _)| (i >> t) % 2 == 1)
            .map(|(_, element)| element.clone())
            .collect()
    }).collect()
}

pub fn are_adjacent(e: &TEdge, f: &TEdge) -> bool {
    (e.0 == f.0) || (e.0 == f.1) || (e.1 == f.0) || (e.1 == f.1)
}

pub fn edges_to_bits(n: u8, edges: &[TEdge; M]) -> u32 {
    let mut res: u32 = 0;
    for e in edges.iter(){
        let i = e2i(n, e.0, e.1);
        res |= 1 << i;
    }
    return res;
}
pub fn bits_to_edges(n: u8, bits: u32) -> Vec<TEdge> { // LALA TEST THIS!
    let mut res: Vec<TEdge> = vec![];
    for i in (0..32_u8).into_iter(){
        if bits & (1 << i) != 0{
            res.push(i2e(n, i));
        }
    }
    return res;
}

pub const fn each_ones() -> [u8; N]{
    let mut preds = [0; N];
    let mut i = 0;
    while i < N{
        preds[i] |= 1 << i;
        i += 1;
    }
    preds
}
pub const fn all_ones() -> u8 {
    (2_u32.pow(N as u32) - 1) as u8
}
