#![allow(dead_code)]
mod tedges;
mod helper;
mod cache;
mod tgraph;
mod algos;

const N: usize = 6;
const M: usize = N * (N - 1) / 2;
const DNMT: u8 = (2 * N - 3) as u8;

use std::process::exit;
use rayon::prelude::*;
use crate::algos::select;
use crate::cache::{Cache, compute_cache};
use crate::tgraph::{TGraph};

fn main() {
    let n: u8 = N as u8;
    println!("n={}", n);
    // let nb = generate(TGraph::new(n), &compute_cache());
    // let nb = count_all(TGraph::new(n));
    // let nb = check_spanners(TGraph::new(n), &compute_cache());
    // let nb = count_nondismountable(TGraph::new(n));
    // let nb = count_nonpivotable(TGraph::new(n));
    // let nb = count_nonboth(TGraph::new(n));
    // let nb = check_spanners_par(TGraph::new(n));
    let nb = generate_par(TGraph::new(n));
    println!("Nombre pour n={}: {}", n, nb);
}


///////////////// ITERATOR /////////////////////

pub struct TGraphs<'a> {
    stack: Vec<TGraph>,
    cache: &'a Cache,
    select: fn(&TGraph, &Cache) -> bool,
}

pub fn descendants(g: TGraph, select: Option<fn(&TGraph, &Cache) -> bool>, cache: &Cache) -> TGraphs {
    if select.is_none() {
        TGraphs { stack: vec![g], cache, select: |_, _| true }
    } else {
        TGraphs { stack: vec![g], cache, select: select.unwrap() }
    }
}

impl<'a> Iterator for TGraphs<'a>{
    type Item = TGraph;
    fn next(&mut self) -> Option<Self::Item> {
        while ! self.stack.is_empty() {
            let g = self.stack.pop().unwrap();
            if (self.select)(&g, &self.cache) {
                if g.has_symmetries() {
                    for h in g.successors_aut(&self.cache) {
                        self.stack.push(h);
                    }
                } else {
                    for h in g.successors_rigid(&self.cache) {
                        self.stack.push(h);
                    }
                }
                return Some(g);
            }
        }
        None
    }
}

pub fn count_nondismountable(g: TGraph) -> u64 {
    let mut nb = 0;
    let cache = compute_cache();
    for h in descendants(g, Some(|g, _| ! g.is_dismountable()), &cache){
        if h.is_clique() {
            nb += 1;
        }
    }
    nb
}

pub fn count_nonpivotable(g: TGraph) -> u64 {
    let mut nb = 0;
    let cache = compute_cache();
    for h in descendants(g, Some(|g, cache| ! g.has_pivot_vertex(cache)), &cache){
        if h.is_clique() {
            nb += 1;
        }
    }
    nb
}

pub fn count_nonboth(g: TGraph) -> u64 {
    let mut nb = 0;
    let cache = compute_cache();
    for h in descendants(g, Some(|g, cache| !g.has_pivot_vertex(cache) && !g.is_dismountable()), &cache){
        if h.is_clique() {
            nb += 1;
        }
    }
    nb
}

pub fn check_spanners(g: TGraph, cache: &Cache) -> u64 {
    let mut nb = 0;
    for h in descendants(g, Some(select), cache){
        if h.is_clique(){
            nb += 1;
            if !h.has_optimal_spanner(1000){
                println!("FAILING ON:");
                println!("{:?}", h.tedges());
                exit(0);
            }
        }
    }
    nb
}

////////////////////// GENERATION /////////////////////////

pub fn count_all(g: TGraph, cache: &Cache) -> u64{
    let mut nb = 0;
    for _ in descendants(g, None, &cache){
        nb += 1;
    }
    nb
}


///////////////// PARALLEL VERSIONS //////////////////////

// Splits the work into chunks to be done in parallel
pub fn get_pool(g: TGraph) -> (Vec<TGraph>, usize) {
    const TCUT: u8 = 5;
    let mut pool = vec![];
    let mut nbdropped = 0;
    let cache = compute_cache();
    for h in descendants(g, Some(|k, _| k.nb_edges <= TCUT), &cache).into_iter(){
        for hh in h.successors(&cache){
            if hh.nb_edges > TCUT{
                pool.push(hh);
            }
        }
        nbdropped += 1;
    }
    (pool, nbdropped)
}

pub fn generate_par(g: TGraph) -> usize {
    let cache: Cache = compute_cache();
    let pool = get_pool(g);
    let npool = pool.0;
    // let pool_size = npool.len();
    // println!("{} pieces in total", pool_size);
    let nb: u64 = npool.into_par_iter().enumerate()
        .map(|(_, h)| { // change _ with i for printing progress
            let mut res = 0;
            for _ in descendants(h, None, &cache){
                res += 1;
                // println!("{} / {}", pool_size, i + 1); // print progress in terminal
            }
            res
        })
        .sum();
    println!();
    return nb as usize + pool.1;
}

pub fn check_spanners_par(g: TGraph) -> usize {
    let cache: Cache = compute_cache();
    let pool = get_pool(g);
    let npool = pool.0;
    // let pool_size = npool.len();
    // println!("{} pieces in total", pool_size);
    let nb: u64 = npool.into_par_iter().enumerate()
        .map(|(_, h)| { // change _ with i for printing
            let res = check_spanners(h, &cache);
            // println!("{} / {}", pool_size, i + 1); // uncomment for progression in terminal
            res
        })
        .sum();
    println!();
    return nb as usize + pool.1;
}
