use std::collections::{HashSet};
use crate::{Cache, DNMT, M, N};
use crate::tgraph::{TGraph};
use rand::{thread_rng};
use rand::seq::SliceRandom;
use crate::helper::{all_ones};

pub fn select(g: &TGraph, cache: &Cache) -> bool{
    ! (g.is_dismountable() || g.has_pivot_vertex(cache) || g.has_optimal_spanner(2))
}

impl TGraph {
    // pub fn get_components(&self) -> Vec<u8>{
    //     // faster than union-find
    //     let mut comps: [u8; N] = [0; N];
    //     for i in 0..self.n as usize{
    //         comps[i] = 1 << i;
    //     }
    //     for e in self.tedges() {
    //         let u = e.0 as usize;
    //         let v = e.1 as usize;
    //         if comps[u] != comps[v] {
    //             comps[u] |= comps[v];
    //             for i in 0..self.n as usize{
    //                 if comps[u] & 1 << i != 0 {
    //                     comps[i] = comps[u];
    //                 }
    //             }
    //         }
    //     }
    //     return comps.into_iter().unique().collect();
    // }

    pub fn is_clique(&self) -> bool {
        self.nb_edges == (M as u8)
    }



    /////////////////////////////////////////////// DISMOUNTING

    pub fn is_dismountable(&self) -> bool {
        let mins = self.dismountability.mins;
        let maxs = self.dismountability.maxs;
        if mins & maxs != 0 {
            return true;
        }

        let missing_ng = self.dismountability.missing_ng;
        for u in 0..N{
            if missing_ng[u]!= 0 && missing_ng[u] & mins == missing_ng[u]{
                return true;
            }
        }
        false
    }



    /////////////////////////////////////////////// PIVOTING

    pub fn is_tc(&self) -> bool {
        self.predecessors().iter().all(|p| *p == all_ones())
    }

    fn is_tc_without(&self, removed: u32, cand_remove: usize) -> bool {
        let mut preds = [0; N];
        for i in 0..self.n as usize {
            preds[i] |= 1 << i;
        }
        for (i, e) in self.tedges().iter().enumerate() {
            if i != cand_remove && removed & 1 << i == 0 {
                preds[e.0 as usize] |= preds[e.1 as usize];
                preds[e.1 as usize] = preds[e.0 as usize];
            }
        }
        preds.iter().all(|p| *p == all_ones())
    }

    pub fn has_isolated_vertex(&self) -> bool {
        self.dismountability.degrees.into_iter().any(|deg| deg == 0)
    }

    pub fn has_pivot_vertex(&self, cache: &Cache) -> bool {
        let mut fpreds2 = self.reachability.preds2.clone();

        // For each non-edge, the two vertices will eventually merge their full predecessors
        // (eventually = because clique)
        for i in 0..M {
            if self.edges_bits & 1 << i == 0{
                let te = cache.edges[i];
                let u = te.0 as usize;
                let v = te.1 as usize;
                fpreds2[u] |= self.reachability.preds2[v];
                fpreds2[v] |= self.reachability.preds2[u];
            }
        }

        let inter = fpreds2.iter().fold(all_ones(), |res, val| res & *val);
        inter != 0
    }

    pub fn greedy_spanner_size(&self) -> u8 {
        let mut removed = 0_u32;
        let mut nb_removed = 0;
        let target_removed = self.tedges().len() as u8 - DNMT;
        for i in 0..self.nb_edges {
            if self.is_tc_without(removed, i as usize) {
                removed |= 1 << i;
                nb_removed += 1;
                if nb_removed == target_removed{
                    break;
                }
            }
        }
        return self.tedges().len() as u8 - nb_removed;
    }

    pub fn random_spanner_size(&self) -> u8 {
        let mut rand_indices: Vec<usize> = (0..self.nb_edges as usize).collect();
        rand_indices.shuffle(&mut thread_rng());
        let mut removed = 0_u32;
        let mut nb_removed = 0;
        let target_removed = self.tedges().len() as u8 - DNMT;
        // for i in (0..self.nb_edges).map(|_| thread_rng().gen_range(0..self.nb_edges)) {
        for i in rand_indices {
            if self.is_tc_without(removed, i as usize) {
                removed |= 1 << i;
                nb_removed += 1;
                if nb_removed == target_removed{
                    break;
                }
            }
        }
        return self.tedges().len() as u8 - nb_removed;
    }

    // Here optimal means 2n-3 or 2n-4 edges
    pub fn has_optimal_spanner(&self, nb_try: u32) -> bool {
        if self.nb_edges < (2 * self.n - 4) || !self.is_tc() {
            return false;
        }

        if self.greedy_spanner_size() <= DNMT{
            return true;
        }

        let mut i = 1;
        while i < nb_try {
            if self.random_spanner_size() <= DNMT {
                return true;
            }
            i += 1;
        }
        return false;
    }
}

/////////// HELPER ////////////////////////////

pub fn range_set_vec() -> Vec<HashSet<u8>> {
    let mut preds: Vec<HashSet<u8>> = Vec::with_capacity(N);
    for i in 0..N as u8{
        let mut set: HashSet<u8> = HashSet::with_capacity(N);
        set.insert(i);
        preds.push(set);
    }
    preds
}
pub fn empty_set_vec() -> Vec<HashSet<u8>> {
    let mut preds: Vec<HashSet<u8>> = Vec::with_capacity(N);
    for _ in 0..N as u8{
        preds.push(HashSet::with_capacity(N));
    }
    preds
}

