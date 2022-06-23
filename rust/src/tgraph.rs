use itertools::{Either, Itertools};
use itertools::Either::Left;
use itertools::Either::Right;
use crate::helper::{all_ones, bits_to_indices, indices_to_bits};
use crate::{M, N};
use crate::cache::{Cache, get_indicators};
use crate::tedges::TEdge;

#[derive(Clone, Debug)]
pub struct Reachability {
    pub(crate) preds: [u8; N], // who reached this vertex
    pub(crate) preds2: [u8; N], // who reached this vertex after pivoting
}

impl Reachability {
    pub fn new(n: u8) -> Reachability {
        let mut preds = [0_u8; N];
        for i in 0..n as usize {
            preds[i] |= 1 << i;
        }
        Reachability{preds, preds2: [0_u8; N]}
    }

    pub fn update(&mut self, e: TEdge){
        let u = e.0 as usize;
        let v = e.1 as usize;

        if self.preds[u] != all_ones() {
            self.preds[u] |= self.preds[v];
            if self.preds[u] == all_ones() { // now u is reached by everybody
                self.preds2[u] |= 1 << u; // it can thus be reached from "sink u"
            }
        }
        if self.preds[v] != all_ones() {
            self.preds[v] |= self.preds[u];
            if self.preds[v] == all_ones() {
                self.preds2[v] |= 1 << v;
            }
        }
        self.preds2[u] |= self.preds2[v];
        self.preds2[v] = self.preds2[u];
    }
}

#[derive(Clone, Debug)]
pub struct Dismountability {
    pub(crate) degrees: [u8; N],
    pub(crate) missing_ng:[u8; N],
    pub(crate) mins: u8,
    pub(crate) maxs: u8,
}

impl Dismountability {
    pub fn new() -> Dismountability {
        let mut missing_ng:[u8; N] = [all_ones(); N];
        for u in 0..N{
            missing_ng[u] ^= 1 << u;
        }
        Dismountability {
            degrees: [0_u8; N], missing_ng,
            mins: 0, maxs: 0
        }
    }

    pub fn update(&mut self, e: TEdge) {
        let u = e.0 as usize;
        let v = e.1 as usize;

        self.degrees[u] += 1;
        self.degrees[v] += 1;
        self.missing_ng[u] ^= 1 << v;
        self.missing_ng[v] ^= 1 << u;

        if self.degrees[u] == 1{
            self.mins |= 1 << v;
        }
        if self.degrees[v] == 1{
            self.mins |= 1 << u;
        }
        if self.degrees[u] == (N - 1) as u8 {
            self.maxs |= 1 << v;
        }
        if self.degrees[v] == (N - 1) as u8 {
            self.maxs |= 1 << u;
        }
    }
}


#[derive(Clone, Debug)]
pub struct TGraph {
    pub(crate) n: u8,
    tmax: u8, // largest timestamp
    pub(crate) times: [u8; M],
    pub(crate) edges: [TEdge; M],
    pub(crate) edges_bits: u32,
    pub(crate) nb_edges: u8,
    pub(crate) cand_bits: u32,
    pub(crate) nb_cand_edges: usize,
    pub(crate) reachability: Reachability,
    pub(crate) dismountability: Dismountability,
    pub(crate) gens: Option<Vec<Vec<u8>>>,
}

impl TGraph {
    pub fn new(n: u8) -> TGraph {
        TGraph {n, tmax: 0, times: [0; M], edges: [TEdge(0,0,0); M], edges_bits: 0,
            nb_edges: 0, cand_bits: !0, nb_cand_edges: M,
            reachability: Reachability::new(n),
            dismountability: Dismountability::new(),
            gens: Some((0..n).permutations(n as usize).collect())}
    }
    pub fn tmax(&self) -> u8{
        self.tmax
    }
    pub fn has_symmetries(&self) -> bool{
        !self.gens.is_none()
    }
    pub fn extends_by(&self, indicators: &u32, cache: &Cache) -> TGraph {
        let mut times = self.times.clone();
        let mut edges = self.edges.clone();
        let mut reachability = self.reachability.clone();
        let mut dismountability = self.dismountability.clone();
        let tmax = self.tmax() + 1;
        let mut edges_bits = self.edges_bits;
        let mut cand_bits: u32 = 0;
        let mut nb_edges = self.nb_edges;
        for i in 0..M{
            if indicators & 1 << i != 0{
                times[i] = tmax;
                edges_bits |= 1 << i;
                cand_bits |= cache.adjacent_bits[i];
                let mut ne = cache.edges[i];
                ne.2 = tmax;
                edges[nb_edges as usize] = ne;
                nb_edges += 1;
                reachability.update(ne);
                dismountability.update(ne);
            }
        }
        cand_bits &= !edges_bits;
        let nb_cand_edges: usize = cand_bits.count_ones() as usize;
        if ! self.gens.is_none() {
            let gens:Vec<Vec<u8>> = self.gens.as_ref().unwrap().iter().filter(|a| {
                is_automorphism(&times, a, cache)
            }).map(|a| a.clone()).collect();
            if gens.len() > 1{
                return TGraph { n: self.n, tmax, times, edges, edges_bits,
                    nb_edges, cand_bits, nb_cand_edges,
                    reachability,
                    dismountability,
                    gens: Some(gens)};
            }
        }
        TGraph {n: self.n, tmax, times, edges, edges_bits,
            nb_edges, cand_bits, nb_cand_edges,
            reachability,
            dismountability,
            gens: None}
    }

    pub fn get_matchings(&self, cache: &Cache) -> Vec<Vec<usize>> {
        let mut matchings: Vec<Vec<usize>> = vec![];
        for bits in get_indicators(self.nb_cand_edges, cache) {
            if self.cand_bits | bits == self.cand_bits {
                let matching = bits_to_indices(*bits);
                if !contains_same_matching_up_to_automorphisms(self, &matchings, &matching, cache) {
                    matchings.push(matching);
                }
            }
        }
        matchings
    }

    pub fn successors_rigid<'a>(&'a self, cache: &'a Cache) -> impl Iterator<Item=TGraph> + 'a {
            get_indicators(self.nb_cand_edges, cache).into_iter()
                .filter(|&&bits| { self.cand_bits | bits == self.cand_bits })
                .map(|bits| { self.extends_by(&bits, cache) })
    }

    pub fn successors_aut<'a>(&'a self, cache: &'a Cache) -> impl Iterator<Item=TGraph> + 'a {
        self.get_matchings(cache).into_iter()
            .map(|m| self.extends_by(&indices_to_bits(&m), cache) )
    }

    pub fn successors<'a>(&'a self, cache: &'a Cache) -> Either<impl Iterator<Item=TGraph> + 'a, impl Iterator<Item=TGraph> + 'a> {
        if self.has_symmetries() {
            Left(self.successors_aut(cache))
        } else {
            Right(self.successors_rigid(cache))
        }
    }

    pub fn tedges(&self) -> &[TEdge] {
        return &self.edges[0..self.nb_edges as usize];
    }

    pub fn predecessors(&self) -> [u8; N] {
        self.reachability.preds
        // let mut preds = [0_32; N];
        // for i in 0..self.n as usize {
        //     preds[i] |= 1 << i;
        // }
        //
        // for e in self.tedges(){
        //     preds[e.0 as usize] |= preds[e.1 as usize];
        //     preds[e.1 as usize] = preds[e.0 as usize];
        // }
        // preds
    }
}

pub fn is_automorphism(times: &[u8; M], p: &[u8], cache: &Cache) -> bool {
    for (i, t) in times.iter().enumerate().rev() {
        let e = cache.edges[i];
        let i = cache.e2i[p[e.0 as usize] as usize][p[e.1 as usize] as usize];
        if times[i] != *t {
            return false;
        }
    }
    true
}

pub fn same_matching_by_perm(vec1 : &[usize], vec2 : &[usize], perm : &[u8], cache: &Cache) -> bool {
    for i in vec1.iter(){
        let u2: usize = perm[cache.edges[*i].0 as usize] as usize;
        let v2: usize = perm[cache.edges[*i].1 as usize] as usize;
        let i2 = cache.e2i[u2][v2];
        if !vec2.contains(&i2){
            return false;
        }
    }
    true
}

pub fn contains_same_matching_up_to_automorphisms(g: &TGraph, matchings: &Vec<Vec<usize>>, m: &Vec<usize>, cache: &Cache) -> bool {
    for m2 in matchings.iter().rev() {
        if m.len() == m2.len() {
            for perm in g.gens.as_ref().unwrap().iter() {
                if same_matching_by_perm(m, m2, perm, cache) {
                    return true;
                }
            }
        }
    }
    return false;
}

