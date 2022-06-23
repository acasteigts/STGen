# STGen (A Generator for Simple Temporal Graphs)
A tool for generating simple temporal graphs up to isomorphism

## Background
In a temporal graph, the edges of the graph are available only at some times.
*Simple* temporal graphs (STGs) impose two restrictions:

1. Every edge is available at a *single* time
2. Two edges incident to a same vertex are available at *different* times

For instance, the following graphs are STGs (on 4 vertices):

```
o---1---o     o---1---o     o---2---o
|       |     |       |     |       |
2       4     2       2     3       5
|       |     |       |     |       |
o---5---o     o---3---o     o---1---o

   G1            G2            G3
```

A *temporal path* is a path that travels along increasing times.
(Non-decreasing times are allowed in some models, but they can't exist in STGs
by definition.)
If there is a temporal path from vertex *u* to vertex *v*, we say that u can reach v. A temporal graph is *temporally connected* is all vertices can reach
all vertices. For example, G3 is temporally connected, but G1 and G2 are not.

It often occurs that two different graphs have the same reachability.
For example, the graphs G1 and G2 above are equivalent in terms of reachability in the sense that their sets of temporal paths are the same up to time distortion. With a slight abuse of terminology, we say that G1 and G2 are *isomorphic*. In fact, we consider G2 as the *representative* of the isomorphism type of G1 and G2 (and infinitely many other graphs), because G2 is the graph with the smallest possible times that offers such a reachability.

The purpose of STGen is to generate *all* such representatives of a given size (which are in *finite* number), so that one can test various reachability conjectures exhaustively. For instance, an open question in [this paper]() asks whether one can always delete all but 2n-3 edges from a temporal clique (STG whose underlying graph is complete) without breaking temporal connectivity.

For more on the mathematics behind generation, in particular, the special properties of automorphisms in STGs, the user is referred to the following talk:

*[Efficient generation of simple temporal graphs up to isomorphism (YouTube)](https://www.youtube.com/watch?v=pgRBl--JJVc)*  
*Arnaud Casteigts, 3rd workshop on Algorithmic Aspects of Temporal Graphs (@ ICALP 2020)*
*[Longer version of the talk, in French](https://visio.u-bordeaux.fr/playback/presentation/2.0/playback.html?meetingId=bfe00d5046e9d24d0c256a9acfb841c176461c85-1599467221221)

This repository contains two versions: Julia and Rust.

## How to use STGen in Julia?

Here is a basic example that generates (and counts) STG representatives for a given number of vertices:

```Julia
include("generation.jl")

function count_graphs(n::Int)
    count = 0
    for g in TGraphs(n)
        count += 1
    end
    return count
end

julia> count_graphs(5)
15378
```

Technically, an STG is characterized by the number of vertices and a list of timed edges (triplet indicating the two endpoints and a time). The former can be accessed through `g.n` and the later through `g.edges`:

```Julia
julia>
for g in TGraphs(3)
    @show g.edges
end

g.tedges = Tuple{Int8,Int8,Int8}[]
g.tedges = Tuple{Int8,Int8,Int8}[(1, 2, 1)]
g.tedges = Tuple{Int8,Int8,Int8}[(1, 2, 1), (1, 3, 2)]
g.tedges = Tuple{Int8,Int8,Int8}[(1, 2, 1), (1, 3, 2), (2, 3, 3)]
```

The package contains a basic type called `TGraph` (the name for an STG in STGen) which offers only the minimal information used for generation, such as computing the group of automorphisms, testing for temporal connectivity, and adding new edges.

The original motivation for this generator was to test conjectures on temporal *cliques*, such as the one mentioned above. Temporal cliques can easily be filtered as:

```Julia
for g in TGraphs(n)
    if isclique(g)
	...
    end
end

julia> count_cliques(5)
4524
```

As shown below, the number of different STGs up to isomorphism explodes at an unreasonable rate:

| # Vertices   |      # STGs      |  # Temporally connected STGs |  # Simple Temporal cliques |
|:----------:|:-------------:|:------:|:------:|
| 1 |  1 | 1 | 1 |
| 2 |  2 | 1 | 1 |
| 3 |  4 | 1 | 1 |
| 4 | 62 | 32 | 20 |
| 5 | 15378 | 10207 | 4524 |
| 6 | 89769096 | 70557834 | 23218501 |
| 7 | 13828417028594 | ? | 3129434545680 |
| 8 | ? | ? | ? |

Therefore, in practice, one should not aim to brute force a conjecture naively above n=6 or perhaps 7.

### Make it parallel

The special features of the automorphism group of an STG make it possible to generate them as a tree whose branches are *independent*
(no comparisons are required between the branches). Thus, the code can easily be parallelized.
A look in the `main.jl` file will show you some ways to do so.
Essentially, the generator can take as parameter a graph rather than a number of vertices.
In this case, it generates only the subtree from that graph (the version with a number of vertices actually just calls this generator with an empty graph).
It is therefore sufficient to generate a bunch of graphs at small depth and map these graphs to different processors/threads as a pool of independent jobs.
Julia offers the `pmap` primitive (among others) to perform thus things.
Again, the interested user is referred to the `main.jl` class for an example.


## How to use STGen in Rust?

This version actually uses a lot of optimization that the Julia version does not, so it would be unfair to Julia to claim that the difference in performance is only due to the language.

Assuming that you have a valid version of rust and cargo installed.

Go into the rust directory, and type `cargo build --release` in order to compile an optimized binary of STGen. Then, type `cargo run --release`.

By default, this launches the parallel version of the generator, which enumerates (although only the count is shown) all STGs on 6 vertices up to isomorphism (in fact, reachability equivalence, see above). The main function in `main.rs` has a bunch of other function calls commented. You can execute any of them by uncommenting the corresponding line (instead of the default call to `generate_par`).

### Number of vertices?

For performance, the code uses no vectors for most of the processing (thus, no dynamic allocation). As such, the size of the arrays are determined at compilation time according to the number of vertices, which is to be set at the top of the `main.rs` file. This is not very convenient, but makes a big difference.