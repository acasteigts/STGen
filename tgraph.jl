mutable struct TGraph
	n::Int8
	tmax::Int8
	tedges::Vector{Tuple{Int8,Int8,Int8}}
	nedges::Vector{Tuple{Int8,Int8}}
	vmax::Vector{Int8}
	rigid::Bool
	TGraph(n, tmax=0, tedges=[], nedges=genpairs(n), vmax=collect(1:n), rigid=false) = new(Int8(n), tmax, tedges, nedges, vmax, rigid)
	# TGraph(g::TGraph) = new(g.n, g.tmax, copy(g.tedges), copy(g.nedges), Int8[], g.rigid)
end

isclique(g::TGraph) = length(g.nedges) == 0

genpairs(n) = [(i,j) for i::Int8 in 1:n-1 for j::Int8 in i+1:n]

include("automorphisms.jl")

# Reference implementation for information, not used (see construct_from())
function construct_from_ref(g::TGraph, new_edges::Vector{Tuple{Int8,Int8}}, t::Int8, rigid = g.rigid)
	time_edges = copy(g.tedges)
	non_edges = copy(g.nedges)
	vmax = Int8[]
	for (u, v) in new_edges
		push!(time_edges, (u, v, t))
		push!(vmax, u, v)
		filter!(e->e≠(u, v), non_edges)
	end
	return TGraph(g.n, t, time_edges, non_edges, vmax, rigid)
end

# UGLY BUT FASTER
function construct_from(g::TGraph, new_edges::Vector{Tuple{Int8,Int8}}, t::Int8, rigid = g.rigid)
	m = length(g.tedges)
	k = length(new_edges)
	time_edges = Vector{Tuple{Int8,Int8,Int8}}(undef, m + k)
	for i in 1:m
		time_edges[i] = g.tedges[i]
	end
	non_edges = Vector{Tuple{Int8,Int8}}(undef, length(g.nedges) - k)
	j = 1
	for i in 1:length(g.nedges)
		(u, v) = g.nedges[i]
		if !((u, v) in new_edges)
			non_edges[j] = (u, v)
			j += 1
		end
	end
	vmax = Vector{Int8}(undef, k*2)
	for i in 1:k
		(u, v) = new_edges[i]
		time_edges[m + i] = (u, v, t)
		vmax[2*i-1] = u
		vmax[2*i] = v
	end
	return TGraph(g.n, t, time_edges, non_edges, vmax, rigid)
end

function add_edges_new_time(g::TGraph, edges::Vector{Tuple{Int8,Int8}}, t::Int8)
	g.tmax = t
	empty!(g.vmax)
	for (u, v) in edges
		push!(g.tedges, (u, v, t))
		push!(g.vmax, u, v)
		filter!(e->e≠(u, v), g.nedges)
	end
end

function add_edge(g::TGraph, u, v, t)
    push!(g.tedges, (u, v, t))
	if t > g.tmax
		g.tmax = t
		empty!(g.vmax)
	end
	push!(g.vmax, u, v)
	filter!(e->e≠(u, v), g.nedges)
end

function neighbors_dict(g::TGraph)
	neighbors = [Dict{Int8,Int8}() for _ in 1:g.n]
	for (u,v,t) in g.tedges::Vector{Tuple{Int8,Int8,Int8}}
		neighbors[u][t] = v
		neighbors[v][t] = u
	end
	return neighbors
end

function has_isolated_vertex(g::TGraph)
    for u in 1:g.n
		deg = 0
        for (v, w, t) in g.tedges
            if v == u || w == u
				deg += 1
			end
		end
		if deg == 0
			return true
		end
	end
    return false
end

function predecessors(g::TGraph)
	preds = [Set{Int8}(u) for u in 1:g.n]
	for (u, v, _) in g.tedges
		union!(preds[u], preds[v])
		union!(preds[v], preds[u])
	end
	return preds
end
