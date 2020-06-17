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

isclique(g) = length(g.nedges) == 0

genpairs(n) = [(i,j) for i::Int8 in 1:n-1 for j::Int8 in i+1:n]

# include("automorphisms.jl")

# Reference implementation for information, not used (see construct_from())
function construct_from_ref(g, new_edges::Vector{Tuple{Int8,Int8}}, t::Int8, rigid = isrigid(g))
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
function construct_from(g, new_edges::Vector{Tuple{Int8,Int8}}, t::Int8, rigid = isrigid(g))
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

function add_edges_new_time(g, edges::Vector{Tuple{Int8,Int8}}, t::Int8)
	g.tmax = t
	empty!(g.vmax)
	for (u, v) in edges
		push!(g.tedges, (u, v, t))
		push!(g.vmax, u, v)
		filter!(e->e≠(u, v), g.nedges)
	end
end

function add_edge(g, u, v, t)
    push!(g.tedges, (u, v, t))
	if t > g.tmax
		g.tmax = t
		empty!(g.vmax)
	end
	push!(g.vmax, u, v)
	filter!(e->e≠(u, v), g.nedges)
end

function edge_index(n::Int8, u::Int8, v::Int8)::Int8
	return (u - 1) * (n - (u / 2)) + (v - u)
end

function non_edges(g)
	return g.nedges
end

function isrigid(g)
	return g.rigid
end

function neighbors_dict(n, tedges)
	neighbors = [Dict{Int8,Int8}() for _ in 1:n]
	for (u,v,t) in tedges::Vector{Tuple{Int8,Int8,Int8}}
		neighbors[u][t] = v
		neighbors[v][t] = u
	end
	return neighbors
end
neighbors_dict(g) = neighbors_dict(g.n, g.tedges)

function get_components(n, tedges)
	# faster than union-find
	comps = [[u] for u::Int8 in 1:n]
	for (u, v, t) in tedges
		if comps[u] != comps[v]
			append!(comps[u],comps[v])
			for w in comps[v]
				comps[w] = comps[u]
			end
		end
	end
	# return unique!(comps) # Ref impl (slower than the following loop)
	final = Vector{Vector{Int8}}()
	for comp in comps
		if ! (comp in final)
			push!(final, comp)
		end
	end
	return final
end
get_components(g) = return get_components(g.n, g.tedges)


function has_isolated_vertex(g)
	isolated = fill(true, g.n)
	for (u, v, t) in g.tedges
		isolated[u] = false
		isolated[v] = false
	end
	return any(isolated)
end

function predecessors(g)
	preds = [Set{Int8}(u) for u in 1:g.n]
	for (u, v, _) in g.tedges
		union!(preds[u], preds[v])
		union!(preds[v], preds[u])
	end
	return preds
end

function is_tc(g)
    return all(length(p) == g.n for p in predecessors(g))
end
