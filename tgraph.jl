struct TGraph
	n::Int8
	tmax::Int8
	tedges::Array{Tuple{Int8,Int8,Int8},1}
	nedges::Array{Tuple{Int8,Int8},1}
	vmax::Array{Int8,1}
	rigid::Bool
	TGraph(n::Int8, tmax=0, tedges=[], nedges=genpairs(n), vmax=collect(1:n), rigid=false) = new(n, tmax, tedges, nedges, vmax, rigid)
	# TGraph(g::TGraph) = new(g.n, g.tmax, copy(g.tedges), copy(g.nedges), Int8[], g.rigid)
end

genpairs(n::Int8) = [(i,j) for i::Int8 in 1:n-1 for j::Int8 in i+1:n]

include("automorphisms.jl")

# Reference implementation for information, not used (see construct_from())
function construct_from_ref(g::TGraph, new_edges::Array{Tuple{Int8,Int8},1}, t::Int8, rigid = g.rigid)
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
function construct_from(g::TGraph, new_edges::Array{Tuple{Int8,Int8},1}, t::Int8, rigid = g.rigid)
	m = length(g.tedges)
	k = length(new_edges)
	time_edges = Array{Tuple{Int8,Int8,Int8}, 1}(undef, m + k)
	for i in 1:m
		time_edges[i] = g.tedges[i]
	end
	non_edges = Array{Tuple{Int8,Int8}, 1}(undef, length(g.nedges) - k)
	j = 1
	for i in 1:length(g.nedges)
		(u, v) = g.nedges[i]
		if !((u, v) in new_edges)
			non_edges[j] = (u, v)
			j += 1
		end
	end
	vmax = Array{Int8, 1}(undef, k*2)
	for i in 1:k
		(u, v) = new_edges[i]
		time_edges[m + i] = (u, v, t)
		vmax[2*i-1] = u
		vmax[2*i] = v
	end
	return TGraph(g.n, t, time_edges, non_edges, vmax, rigid)
end

function add_edges_new_time(g::TGraph, edges::Array{Tuple{Int8,Int8},1}, t::Int8)
	g.tmax = t
	empty!(g.vmax)
	for (u, v) in edges
		push!(g.tedges, (u, v, t))
		push!(g.vmax, u, v)
		filter!(e->e≠(u, v), g.nedges)
	end
end

# function add_edge(g::TGraph, u::Int8, v::Int8, t::Int8)
#     push!(g.tedges, (u, v, t))
# 	if t > g.tmax
# 		g.tmax = t
# 		empty!(g.vmax)
# 	end
# 	push!(g.vmax, u, v)
# 	filter!(e->e≠(u, v), g.nedges)
# end

function neighbors_dict(g::TGraph)
	neighbors = [Dict{Int8,Int8}() for _ in 1:g.n]
	for (u,v,t) in g.tedges::Array{Tuple{Int8,Int8,Int8},1}
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

function get_components(g::TGraph)
	# faster than union-find
	comps = [[u] for u::Int8 in 1:g.n]
	for (u, v, t) in g.tedges
		if comps[u] != comps[v]
			append!(comps[u],comps[v])
			for w in comps[v]
				comps[w] = comps[u]
			end
		end
	end
	# return unique!(comps) # Ref impl (slower than the following loop)
	final = Array{Array{Int8,1},1}()
	for comp in comps
		if ! (comp in final)
			push!(final, comp)
		end
	end
	return final
end

function are_adjacent(e::Tuple{Int8,Int8}, f::Tuple{Int8,Int8})
	return e[1]==f[1] || e[1]==f[2] || e[2]==f[1] || e[2]==f[2]
end

function valid_subsets(lst::Array{Tuple{Int8,Int8},1})::Array{Array{Tuple{Int8,Int8}, 1}, 1}
    if length(lst) == 0
        return [Tuple{Int8, Int8}[]]
	end
    if length(lst) == 1
    	return [Tuple{Int8,Int8}[lst[1]],Tuple{Int8,Int8}[]]
	end
	head = popfirst!(lst)
	non_adjacent = [e for e in lst if !are_adjacent(head, e)]
	subsets = valid_subsets(non_adjacent)
	with_it = Array{Array{Tuple{Int8,Int8}, 1}, 1}(undef, length(subsets))
	head_tab = [head]
    for i in 1:length(subsets)
		with_it[i] = [head_tab; subsets[i]]
	end
    res = [with_it; valid_subsets(lst)]
	return res
end

function get_matchings_rigid(g::TGraph)
	edges = filter(e -> e[1] in g.vmax || e[2] in g.vmax, g.nedges)
	res = valid_subsets(edges)
	pop!(res)
	return res
end


function get_matchings_aut(g::TGraph, gens)
	matchings = Array{Tuple{Int8,Int8}, 1}[]
	nmatchings = Array{Tuple{Int8,Int8}, 1}[]
	init_orbits = edge_orbits(g, gens)
	for i in 1:Int8(floor(g.n / 2))
		nmatchings = extend_matchings_aut(g, init_orbits, nmatchings)
		if !isempty(nmatchings)
			append!(matchings, nmatchings)
		else
			break
		end
	end
	return matchings
end


function extensions(g::TGraph)
	rigid = g.rigid
	if rigid
		matchings = get_matchings_rigid(g)
	else
		gens = automorphism_group(g)
		rigid = isempty(gens)
		if rigid
			matchings = get_matchings_rigid(g)
		else
			matchings = get_matchings_aut(g, gens)
		end
	end
	t = Int8(g.tmax + 1)
	succ = Array{TGraph, 1}(undef, length(matchings))
	for i in 1:length(matchings)
		h = construct_from(g, matchings[i], t, rigid)
		succ[i] = h
	end
	return succ
end
