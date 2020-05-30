mutable struct TGraph
	n::Int8
	tmax::Int8
	tedges::Array{Tuple{Int8,Int8,Int8},1}
	nedges::Array{Tuple{Int8,Int8},1}
	vmax::Array{Int8,1}
	rigid::Bool
	TGraph(n::Int8) = new(n, 0, [], genpairs(n), collect(1:n), false)
	TGraph(g::TGraph) = new(g.n, g.tmax, copy(g.tedges), copy(g.nedges), Int8[], g.rigid)
end

genpairs(n::Int8) = [(i,j) for i::Int8 in 1:n-1 for j::Int8 in i+1:n]

include("automorphisms.jl")

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
	comp = [[u] for u::Int8 in 1:g.n]
	for (u, v, t) in g.tedges
		if comp[u] != comp[v]
			append!(comp[u],comp[v])
			for w in comp[v]
				comp[w] = comp[u]
			end
		end
	end
	final = Array{Array{Int8,1},1}()
	for c1 in comp
		if ! (c1 in final)
			push!(final,c1)
		end
	end
	return final
end

function are_adjacent(e::Tuple{Int8,Int8}, f::Tuple{Int8,Int8})
	return e[1]==f[1] || e[1]==f[2] || e[2]==f[1] || e[2]==f[2]
end

function valid_subsets(lst::Array{Tuple{Int8,Int8},1})::Array{Array{Tuple{Int8,Int8}, 1}, 1}
    if length(lst) == 0
		res = Array{Tuple{Int8,Int8}, 1}[]
		push!(res, Tuple{Int8,Int8}[])
        return res
	end
    if length(lst) == 1
		res = [Tuple{Int8,Int8}[lst[1]],Tuple{Int8,Int8}[]]
    	return res
	end
	head = popfirst!(lst)
    with_it = Array{Tuple{Int8,Int8}, 1}[]
	non_adjacent = [e for e in lst if !are_adjacent(head, e)]
	subsets = valid_subsets(non_adjacent)
    for s in subsets
		append!(with_it, [vcat([head], s)])
	end
    res = vcat(with_it, valid_subsets(lst))
	return res
end

function get_matchings_rigid(g::TGraph)
	vt = g.vmax
	edges = [(u, v) for (u, v) in g.nedges if u in vt || v in vt]
	res = valid_subsets(edges)
	splice!(res,length(res))
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
	succ = TGraph[]
	for matching in matchings
		h = TGraph(g)
		h.rigid = rigid
		# for (u, v) in matching
		# 	add_edge(h, u, v, t)
		# end
		add_edges_new_time(h, matching, t)
		push!(succ, h)
	end
	return succ
end
