include("tgraph.jl")
using DataStructures

function are_adjacent(e::Tuple{Int8,Int8}, f::Tuple{Int8,Int8})
	return e[1]==f[1] || e[1]==f[2] || e[2]==f[1] || e[2]==f[2]
end

function valid_subsets_small(edges::Vector{Tuple{Int8,Int8}})::Vector{Vector{Tuple{Int8,Int8}}}
	if isempty(edges)
        return [Tuple{Int8, Int8}[]]
	elseif length(edges) == 1
    	return [[edges[1]], Tuple{Int8,Int8}[]]
	elseif length(edges) == 2
		if are_adjacent(edges[1], edges[2])
			return [[edges[1]], [edges[2]], Tuple{Int8,Int8}[]]
		else
			return [[edges[1]], [edges[2]], [edges[1], edges[2]], Tuple{Int8,Int8}[]]
		end
	elseif length(edges) == 3
	 	e1, e2, e3 = edges[1], edges[2], edges[3]
	 	if are_adjacent(e1, e2)
	 		if are_adjacent(e1, e3)
	 			if are_adjacent(e2, e3)
	 				return [[e1], [e2], [e3], Tuple{Int8,Int8}[]]
	 			else
	 				return [[e1], [e2], [e3], [e2, e3], Tuple{Int8,Int8}[]]
	 			end
	 		else
	 			if are_adjacent(e2, e3)
	 				return [[e1], [e2], [e3], [e1, e3], Tuple{Int8,Int8}[]]
	 			else
	 				return [[e1], [e2], [e3], [e1, e3], [e2, e3], Tuple{Int8,Int8}[]]
	 			end
	 		end
	 	else
	 		if are_adjacent(e1, e3)
	 			if are_adjacent(e2, e3)
	 				return [[e1], [e2], [e3], [e1, e2], Tuple{Int8,Int8}[]]
	 			else
	 				return [[e1], [e2], [e3], [e1, e2], [e2, e3], Tuple{Int8,Int8}[]]
	 			end
	 		else
	 			if are_adjacent(e2, e3)
	 				return [[e1], [e2], [e3], [e1, e2], [e1, e3], Tuple{Int8,Int8}[]]
	 			else
	 				return [[e1], [e2], [e3], [e1, e2], [e1, e3], [e2, e3], [e1, e2, e3], Tuple{Int8,Int8}[]]
	 			end
	 		end
	 	end
	end
end

function valid_subsets(edges::Vector{Tuple{Int8,Int8}})::Vector{Vector{Tuple{Int8,Int8}}}
	if length(edges) < 4
		return valid_subsets_small(edges)
	end
	head = popfirst!(edges)
	non_adjacent = filter(e -> !are_adjacent(head, e), edges)
	subsets = valid_subsets(non_adjacent)
	with_it = Vector{Vector{Tuple{Int8,Int8}}}(undef, length(subsets))
	head_tab = [head]
    for i in 1:length(subsets)
		with_it[i] = [head_tab; subsets[i]]
	end
	without_it = valid_subsets(edges)
    res = [with_it; without_it]
	return res
end


function get_matchings_rigid(g::TGraph)
	edges = filter(e -> e[1] in g.vmax || e[2] in g.vmax, g.nedges)
	res = valid_subsets(edges)
	pop!(res)
	return res
end


function get_matchings_aut(g::TGraph, gens)
	matchings = Vector{Tuple{Int8,Int8}}[]
	nmatchings = Vector{Tuple{Int8,Int8}}[]
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


function get_noncomponents(g::TGraph)::Vector{Vector{Int8}}
	# faster than union-find
	comps = [[u] for u::Int8 in 1:g.n]
	for (u, v) in g.nedges
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
		if length(comp) > 1 && ! (comp in final)
			push!(final, comp)
		end
	end
	return final
end

# At least one edge must be selected in each component
function filter_dead_end(m::Vector{Tuple{Int8, Int8}}, components::Vector{Vector{Int8}})
    if length(m) < length(components)
        return false
	end
    vertices = collect(Iterators.flatten(m))
	return all(!isempty(intersect(vertices, comp)) for comp in components)
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

    # The following is only relevant for generating cliques with n >= 8
    if g.n >= 8
        noncomponents = get_noncomponents(g)
		nb_noncomp = length(noncomponents)
        if nb_noncomp > 1
			filter!(m -> filter_dead_end(m, noncomponents), matchings)
		end
	end

	succ = Vector{TGraph}(undef, length(matchings))
	t = Int8(g.tmax + 1)
	for i in 1:length(matchings)
		h = construct_from(g, matchings[i], t, rigid)
		succ[i] = h
	end
	return succ
end


#######################################
# ITERATORS

struct TGraphs
	root::TGraph
	select_condition::Function
	TGraphs(root::TGraph, select_condition = (g) -> true) = new(root, select_condition)
	TGraphs(n::Int, select_condition = (g) -> true) = new(TGraph(n), select_condition)
end

function Base.iterate(graphs::TGraphs)
	stack = Stack{TGraph}()
	push!(stack, graphs.root)
	return iterate(graphs, stack)
end

function Base.iterate(graphs::TGraphs, stack)
	if isempty(stack)
		return
	else
		g = pop!(stack)
		if length(g.nedges) > 0
	        for s in extensions(g)
				if graphs.select_condition(s)
					push!(stack, s)
				end
			end
		end
		return (g, stack)
	end
end

function testt()
	c = 0
	for g in TGraphs(6)
		if isclique(g)
			c += 1
		end
	end
	println(c)
end
