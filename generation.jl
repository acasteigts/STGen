include("tgraph.jl")
include("automorphisms.jl")
using DataStructures

function are_adjacent(e, f)
	return e[1]==f[1] || e[1]==f[2] || e[2]==f[1] || e[2]==f[2]
end

function valid_subsets_small(edges)
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

function valid_subsets(edges)::Vector{Vector{Tuple{Int8,Int8}}}
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


function get_matchings_rigid(g)
	edges = filter(e -> e[1] in g.vmax || e[2] in g.vmax, non_edges(g))
	res = valid_subsets(edges)
	pop!(res)
	return res
end


function extend_matchings_aut(g, init_orbits, matchings)
	output = Vector{Tuple{Int8,Int8}}[]
	if isempty(matchings)
		# Nothing added so far: add first edge of every orbit (separately)
		for orbit in init_orbits
			e = orbit[1]
			if e in non_edges(g)
				if e[1] in g.vmax || e[2] in g.vmax
					push!(output,[e])
				end
			end
		end
		return output
	end

	for m in matchings
		t = Int8(g.tmax + 1)
		h = construct_from(g, m, t)
		orbits = edge_orbits_from_gens(h, h.gens)
		for orbit in orbits
			e = orbit[1]
			if e in non_edges(h)
				if e[1] in g.vmax || e[2] in g.vmax
					if ! (e[1] in h.vmax) && ! (e[2] in h.vmax)
						if orbit_of_edge(e, init_orbits) >= orbit_of_edge(m[end], init_orbits)
							push!(output, [m; [e]])
						end
					end
				end
			end
		end
	end
	return output
end

function get_matchings_aut(g, gens)
	matchings = Vector{Tuple{Int8,Int8}}[]
	nmatchings = Vector{Tuple{Int8,Int8}}[]
	init_orbits = edge_orbits_from_gens(g, gens)
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


function extensions(g) 
	if isrigid(g)
		matchings = get_matchings_rigid(g)
	else
		matchings = get_matchings_aut(g, g.gens)
	end

	exts = Vector{TGraph}(undef, length(matchings))
	t = Int8(g.tmax + 1)
	for i in 1:length(matchings)
		h = construct_from(g, matchings[i], t)
		exts[i] = h
	end
	return exts
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
		if !isclique(g)
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
