include("tgraph.jl")
using DataStructures

function are_adjacent(e::Tuple{Int8,Int8}, f::Tuple{Int8,Int8})
	return e[1]==f[1] || e[1]==f[2] || e[2]==f[1] || e[2]==f[2]
end

function valid_subsets(edges::Vector{Tuple{Int8,Int8}})::Vector{Vector{Tuple{Int8,Int8}}}
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
	end
	head = popfirst!(edges)
	non_adjacent = filter(e -> !are_adjacent(head, e), edges)
	subsets = valid_subsets(non_adjacent)
	with_it = Vector{Vector{Tuple{Int8,Int8}}}(undef, length(subsets))
	head_tab = [head]
    for i in 1:length(subsets)
		with_it[i] = [head_tab; subsets[i]]
	end
    res = [with_it; valid_subsets(edges)]
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
	succ = Vector{TGraph}(undef, length(matchings))
	t = Int8(g.tmax + 1)
	for i in 1:length(matchings)
		h = construct_from(g, matchings[i], t, rigid)
		succ[i] = h
	end
	return succ
end

function Base.iterate(g::TGraph)
	stack = Stack{TGraph}()
	push!(stack, g)
	return iterate(g, stack)
end

function Base.iterate(g::TGraph, stack)
	if isempty(stack)
		return
	else
		h = pop!(stack)
		if length(h.nedges) > 0
	        for s in extensions(h)
				push!(stack, s)
			end
		end
		return (h, stack)
	end
end

using ResumableFunctions
@resumable function temporal_cliques(g::TGraph)
	stack = Stack{TGraph}()
	push!(stack, g)
    while !isempty(stack)
		h = pop!(stack)
        for s in extensions(h)
            if isempty(s.nedges)
            	@yield s
            else
				push!(stack, s)
			end
		end
	end
end

temporal_cliques(n::Int64) = temporal_cliques(TGraph(n))

function test_yield(n::Int64)
	clique_number = 0
	for g in Iterators.filter(isclique, TGraph(n)) # temporal_cliques(n)
		#if isclique(g)
			clique_number += 1
		#end
	end
	println(clique_number)
end
