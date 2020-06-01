include("tgraphs.jl")

function are_adjacent(e::Tuple{Int8,Int8}, f::Tuple{Int8,Int8})
	return e[1]==f[1] || e[1]==f[2] || e[2]==f[1] || e[2]==f[2]
end

function valid_subsets(lst::Vector{Tuple{Int8,Int8}})::Vector{Vector{Tuple{Int8,Int8}}}
    if length(lst) == 0
        return [Tuple{Int8, Int8}[]]
	end
    if length(lst) == 1
    	return [[lst[1]],Tuple{Int8,Int8}[]]
	end
	head = popfirst!(lst)
	non_adjacent = filter(e -> !are_adjacent(head, e), lst)
	subsets = valid_subsets(non_adjacent)
	with_it = Vector{Vector{Tuple{Int8,Int8}}}(undef, length(subsets))
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
	for i in 1:length(matchings)
		h = construct_from(g, matchings[i], rigid)
		succ[i] = h
	end
	return succ
end
