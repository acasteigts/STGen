
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
	final = Vector{Vector{Int8}}()
	for comp in comps
		if ! (comp in final)
			push!(final, comp)
		end
	end
	return final
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

function get_root(parents, i)::Int8
	res = parents[i] < 0 ? i : get_root(parents, parents[i])
	return res
end

function edge_image(u, v, perm)::Tuple{Int8, Int8}
	res = perm[u] < perm[v] ? (perm[u], perm[v]) : (perm[v], perm[u])
	return res
end

function edge_orbits(g::TGraph, gens)
    n = g.n
    m = Int8(n * (n - 1) / 2)
	epairs = genpairs(n)
    parents = fill(Int8(-1), m)

    for i in 1:length(epairs)
        if parents[i] < 0 # if root
            (u, v) = epairs[i]
            for perm in gens
                u2, v2 = u, v
                while true
                    u2, v2 = edge_image(u2, v2, perm)
                    if u2 == u && v2 == v
                        break
                    else
						ind = Int8(findfirst(x -> x==(u2,v2), epairs))
                        root2 = get_root(parents, ind)
                        if i != root2
							parents[i] += parents[root2]
							parents[root2] = i
						end
					end
				end
			end
		end
	end

    for i::Int8 in 1:m
        if parents[i] >= 0
            parents[i] = get_root(parents, i)
		end
	end
    for i::Int8 in 1:m
        if parents[i] < 0
            parents[i] = i
		end
	end
	sparents = Set(parents)
	res = [[epairs[i] for (i,e) in enumerate(parents) if e == r] for r in sparents]
	return res
end

function orbit_of_edge(e, orbits)
    for (i, orbit) in enumerate(orbits)
        if e in orbit
            return Int8(i)
		end
	end
end

function extend_matchings_aut(g::TGraph, init_orbits, matchings)
	output = Vector{Tuple{Int8,Int8}}[]
	if isempty(matchings)
		for orbit in init_orbits
			e = orbit[1]
			if e in g.nedges
				if e[1] in g.vmax || e[2] in g.vmax
					push!(output,[e])
				end
			end
		end
		return output
	end

	for m in matchings
		t = Int8(g.tmax + 1)
		h = construct_from(g, m, t, false)
		gens = automorphism_group(h)
		orbits = edge_orbits(h, gens)
		for orbit in orbits
			e = orbit[1]
			if e in g.nedges
				if e[1] in g.vmax || e[2] in g.vmax
					if ! (e[1] in h.vmax) && ! (e[2] in h.vmax)
						if orbit_of_edge(e, init_orbits) >= orbit_of_edge(m[end], init_orbits)
							push!(output, vcat(m, [e]))
						end
					end
				end
			end
		end
	end
	return output
end

function find_gens_intra_comp(neighbors, comp) # restricted to given component
    n = length(neighbors)
    perms = Vector{Vector{Int8}}() # generators for this component, trivially extended to V
    u, tail = Iterators.peel(comp)
    for v in tail
        perm = fill(Int8(-1), n)
        if find_candidate_aut(perm, neighbors, u, v)
            if -1 in perm
                for i::Int8 in 1:n
                    if perm[i] == -1
                        perm[i] = i
					end
				end
			end
            if is_automorphism(perm, neighbors)
                push!(perms, perm)
			end
		end
	end
    return perms
end

function find_gens_inter_comp(neighbors, comps)
    n = length(neighbors)
    perms = Vector{Vector{Int8}}()
    for i in 1:(length(comps)-1)
        u = comps[i][1]
        for j in (i + 1):length(comps)
            for v in comps[j]
                perm = fill(Int8(-1), n)
                if find_candidate_aut(perm, neighbors, u, v)
                    for i::Int8 in 1:n
                        if perm[i] != -1
                            perm[perm[i]] = i # swaps
						end
					end
                    for i::Int8 in 1:n
                        if perm[i] == -1
                            perm[i] = i
						end
					end
                    if is_automorphism(perm, neighbors)
                        push!(perms,perm)
					end
				end
			end
		end
	end
    return perms
end

# Faire un DFS parmi les voisins de u,
# à chaque étape,
function find_candidate_aut(perm, neighbors, u, v)
    if fit(neighbors, u, v)
		perm[u] = v
        for t in keys(neighbors[u])
            w = neighbors[u][t]
            if perm[w] == -1
                if ! find_candidate_aut(perm, neighbors, w, neighbors[v][t])
                    return false
				end
			end
		end
        return true
	end
    return false
end

function is_automorphism(perm, neighbors)
    for u::Int8 in 1:length(perm)
        u2 = perm[u]
        for t in keys(neighbors[u])
            v = neighbors[u][t]
            v2 = neighbors[u2][t]
            if perm[v] != v2
                return false
			end
		end
	end
    return true
end

function fit(neighbors, u, v)
    return Set(keys(neighbors[u])) == Set(keys(neighbors[v]))
end

function automorphism_group(g::TGraph)
    neighbors = neighbors_dict(g)
    components = get_components(g)
    gens = Vector{Vector{Int8}}()
    for comp in components
        append!(gens, find_gens_intra_comp(neighbors, comp))
	end
    if length(components) > 1
        append!(gens, find_gens_inter_comp(neighbors, components))
	end
    return gens
end
