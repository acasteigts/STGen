
function get_root(parents, i)::Int8
	res = parents[i] < 0 ? i : get_root(parents, parents[i])
	return res
end

function edge_image(u, v, perm)::Tuple{Int8, Int8}
	u2 = perm[u]
	v2 = perm[v]
	return min((u2, v2), (v2, u2))
end

function edge_orbits_from_gens(g::TGraph, gens)
    n = g.n
    m = Int8(n * (n - 1) / 2)
	edges = genpairs(n)
    parents = fill(Int8(-1), m)

    for (i, (u, v)) in enumerate(edges)
        if parents[i] < 0 # if root
            for perm in gens
                u2, v2 = u, v
                while true
                    u2, v2 = edge_image(u2, v2, perm)
                    if u2 == u && v2 == v
                        break
                    else
						ind = edge_index(n, u2, v2)
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
	res = [[edges[i] for (i,e) in enumerate(parents) if e == r] for r in sparents]
	# return filter!(e -> e[1] in non_edges(g), res)
	return res
end

function orbit_of_edge(e, orbits)
    for (i, orbit) in enumerate(orbits)
        if e in orbit
            return Int8(i)
		end
	end
end



function find_gens_intra_comp(neighbors, comp)
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

function automorphism_group(n, tedges)
    neighbors = neighbors_dict(n, tedges)
    components = get_components(n, tedges)
    gens = Vector{Vector{Int8}}()
    for comp in components
        append!(gens, find_gens_intra_comp(neighbors, comp))
	end
    if length(components) > 1
        append!(gens, find_gens_inter_comp(neighbors, components))
	end
    return gens
end
automorphism_group(g) = automorphism_group(g.n, g.tedges)
