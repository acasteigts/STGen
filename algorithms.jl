function is_dismountable(g::TGraph)
    n = g.n
    vertices = 1:n
    is_min = [false for i in 1:n]
    is_max = [false for i in 1:n]

    missing_max = Int8[]

    neighbors = neighbors_dict(g)
    for v in vertices
        ng_v = neighbors[v]
        degree = length(ng_v)
        if degree > 0
            min_ng = ng_v[minimum(keys(ng_v))]
            is_min[min_ng] = true
        end
        if degree == n - 1
            max_ng = ng_v[maximum(keys(ng_v))]
            is_max[max_ng] = true
        else
            push!(missing_max, v)
        end
    end

    for v in vertices
        if is_min[v] && is_max[v]
            return true
        end
    end

    for v in missing_max
        missing_ng = Set{Int8}(vertices)
        setdiff!(missing_ng, values(neighbors[v]))
        setdiff!(missing_ng, v)
        broke = false
        for u in missing_ng  # if all missing are min, one is dismountable
            if !is_min[u]
                broke = true
                break;
            end
        end
        if !broke
            return true
        end
    end
    return false
end




######## PIVOTING ###

function isolated_trick(g)
    if length(g.vmax) == 2
        u, v, _ = g.tedges[end]
        preds = predecessors(g)
        if length(preds[u]) == g.n - 1 || length(preds[v]) == g.n - 1
            return true
        end
    end
    return false
end


function has_pivot_vertex(g)
    if has_isolated_vertex(g)
        return isolated_trick(g)
    end

    n = g.n
    preds = [Set{Int8}(i) for i in 1:n] # who reached this vertex
    preds2 = [Set{Int8}() for _ in 1:n] # who reached this vertex after pivoting
    fpreds2 = [Set{Int8}() for _ in 1:n] # completion with non-edges

    # if g.is_incremental():
    #     preds, preds2 = g.predecessors(), g.predecessors2()
    # else:
    if true
        for (u, v, t) in g.tedges
            if length(preds[u]) < n # not reached yet
                union!(preds[u], preds[v])
                if length(preds[u]) == n
                    push!(preds2[u], u)
                end
            end
            if length(preds[v]) < n
                union!(preds[v], preds[u])
                if length(preds[v]) == n
                    push!(preds2[v], v)
                end
            end

            union!(preds2[u], preds2[v])
            union!(preds2[v], preds2[u])
        end
    end

    for u in 1:n
        union!(fpreds2[u], preds2[u])
    end

    for (u, v) in g.nedges
        union!(fpreds2[u], preds2[v])
        union!(fpreds2[v], preds2[u])
    end

    inter = intersect(fpreds2...)
    return length(inter) > 0
end




########## RANDOM SPANNER ###
using Random
Random.seed!(0)

function is_tc(g)
    return all(length(p) == g.n for p in predecessors(g))
end

function is_tc_without(n, edges, e)
    if length(edges) < 2*n - 3
        return false
    end
    preds = [Set{Int8}(u) for u in 1:n]
    for (u,v,t) in edges
        if (u,v,t) != e
            union!(preds[u], preds[v])
            empty!(preds[v]) # speeds things up
            union!(preds[v], preds[u])
        end
    end
    return all(length(p) == n for p in preds)
end

function random_minimal_spanner(g)
    rand_edges = shuffle(g.tedges)
    spanner = copy(g.tedges)
    for e in rand_edges
        if is_tc_without(g.n, spanner, e)
            filter!(edge -> edge != e, spanner)
        end
    end
    return spanner
end

function has_optimal_spanner(g, nb_try = 1)
    if length(g.tedges) < (2*g.n - 4) || !is_tc(g)
        return false
    end

    i = 1
    while true
        i += 1
        spanner = random_minimal_spanner(g)
        if length(spanner) <= 2 * g.n - 3
            return true
        end
        if i > nb_try
            return false
        end
    end
    return true
end


function brute_force(g)
    if !has_optimal_spanner(g, 1000)
        print("FAILING ON ")
        print(g.tedges)
        readline()
        return false
    end
    return true
end

function select(g::TGraph)
	if is_dismountable(g)
		return false
	elseif has_pivot_vertex(g)
		return false
	elseif has_optimal_spanner(g, 1)
		return false
	else
		return true
	end
end
