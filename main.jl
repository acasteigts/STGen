include("tgraph.jl")
include("generation.jl")
include("algorithms.jl")


function count_from(root::TGraph, predicate::Function = isclique)
	count = 0
	for g in root # syntax equiv. descendants(root)
        if predicate(g)
			count += 1
		end
	end
	return count
end

function exploreDFS_check(root::TGraph)
	stack = Stack{TGraph}()
	push!(stack, root)
	nb_cliques = 0
    while !isempty(stack)
		g = pop!(stack)
        for s in extensions(g)
			if select(s)
	            if isclique(s)
					if !has_optimal_spanner(s, 100)
	                	nb_cliques += 1
					end
	            else
					push!(stack, s)
				end
			end
		end
	end
	return nb_cliques
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

# Generates a set of instances that covers all branches (used for parallelism)
# Each instance has at least depth edges
function branches(n, limit=Inf)
	branches = [TGraph(Int8(n))]
	if n < 5
		return branches
	elseif n == 5
		depth = 3
	elseif n == 6
		depth = 5
	elseif n > 6
		depth = 8
	end
	println("Creating sub-branches...")
	ok = false
	while !ok
		ok = true
		for i in 1:length(branches)
			s = branches[i]
			if length(s.tedges) < depth
				splice!(branches, i, extensions(s))
				ok = false
				break
			end
		end
	end
	println(length(branches), " sub-branches created.")
	return shuffle!(branches)
end


#########################################################################
# SEQUENTIAL VERSION

function gen(n, check::Bool = false)
    g = TGraph(Int8(n))
	if check
		nb_cliques = exploreDFS_check(g)
	else
		nb_cliques = count_from(g, isclique)
	end
	println(nb_cliques, " cliques generated")
end


#########################################################################
# PARALLEL VERSION

using Distributed
using ProgressMeter

function gen_par(n, check::Bool = false)
	Random.seed!(0)
	bra = branches(n)
	if n == 7
		bra = bra[1:280]
	end
	if check
		results = @showprogress 1 "Computing..." pmap(exploreDFS_check, bra)
	else
		results = @showprogress 1 "Computing..." pmap(exploreDFS_nocheck, bra)
	end
	nb_cliques = sum(results)
	println(nb_cliques)
end
