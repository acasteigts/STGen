try
	g = TGraph(Int8(4))
catch
	include("tgraph.jl")
end
include("algorithms.jl")

using DataStructures
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
#			return iterate(h, stack)
		end
		return (h, stack)
	end
end

function test_iter(n::Int)
	nb_graphs = 0
	for g in TGraph(Int8(n))
		println(g)
		break;
		# if isempty(g.nedges)
			nb_graphs += 1
		# end
	end
	println(nb_graphs)
end

function exploreDFS_nocheck(g::TGraph)
	stack = Stack{TGraph}()
	push!(stack, g)
	nb_cliques = 0
    while !isempty(stack)
		h = pop!(stack)
        for s in extensions(h)
            if isempty(s.nedges)
            	nb_cliques += 1
            else
				push!(stack, s)
			end
		end
	end
	return nb_cliques
end

function exploreDFS_check(g::TGraph)
	stack = Stack{TGraph}()
	push!(stack, g)
	nb_cliques = 0
    while !isempty(stack)
		h = pop!(stack)
        for s in extensions(h)
			if select(s)
	            if isempty(s.nedges)
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
		nb_cliques = exploreDFS_nocheck(g)
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
