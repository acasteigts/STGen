include("tgraph.jl")
include("generation.jl")
include("algorithms.jl")


function count_cliques(root::TGraph)
	count = 0
	for g in TGraphs(root)
        if isclique(g)
			count += 1
		end
	end
	return count
end

function check_spanners(root::TGraph)
	for g in TGraphs(root, g -> select(g))
        if isclique(g)
			if !has_optimal_spanner(g, 1000)
		        print("FAILING ON ")
		        print(g.tedges)
				return false
		    end
		end
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
    root = TGraph(Int8(n))
	if check
		always_admit = check_spanners(root)
		@show always_admit
	else
		nb_cliques = count_cliques(root)
		@show nb_cliques
	end
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
		results = @showprogress 1 "Computing..." pmap(check_spanners, bra)
		always_admit = all(results)
		@show always_admit
	else
		results = @showprogress 1 "Computing..." pmap(count_cliques, bra)
		nb_cliques = sum(results)
		@show nb_cliques
	end
end
