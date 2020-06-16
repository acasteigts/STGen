include("tgraph.jl")
include("generation.jl")
include("spanners.jl")


function count_cliques(root::TGraph)
	count = 0
	for g in TGraphs(root)
        if isclique(g)
			count += 1
		end
	end
	return count
end
count_cliques(n::Int64) = count_cliques(TGraph(n))

function count_graphs(root::TGraph)
	count = 0
	for g in TGraphs(root)
		count += 1
	end
	return count
end
count_graphs(n::Int64) = count_graphs(TGraph(n))

function count_symmetric(root::TGraph) # Should be 14 instead of 18 for n=4
	count = 0
	for g in TGraphs(root, g -> ! isrigid(g))
		count += 1
	end
	return count
end
count_symmetric(n::Int64) = count_symmetric(TGraph(n))

function check_spanners(root::TGraph)
	Random.seed!(0)
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
check_spanners(n::Int64) = check_spanners(TGraph(n))

# Generates a set of instances that covers all branches (used for parallelism)
# Each instance has at least depth edges
function branches(n, limit=Inf)
	branches = [TGraph(Int8(n))]
	count_skipped = 0
	if n < 5
		return branches, 0
	elseif n == 5
		depth = 3
	elseif n == 6
		depth = 5
	elseif n >= 7
		depth = 8
	end
	println("Creating sub-branches...")
	ok = false
	while !ok && length(branches) < limit
		ok = true
		for i in 1:length(branches)
			s = branches[i]
			if length(s.tedges) < depth
				splice!(branches, i, extensions(s))
				count_skipped += 1
				ok = false
				break
			end
		end
	end
	println(length(branches), " sub-branches created.")
	return shuffle!(branches), count_skipped
end



#########################################################################
# PARALLEL VERSIONS

using Distributed
using ProgressMeter

Random.seed!(0)
function count_graphs_par(n)
	bra, skipped = branches(n)
	results = @showprogress 1 "Computing..." pmap(count_graphs, bra)
	nb_graphs = sum(results) + skipped
	return nb_graphs
end

function count_cliques_par(n)
	bra, _ = branches(n)
	results = @showprogress 1 "Computing..." pmap(count_cliques, bra)
	nb_cliques = sum(results)
	return nb_cliques
end

function check_spanners_par(n)
	bra, _ = branches(n)
	results = @showprogress 1 "Computing..." pmap(check_spanners, bra)
	always_admit = all(results)
	return always_admit
end
