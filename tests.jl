include("main.jl")

function dowork(n::Int, results::Channel)
	return n*2
end

import Base.Threads.@spawn
function partest()
	results = Channel(100)
	for i in 1:100
		res = @spawn dowork(results, i)
	end
	total = 0
	for i in 1:100
		n = take!(results)
		total += n
	end
	println(total)
end

function testt()
	g = TGraph(Int8(4))
	add_edge(g, Int8(1), Int8(2), Int8(10))
	gens = get_components(g) # ext = extensions(g)
	return nothing
end

function test2()
	a = Array{Int64, 1}([1,4,5])
	for v in a
		# anything
	end
end


function test_iter(n::Int)
	nb_graphs = 0
	for g in TGraph(n)
		if isclique(g)
			nb_graphs += 1
		end
	end
	println(nb_graphs)
end
