function serialize(g::TGraph)
	edges = copy(g.tedges)
	for (u, v) in g.nedges
		push!(edges, (u, v, 0))
	end
	sort!(edges)
	base = Int((g.n * (g.n-1)) / 2) + 1
	result = 0;
	for (i, (u, v, t)) in enumerate(edges)
		result += t * (base^i)
	end
	return string(result)
end


function unserialize(val::String, n)
	base = Int128((n * (n-1)) / 2) + 1
	tedges = Vector{Tuple{Int8, Int8, Int8}}()
	nedges = Vector{Tuple{Int8, Int8}}()
	tmax = 0
	valint = parse(Int, val)
	for (i, (u, v)) in enumerate(genpairs(n))
		b = base^(i+1)
		mod = valint % b
		t = floor(mod / (base^i))
		if t == 0
			push!(nedges, (u, v))
		else
			push!(tedges, (u, v, t))
		end
		if t > tmax
			tmax = t
		end
	end
	sort!(tedges, by=e -> e[3])
	return TGraph(n, tmax, tedges, nedges)
end

function test_serialize()
	g = TGraph(5)
	add_edge(g, 1, 3, 1)
	add_edge(g, 1, 2, 2)
	add_edge(g, 2, 4, 3)
	add_edge(g, 4, 5, 5)
	add_edge(g, 3, 5, 5)
	add_edge(g, 1, 4, 6)
	println(g.tedges)
	s = serialize(g)
	println(s)
	h = unserialize(s, g.n)
	println(h.tedges)
end

function write_cliques(n::Int64)
	open("CLIQUES-" * string(n), "w") do io
		nb_clique = 0
		for g in TGraph(n)
			if isclique(g)
				nb_clique += 1
				write(io, serialize(g) * "\n")
			end
		end
		println(nb_clique)
	end;
end

function read_file(filename::String, n::Int)
	open(filename) do io
		nb_clique = 0
		for line in eachline(io)
			nb_clique += 1
			g = unserialize(line, n)
		end
		println(nb_clique)
	end;
end
