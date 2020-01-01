
# This file contains function definitions to create and manipulate
# collaboration networks

## set up packages

using SparseArrays

## The functions

"""
    histogram(s)

Count how many times different items occur in `s`.


This function is from https://benlauwens.github.io/ThinkJulia.jl/latest/book.html#dictionary_collection_counters
"""
function histogram(s)
	d = Dict()
	for c in s
		if c ∉ keys(d)
			d[c] = 1
		else
			d[c] += 1
		end
	end
	return d
end

"""
    linbins(k::Array{Int,1})

Convert the dictionary produced by `histogram(s)` into an sparse array
in increasing order.
"""
function linbins(k::Array{Int,1})
	Nk = histogram(k)
	N = length(k)
	max_k, b = maximum(Nk)
	pk = spzeros(max_k)
	for i in sort(collect(keys(Nk)))
		i <= 0.0 && continue
		pk[i] = Nk[i]/N
	end
	return pk
end

"""
    loglogbins(k)

Bin the values in `k` in a logarithmic scale.

After [Barabási](http://networksciencebook.com/chapter/4#advanced-b).
"""
function loglogbins(k)
	Nk = histogram(k)
	b_max = ceil(log2(maximum(k)))
	bin = 2 .^ collect(0:b_max)
	l_bin = length(bin)
	binval = zeros(l_bin-1)
	binsval = zeros(l_bin-1)
	for i in keys(Nk)
		for b in 1:(l_bin-1)
			if i < bin[b+1]
				binval[b] += Nk[i]
				binsval[b] += Nk[i] * i
				break
			end
		end
	end
	binsize = bin[2:end] .- bin[1:(l_bin-1)]
	pk = (binval ./ binsize) ./ length(k)
	kn = binsval ./ binval
	return Dict("pk" => pk, "kn" => kn)
end

"""
    collaborationmatrix(pubmat)

Create a collaboration matrix (showing which author publishes which
other authors) from a publication matrix (showing which paper
(co)authored by which authors).
"""
function collaborationmatrix(pubmat::SparseMatrixCSC{Float64,Int64})
	A = size(pubmat, 2)
	colnet = spzeros(A,A)
	for i in 1:(A-1)
		for j in (i+1):A
			icapj = sum(pubmat[:, i] .* pubmat[:, j])
			icapj == 0.0 && continue
			icupj = sum(pubmat[:,i]) + sum(pubmat[:,j]) - icapj
			wij = icapj/icupj
			colnet[i,j] = wij
			#colnet[j,i] = wij
		end
	end
	return colnet
end


"""
    ind2id(ind)

Produce an id from an index.
"""
function ind2id(ind)
	string(ind, base=36, pad=5)
end

"""
    id2ind(id)
	
Recover the index from an id.
"""
function id2ind(id)
	parse(Int, id, base=36)
end


"""
    collaborationgraph(colmat)

Create a collaboration graph from the collaboration matrix.

The nodes store ids created from the indexes of the nodes in the
collaboration matrix.
"""
function collaborationgraph(colmat)
	g=MetaGraph()
	for v in 1:size(colmat, 1)
		hh = ind2id(v)
		add_vertex!(g, :id, hh)
	end
	src, dst, vals = findnz(colmat)
	for i in 1:length(src)
		add_edge!(g, src[i], dst[i], :weight, vals[i])
	end
	return g
end

"""
    Weights(colnet)

Return the weights of the edges in `colnet`, a weighted MetaGraph.
"""
function Weights(colnet)
	map(ed -> LightGraphs.weights(colnet)[src(ed), dst(ed)],
			edges(colnet))
end

"""
    subnet(colnet, cutoff)

Return a subnet of `colnet`.

The subnet is formed by those vertices which is connected by edges with
weights greater than `cutoff`.
"""
function subnet(colnet, cutoff)
	heavyW = filter_edges(colnet, (g, e) ->
												LightGraphs.weights(g)[src(e), dst(e)] > cutoff)
	return colnet[heavyW]
end

"""
    cartels(colnet)

Return an array of the sizes of the connected components of `colnet`.
"""
function cartels(colnet)
	ca = connected_components(colnet)
	map(length, ca)
end

"""
    expectedpapers(A, p, gamma, k_sat)

Generate an array of expected number of papers an author produces
according to a saturated power law distribution. Computation logic
follows [Barabási's Network Science book (BOX
4.7)](http://networksciencebook.com/chapter/4#generating-networks).

# Parameters

- `A`: number of authors
- `p`: number of possible papers
- `gamma`: the exponent of the power law distribution
- `k_sat`: the saturation point
	- the power law distribution is flattened for k < k_sat
"""
function expectedpapers(A, p, gamma, k_sat)
	k = collect(1:p)
	pk = float.(k) .^(-gamma)
	pk[1:k_sat] .= pk[k_sat+1]
	pk = pk ./ sum(pk)
	Drev = cumsum(reverse(pk))
	N = rand(A)
	k_rand = zeros(A)
	for i in 1:A
		k_rand[i] = sum(Drev .> N[i])
	end
	return k_rand
end

"""
    saturatedexpectedpapers(A, p, gamma, k_sat, k_cut)

Generate an array of expected number of papers an author produces
according to a saturated power law distribution. Computation logic
follows [Barabási's Network Science book (BOX
4.7)](http://networksciencebook.com/chapter/4#generating-networks).

# Parameters

- `A`: number of authors
- `p`: number of possible papers
- `gamma`: the exponent of the power law distribution
- `k_sat`: the saturation point
	- the power law distribution is flattened for k < k_sat
- `k_cut`: the distribution drops when k > k_cut
"""
function saturatedexpectedpapers(A, p, gamma, k_sat, k_cut, alpha=1.0)
	k = float.(collect(1:p))
	px = (alpha .* (k .+ k_sat) .^ -gamma) .* exp.(k ./(-k_cut))
	#px = @. (alpha * (k + k_sat) ^ -gamma) * exp(k /(-k_cut))
	#px[1:k_sat] .= px[k_sat+1]
	px = px ./ sum(px)
	Drev = cumsum(reverse(px))
	N = rand(A)
	k_rand = zeros(A)
	for i in 1:A
		k_rand[i] = sum(Drev .> N[i])
	end
	return k_rand
end


"""
    generate_publicationmatrix(k_rand, scaling=1.0)

Generate a random publication network.

Publication networks are bipartitite networks where authors are one
kind of nodes, while papers are the other kind. An author is connected
to a paper if the author (co)writes that paper.

# Arguments

- `k_rand`: is the number of papers each author produces.
- `scaling`: scales the total number of papers the authors produce.
"""
function generate_publicationmatrix(k_rand, scaling=1.0)
	A = length(k_rand)
	p = Int(ceil(scaling * sum(k_rand)))
	if p < maximum(k_rand)
		p = Int(maximum(k_rand))
	end
	pubnet = spzeros(p, A)
	for i in 1:A
		papers = sample(1:p, Int(k_rand[i]), replace=false)
		for j in papers
			pubnet[j, i] = 1
		end
	end
	return pubnet
end

"""
    rnd_collaborationnetwork(A, gamma, k_sat, k_cut, scaling, p)

Create a random collaboration network.
"""
function rnd_collaborationnetwork(A; gamma=2.5, k_sat=10, k_cut=450,
																	scaling=1.0, p=100000)
	k_rand = saturatedexpectedpapers(A, p, gamma, k_sat, k_cut)
	pubmat = generate_publicationmatrix(k_rand, scaling)
	colmat = collaborationmatrix(pubmat)
	colnet = collaborationgraph(colmat)
	return colnet
end

"""
    describecartels(colnet, cutoff)

Gives basic information about the cartels in `colnet` deliniated by
strong links defined by `cutoff`.

Returns a Dict with the following fields:
- `no_nodes`: the number of nodes in the collaboration network,
- `cartel_nodes`,: the number of nodes concerned in cartels,
- `no_edges`: number of edges in `colnet`,
- `no_strongedges`: the number of strong edges (i.e. edges with weight >
`cutoff`),
- `cartel_sizes`: frequency histogram of cartel sizes.
"""
function describecartels(colnet, cutoff=0.4)
	W = Weights(colnet)
	sub_cn = subnet(colnet, cutoff)
	carts = cartels(sub_cn)
	histcarts = histogram(carts)
	no_v = nv(colnet)
	Dict("no_nodes" => no_v, "cartel_nodes" => nv(sub_cn),
			 "no_edges" => length(W), "no_strongedges" => sum(W .> cutoff),
			 "cartel_sizes" => histcarts)
end

"""
    sim_collaborationnetwork(A, cutoff, gamma, k_sat, scaling)

Simulate a random collaboration network and return descriptors of
cartels.

# Parameters

- `A`: number of authors
- `cutoff`: deliniate strong links between authors
- `gamma`: the exponent of the power law distribution
- `k_sat`: the saturation point
	- the power law distribution is flattened for k < k_sat
- `scaling`: scales the total number of papers the authors produce.
"""
function sim_collaborationnetwork(A, cutoff=0.4, gamma=2.5, k_sat=5,
																		scaling=1.0)
	cn = rnd_collaborationnetwork(A, gamma, k_sat, scaling, 100000)
	return describecartels(cn, cutoff)
end

"""
    run_sims(A, n_iter, cutoff, gamma, k_sat, scaling)

Run a series of simulation to collect data on cartels over `n_iter`
random collaboration network.
"""
function run_sims(A; n_iter=10, cutoff=0.4, gamma=2.5, k_sat=5, scaling=1.0)
	res = Dict()
	for i in 1:n_iter
		res[i] = sim_collaborationnetwork(A, cutoff, gamma, k_sat, scaling)
	end
	return res
end

"""
    plotloglog(x)

Plot the distribution of `x` on a log-log scale.
"""
function plotloglog(x)
	d = loglogbins(x)
	plot(log10.(d["kn"]), log10.(d["pk"]))
	scatter(log10.(d["kn"]), log10.(d["pk"]))
	title("Number of papers per author")
	xlabel(L"$\log_{10}(k)$")
	ylabel(L"$\log_{10}(p_k)$")
	PyPlot.grid(true)
	tight_layout()
end

"""
    coauthorsnumber(pubmat)

Calculate the number of authors a paper has from the publication matrix
`pubmat`.
"""
function coauthorsnumber(pubmat)
	n_papers = size(pubmat, 1)
	n_authors = zeros(n_papers)
	for i in 1:n_papers
		n_authors[i] = sum(pubmat[i,:])
	end
	return n_authors
end

"""
    graphplot(g)

Produce a pdf plot of graph `g`
"""
function graphplot(g; cutoff=0.4, backend=PDF, filename="proba",
									 layout=random_layout)
	backend == PDF && (filename *= ".pdf")
	backend == PNG && (filename *= ".png")
	W = Weights(g)
	edgewidth = ones(length(W))
	edgewidth[W .> cutoff] .= 5
	edgecolor = [colorant"lightgrey" for i in 1:length(W)]
	edgecolor[W .> cutoff] .= colorant"darkgrey"
	nodesize = degree(g)
	#nodefillc = distinguishable_colors(nv(g), colorant"blue")
	nodefillc = range(colorant"lightsalmon", stop=colorant"darksalmon",
										length=nv(g))
	#nodefillc = colorant"turquoise"
	Compose.draw(backend(filename, 50cm, 50cm),
							 gplot(g, layout=layout, edgelinewidth=edgewidth,
										 EDGELINEWIDTH=2, nodesize=nodesize,
										 nodefillc=nodefillc, edgestrokec=edgecolor))
end

"""
    addcartel!(pubmat, cartel)

Create publication cartel in the publication matrix. A cartel is a set
of authors who take part each others publications.
"""
function addcartel!(pubmat, cartel)
	cartpub = pubmat[:, cartel]
	i = sum(cartpub, dims=2) .> 0
	i = reshape(i, size(pubmat, 1))
	println(typeof(i))
	println(size(i))
	println(length(i))
	ii = collect(1:size(pubmat, 1))
	print(ii)
	ii = ii[i]
	pubmat[ii, cartel] .= 1
end

## Handle MTMT records

"""
    read_MTMT(file)

Read the MTMT records of an author from JSON file downloaded from
mtmt.hu.
"""
function read_MTMT(file)
	local r
	try
		r = JSON.parsefile(file)
	catch
		return nothing
	end
	return r["content"]
end
