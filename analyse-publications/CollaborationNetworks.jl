
# This file contains function definitions to create and manipulate
# collaboration networks

## set up packages

using SparseArrays
using JSON
using Colors

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
	tocutoff = W .> cutoff
	edgewidth = ones(length(W))
	edgewidth[tocutoff] .= 5
	edgecolor = [colorant"lightgrey" for i in 1:length(W)]
	edgecolor[tocutoff] .= colorant"black"
	nodesize = degree(g)
	#nodefillc = distinguishable_colors(nv(g), colorant"blue")
	#nodefillc = range(colorant"lightsalmon", stop=colorant"darksalmon",
										#length=nv(g))
	nodefillc = colorant"turquoise"
	if sum(tocutoff) == 0
		maxlinewidth = 0.2
	else
		maxlinewidth = 2
	end
	Compose.draw(backend(filename, 50cm, 50cm),
							 gplot(g, layout=layout, edgelinewidth=edgewidth,
										 EDGELINEWIDTH=maxlinewidth, nodesize=nodesize,
										 nodefillc=nodefillc, edgestrokec=edgecolor))
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
