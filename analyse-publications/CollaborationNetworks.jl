
# This file contains function definitions to create and manipulate
# collaboration networks

## set up packages

using SparseArrays
using Random, CSV, DataFrames
using LightGraphs, MetaGraphs

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
    labels(g)

Return the vertex labels in graph `g`
"""
function labels(g)
	map((v) -> get_prop(g, v, :id), vertices(g))
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
    cartel_sizes(colnet)

Return an array of the sizes of the connected components of `colnet`.
"""
function cartel_sizes(colnet)
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
	carts = cartel_sizes(sub_cn)
	histcarts = histogram(carts)
	no_v = nv(colnet)
	Dict("no_nodes" => no_v, "cartel_nodes" => nv(sub_cn),
			 "no_edges" => length(W), "no_strongedges" => sum(W .> cutoff),
			 "cartel_sizes" => histcarts)
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
    write_spmatrix(mat)

Write a sparse matrix `mat` to a csv file.
"""
function write_spmatrix(file, mat::SparseMatrixCSC{Float64,Int64})
	matcp = copy(mat)
	s1, s2 = size(matcp)
	matcp[1,1] <= 0.0 && (matcp[1,1] = -9999.0)
	matcp[s1,s2] <= 0.0 && (matcp[s1,s2] = -9999.0)
	rowis, colis, vals = findnz(matcp)
	df = DataFrame(Dict(:rowis => rowis, :colis => colis, :vals => vals))
	CSV.write(file, df)
	return nothing
end

"""
    read_spmatrix(file)

Read a sparse matrix from csv file `file`.
"""
function read_spmatrix(file)
	df = CSV.read(file)
	mat = sparse(df[:,:rowis], df[:,:colis], df[:,:vals])
	s1, s2 = size(mat)
	mat[1,1] == -9999.0 && (mat[1,1] = 0.0)
	mat[s1,s2] == -9999.0 && (mat[s1,s2] = 0.0)
	return dropzeros!(mat)
end

"""
    rewire!(mat, niter)

Randomly rewire in-place the bipartite graph represented by the `mat`
adjacency matrix, `niter` times.
"""
function rewire!(mat::SparseMatrixCSC{Float64, Int64}, niter=100)
	countwiring = 0
	while countwiring < niter
		rowind, colind, vals = findnz(mat)
		n_entries = length(rowind)
		ri = randperm(n_entries)
		wi = 2
		#println("outer: ", countwiring)
		while wi < n_entries
			i1 = rowind[ri[wi-1]]
			i2 = rowind[ri[wi]]
			j1 = colind[ri[wi-1]]
			j2 = colind[ri[wi]]
			if mat[i1, j2] == 0 && mat[i2, j1] == 0
				mat[i1, j1] = 0
				mat[i2, j2] = 0
				mat[i1, j2] = 1
				mat[i2, j1] = 1
				countwiring += 1
				#println("inner: ", countwiring)
			end
			wi += 2
		end
		dropzeros!(mat)
	end
	return countwiring
end

"""
    rewire(mat, niter)

Randomly rewire a copy of `mat` `niter` times.
"""
function rewire(mat::SparseMatrixCSC{Float64, Int64}, niter=100)
	matcp = copy(mat)
	rewire!(matcp, niter)
	return matcp
end

