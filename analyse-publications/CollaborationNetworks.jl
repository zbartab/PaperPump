
# This file contains function definitions to create and manipulate
# collaboration networks

## set up packages

using SparseArrays
using Random, CSV, DataFrames
using LightGraphs, MetaGraphs

## The data structures


abstract type ScienceMat end

mutable struct PubMat <: ScienceMat
	mat::SparseMatrixCSC{Float64,Int64}
	authorIDs::Dict{String,Int64}
	paperIDs::Dict{String,Int64}
end

mutable struct ColMat <: ScienceMat
	mat::SparseMatrixCSC{Float64,Int64}
	authorIDs::Dict{String,Int64}
end

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
function collaborationmatrix(pubmat::PubMat)
	A = size(pubmat.mat, 2)
	colmat = spzeros(A,A)
	for i in 1:(A-1)
		for j in (i+1):A
			icapj = sum(pubmat.mat[:, i] .* pubmat.mat[:, j])
			icapj == 0.0 && continue
			icupj = sum(pubmat.mat[:,i]) + sum(pubmat.mat[:,j]) - icapj
			wij = icapj/icupj
			colmat[i,j] = wij
		end
	end
	return ColMat(colmat, pubmat.authorIDs)
end


"""
    collaborationgraph(colmat)

Create a collaboration graph from the collaboration matrix.

The nodes store ids created from the indexes of the nodes in the
collaboration matrix.
"""
function collaborationgraph(colmat::ColMat)
	g=MetaGraph()
	for v in 1:size(colmat.mat, 1)
		idstring = findfirst((x) -> x == v, colmat.authorIDs)
		add_vertex!(g, :id, idstring)
	end
	src, dst, vals = findnz(colmat.mat)
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
    write_spmatrix(f, mat)

Write the sparse matrix part of ScienceMat into IOStream `f`.
"""
function write_spmatrix(f::IOStream, mat::SparseMatrixCSC{Float64, Int64})
	rowis, colis, vals = findnz(mat)
	for i in 1:length(rowis)
		println(f, rowis[i], ",", colis[i], ",", vals[i])
	end
	return nothing
end

"""
    write_IDs(f, IDs)

Write the IDs part of ScienceMat in to IOStream `f`.
"""
function write_IDs(f::IOStream, IDs::Dict{String,Int64})
	for id in keys(IDs)
		println(f, id, ",", IDs[id])
	end
	return nothing
end

"""
    write_scimat(file, mat)

Write a representation of ScienceMat `mat` to `file`.
"""
function write_scimat(file::String, mat::ScienceMat)
	smat = size(mat.mat)
	f = open(file, "w")
	if typeof(mat) == PubMat
		print(f, "#### pubmat,")
	else
		print(f, "#### colmat,")
	end
	println(f, smat[1], ",", smat[2], ",", nnz(mat.mat))
	write_spmatrix(f, mat.mat)
	println(f, "#### authorIDs,", length(mat.authorIDs))
	write_IDs(f, mat.authorIDs)
	if typeof(mat) == PubMat
		println(f, "#### paperIDs,", length(mat.paperIDs))
		write_IDs(f, mat.paperIDs)
	end
	close(f)
	return nothing
end

"""
    read_scimat(file)

Read `ScienceMat` data from file `file`.
"""
function read_scimat(file::String)
	f = open(file)
	lines = readlines(f)
	close(f)
	i = 1
	matrix_type, nr, nc, ne = split(lines[i], ",")
	nrow = parse(Int, nr)
	ncol = parse(Int, nc)
	nrec = parse(Int, ne)
	ris = Array{Int,1}(undef, nrec)
	cis = Array{Int,1}(undef, nrec)
	vis = Array{Float64,1}(undef, nrec)
	for j in (i+1):(nrec+i)
		jj= j-i
		ri, ci, vi = split(lines[j], ",")
		ris[jj] = parse(Int, ri)
		cis[jj] = parse(Int, ci)
		vis[jj] = parse(Float64, vi)
	end
	mat = sparse(ris, cis, vis, nrow, ncol)
	i += nrec+1
	id_type, ne = split(lines[i], ",")
	nrec = parse(Int, ne)
	authorIDs = Dict{String, Int64}()
	for j in (i+1):(i+nrec)
		jj = j-i
		id, ind = split(lines[j], ",")
		ni = parse(Int, ind)
		authorIDs[id] = ni
	end
	if occursin(r"colmat$", matrix_type)
		return ColMat(mat, authorIDs)
	else
		i += nrec+1
		id_type, ne = split(lines[i], ",")
		nrec = parse(Int, ne)
		paperIDs = Dict{String, Int64}()
		for j in (i+1):(i+nrec)
			jj = j-i
			id, ind = split(lines[j], ",")
			ni = parse(Int, ind)
			paperIDs[id] = ni
		end
		return PubMat(mat, authorIDs, paperIDs)
	end
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
function rewire(pubmat::PubMat, niter=100)
	matcp = deepcopy(pubmat)
	rewire!(matcp.mat, niter)
	return matcp
end

#=
"""
    samplepapers(npapers, pm)

Sample `npapers` papers randomly from publication matrix `pm`. It
removes authors with zero papers from the resulting matrix.
"""
function samplepapers(pm::PubMat, npapers::Int)
	spm =  size(pm)
	npapers > spm[1] && (npapers = spm[1])
	npapers = sample(1:spm[1], npapers, replace=false)
	np = papernumbers(pm[npapers,:])
	return pm[npapers, np .> 0]
end

"""
    sampleauthors(nauthors, pm)

Sample `nauthors` authors randomly from publication matrix `pm`. It
removes papers with zero authors from the resulting matrix.
"""
function sampleauthors(nauthors::Int, pm::SparseMatrixCSC{Float64,Int64})
	spm =  size(pm)
	nauthors > spm[2] && (nauthors = spm[2])
	nauthors = sample(1:spm[2], nauthors, replace=false)
	na = authornumbers(pm[:, nauthors])
	return pm[na .> 0, nauthors]
end
=#

"""
    papernumbers(pm)

Returns an array with the number of papers each author in publication
matrix `pm` authored.
"""
function papernumbers(pm::PubMat)
	return papernumbers(pm.mat)
end
function papernumbers(mat::SparseMatrixCSC{Float64,Int64})
	no_papers = sum(mat, dims=1)
	no_papers = reshape(no_papers, size(no_papers, 2))
	return no_papers
end

"""
    authornumbers(pm)

Returns an array with the number of authors each paper in publication
matrix `pm` has.
"""
function authornumbers(pm::PubMat)
	return authornumbers(pm.mat)
end
function authornumbers(mat::SparseMatrixCSC{Float64,Int64})
	no_authors = sum(mat, dims=2)
	no_authors = reshape(no_authors, size(no_authors, 1))
	return no_authors
end

"""
"""
function updateIDs(IDs::Dict{String,Int64}, index::BitArray{1})
	@assert length(IDs) == length(index) "IDs length does not match index length!"
	id_s = Array{String,1}(undef, length(IDs))
	for k in keys(IDs)
		id_s[IDs[k]] = k
	end
	id_s = id_s[index]
	newIDs = Dict{String,Int64}()
	for i in 1:length(id_s)
		newIDs[id_s[i]] = i
	end
	return newIDs
end

"""
    selectauthors(pm, npapers)

Select only those authors from publication matrix `pm` who have more
papers than `npapers`. Papers with zero authors also removed from the
resulting matrix.
"""
function selectauthors(pm::PubMat, npapers=-Inf)
	np = papernumbers(pm)
	i = np .> npapers
	pmred = pm.mat[:, i]
	na = authornumbers(pmred)
	j = na .> 0
	pmred = pmred[j, :]
	return PubMat(pmred,
								updateIDs(pm.authorIDs, i),
								updateIDs(pm.paperIDs, j))
end

"""
    selectpapers(pm, nauthors)

Select only those papers from publication matrix `pm` which have more
authors than `nauthors`. Papers with zero authors also removed from the
resulting matrix.
"""
function selectpapers(pm::PubMat, nauthors=-Inf)
	np = authornumbers(pm)
	j = np .> nauthors
	pmred = pm.mat[j, :]
	na = papernumbers(pmred)
	i = na .> 0
	pmred = pmred[:, i]
	return PubMat(pmred,
								updateIDs(pm.authorIDs, i),
								updateIDs(pm.paperIDs, j))
end

