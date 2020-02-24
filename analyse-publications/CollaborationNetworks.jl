
# This file contains function definitions to create and manipulate
# collaboration networks

## set up packages

import Base.size # to extend for returning the size of the matrix component

using SparseArrays
using Random, CSV, DataFrames
using LightGraphs, MetaGraphs

#include("StatCollaborationNetworks.jl")
#include("PlotCollaborationNetworks.jl")
#include("simulations/RandomPublicationNetworks.jl")

## The data structures

const HierarchyType = Dict{Symbol, Dict{Int64, Array{Int64, 1}}}
const HierarchyIDType = Dict{Symbol, Dict{Int64, Array{String, 1}}}

abstract type ScienceMat end

mutable struct PubMat <: ScienceMat
	mat::SparseMatrixCSC{Int64,Int64}
	authorIDs::Dict{String,Int64}
	paperIDs::Dict{String,Int64}
end

mutable struct ColMat <: ScienceMat
	mat::SparseMatrixCSC{Float64,Int64}
	authorIDs::Dict{String,Int64}
end

mutable struct PubNet
	puma::PubMat
	coma::ColMat
	coga::MetaGraph{Int64, Float64}
	nauthors::Array{Int,1}
	npapers::Array{Int,1}
	wpapers::Array{Float64,1}
	strengthes::Array{Float64,1}
	weights::Array{Float64,1}
	degrees::Array{Int,1}
	clustcoefs::Array{Float64,1}
	Q::Float64
	hierarchy::HierarchyIDType
	cutoff::Float64
	cartels::Array{Array{String,1},1}
	function PubNet(puma::PubMat, coma::ColMat, cutoff::Float64=0.4)
		coga = collaborationgraph(coma)
		nauthors = authornumbers(puma)
		npapers = papernumbers(puma)
		wpapers = weightedpapers(puma)
		strengthes = strength(coma)
		weights = Weights(coga)
		degrees = degree(coga)
		clustcoefs = local_clustering_coefficient(coga)
		file2 = "tmp.mat"
		write_scimat(file2, coma)
		treefile = replace(file2, r"\.mat$" => ".tree")
		Qfile = replace(file2, r"\.mat$" => ".Q")
		if isfile(treefile) && mtime(treefile) > mtime(file2)
			H = read_louvain_tree(treefile)
			Q = read_louvain_Q(Qfile)
		else
			run(`./comm_detection.sh $(file2)`)
			H = read_louvain_tree(treefile)
			Q = read_louvain_Q(Qfile)
		end
		H = converthierarchy(coga, H)
		carts = cartels(coga, cutoff)
		return new(puma, coma, coga, nauthors, npapers, wpapers, strengthes,
							 weights, degrees, clustcoefs, Q, H, cutoff,
							carts)
	end
end

## The functions

"""
    PubNet(file, cutoff)

Outer constructor for `PubNet` reading data from files.
"""
function PubNet(file::String, cutoff::Float64=0.4)
	file2 = replace(file, "pubmat" => "colmat")
	puma = read_scimat(file)
	coma = read_scimat(file2)
	return PubNet(puma, coma, cutoff)
end

"""
    collaborationmatrix(pubmat)

Create a collaboration matrix (showing which author publishes which
other authors) from a publication matrix (showing which paper
(co)authored by which authors).
"""
function collaborationmatrix(pubmat::PubMat)
	m = pubmat.mat
	A = size(m, 2)
	colmat = spzeros(A,A)
	for i in 1:(A-1)
		mi = m[:,i]
		for j in (i+1):A
			mj = m[:,j]
			icapj = sum(mi .* mj)
			icapj == 0.0 && continue
			icupj = sum(mi) + sum(mj) - icapj
			wij = icapj/icupj
			colmat[i,j] = wij
		end
	end
	return ColMat(colmat, pubmat.authorIDs)
end
function collaborationmatrix2(pubmat::PubMat)
	A = size(pubmat.mat, 2)
	colmat = spzeros(A,A)
	for i in 1:(A-1)
		mi = pubmat.mat[:,i]
		for j in (i+1):A
			mj = pubmat.mat[:,j]
			icapj = sum(mi .* mj)
			icapj == 0.0 && continue
			icupj = sum(mi) + sum(mj) - icapj
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
function subnet(colnet::MetaGraph{Int64, Float64}, cutoff::Float64)
	heavyW = filter_edges(colnet, (g, e) ->
												LightGraphs.weights(g)[src(e), dst(e)] > cutoff)
	return colnet[heavyW]
end

"""
    getauthorIDs(colnet)

Return the author IDs in the collaboration graph `colnet`.
"""
function getauthorIDs(colnet::MetaGraph{Int64, Float64})
	ids = String[]
	for v in vertices(colnet)
		id = get_prop(colnet, v, :id)
		push!(ids, id)
	end
	return ids
end
function getauthorIDs(colnet::MetaGraph{Int64, Float64},
											group::Array{Int,1})
	ids = getauthorIDs(colnet)
	return ids[group]
end

"""
    getauthorindex(coma, ids)

Return the indices in collaboration matrix `coma` of authors listed in
`ids`.
"""
function getauthorindex(coma::ScienceMat, ids::Array{String,1})
	indices = Int[]
	for id in ids
		push!(indices, coma.authorIDs[id])
	end
	return indices
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
	rowis, colis, vals = findnz(mat.mat)
	for i in 1:length(rowis)
		println(f, rowis[i], ",", colis[i], ",", vals[i])
	end
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
		mat = sparse(ris, cis, vis, nrow, ncol)
		return ColMat(mat, authorIDs)
	else
		mat = sparse(ris, cis, Int.(vis), nrow, ncol)
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
function updateIDs(IDs::Dict{String,Int64}, index::Array{Int64,1})
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
    selectauthors(pm, authors)

Return the part of publication matrix `pm` which belongs to the authors
in `authors`.
"""
function selectauthors(pm::PubMat, authors::Array{Int,1})
	pmred = pm.mat[:, authors]
	na = authornumbers(pmred)
	j = na .> 0
	pmred = pmred[j, :]
	return PubMat(pmred,
								updateIDs(pm.authorIDs, authors),
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

"""
    read_louvain_tree(file)

Read the results of louvain community detection from `file`.
"""
function read_louvain_tree(file::String)
	f = open(file)
	lines = readlines(f)
	close(f)
	#top = Dict{Symbol, Dict{Int64, Array{Int64,1}}}()
	top = HierarchyType()
	level = 0
	levsymb = :a
	prevlevsymb = :a
	for l in lines
		n, c = split(l, " ")
		node = parse(Int, n)
		community = parse(Int, c)
		if node == 0 && community == 0
			level += 1
			levsymb = Symbol("lev$(level)")
			prevlevsymb = Symbol("lev$(level-1)")
			top[levsymb] = Dict{Int, Array{Int,1}}()
		else
			if !haskey(top[levsymb], community)
				top[levsymb][community] = Int64[]
			end
			if level == 1
				push!(top[levsymb][community], node)
			else 
				a = top[levsymb][community]
				top[levsymb][community] = vcat(a, top[prevlevsymb][node])
			end
		end
	end
	return top
end

"""
    converthierarchy(coga, hierarchy)

Convert the community hierarchy found for `coga` from Int index based
respresentation to the String based ID representation.
"""
function converthierarchy(coga::MetaGraph{Int64, Float64},
													hierarchy::HierarchyType)
	auIDs = getauthorIDs(coga)
	H = HierarchyIDType()
	for l in keys(hierarchy)
		H[l] = Dict{Int64, Array{String, 1}}()
		for i in keys(hierarchy[l])
			H[l][i] = auIDs[hierarchy[l][i]]
		end
	end
	return H
end

"""
    combine_sparsematrices(m1, m2)

Combine the two matrices given as argument.
"""
function combine_sparsematrices(m1::SparseMatrixCSC{Int64,Int64},
																m2::SparseMatrixCSC{Int64,Int64})
	r1, c1 = size(m1)
	r2, c2 = size(m2)
	ri1, ci1, vi1 = findnz(m1)
	ri2, ci2, vi2 = findnz(m2)
	ri2 .+= r1
	ci2 .+= c1
	r = vcat(ri1, ri2)
	c = vcat(ci1, ci2)
	v = vcat(vi1, vi2)
	return sparse(r, c, v, r1+r2, c1+c2)
end

"""
    size(m)

Return the size of the matrix component of a `ScienceMat` object.
"""
function size(m::ScienceMat)
	return size(m.mat)
end

"""
    read_louvain_Q(file)

Return network modularity Q read from `file`. If `file` cannot be read
return `false` and print an error meassesge.
"""
function read_louvain_Q(file)
	if isfile(file)
		f = open(file)
		Q = readlines(f)
		close(f)
		Q = parse(Float64, Q[1])
		return Q
	else
		println("ERROR: $(file) cannot be read!")
		return false
	end
end

"""
    processpubnet(file)

Load the publication network specified by `file` and calculate several
network measures for it.
"""
function processpubnet(file)
	file2 = replace(file, "pubmat" => "colmat")
	puma = read_scimat(file)
	coma = read_scimat(file2)
	coga = collaborationgraph(coma)
	pubnet = Dict()
	pubnet[:puma] = puma
	pubnet[:coma] = coma
	pubnet[:coga] = coga
	pubnet[:nauthors] = authornumbers(pubnet[:puma])
	pubnet[:npapers] = papernumbers(pubnet[:puma])
	pubnet[:wpapers] = weightedpapers(pubnet[:puma])
	pubnet[:strength] = strength(pubnet[:coma])
	pubnet[:weights] = Weights(pubnet[:coga])
	pubnet[:degree] = degree(pubnet[:coga])
	pubnet[:clustcoef] = local_clustering_coefficient(pubnet[:coga])
	treefile = replace(file2, r"\.mat$" => ".tree")
	Qfile = replace(file2, r"\.mat$" => ".Q")
	if isfile(treefile) && mtime(treefile) > mtime(file2)
		H = read_louvain_tree(treefile)
		Q = read_louvain_Q(Qfile)
	else 
		run(`./comm_detection.sh $(file2)`)
		H = read_louvain_tree(treefile)
		Q = read_louvain_Q(Qfile)
	end
	pubnet[:Q] = Q
	pubnet[:hierarchy] = H
	return pubnet
end

