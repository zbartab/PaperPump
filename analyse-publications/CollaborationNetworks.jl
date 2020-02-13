
# This file contains function definitions to create and manipulate
# collaboration networks

## set up packages

import Base.size # to extend for returning the size of the matrix component

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
    collaborationmatrix(pubmat)

Create a collaboration matrix (showing which author publishes which
other authors) from a publication matrix (showing which paper
(co)authored by which authors).
"""
function collaborationmatrix(pubmat::PubMat)
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
    Weights(colnet)

Return the weights of the edges in `colnet`, a weighted MetaGraph.
"""
function Weights(colnet)
	map(ed -> LightGraphs.weights(colnet)[src(ed), dst(ed)],
			edges(colnet))
end

"""
		strength(colmat)

Calculate the strength of each node in `colmat`. Strength of node i is
the sum of weights of edges connected to node i.
"""
function strength(colmat::ColMat)
	m = colmat.mat
	n = size(m, 1)
	s = Array{Float64,1}(undef, n)
	for i in 1:n
		s[i] = sum(m[i,:]) + sum(m[:,i])
	end
	return s
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

"""
    getauthorindex(ids, coma)

Return the indices in collaboration matrix `coma` of authors listed in
`ids`.
"""
function getauthorindex(ids::Array{String,1}, coma::ScienceMat)
	indices = Int[]
	for id in ids
		push!(indices, coma.authorIDs[id])
	end
	return indices
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
		while wi < n_entries && countwiring < niter
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

"""
"""
function rewire_community(pm::PubMat, members::Array{Int,1},
													piter::Float64=0.1)
	nonmembers = setdiff(collect(1:size(pm.mat, 2)), members)
	cpm = deepcopy(pm)
	m_nm = cpm.mat[:,nonmembers]
	m_cm = cpm.mat[:,members]
	r, c, v = findnz(m_cm)
	niter = Int(round(sum(v .> 0.0) * piter))
	niter == 0 && return pm
	rewire!(m_cm, niter)
	m = hcat(m_nm, m_cm)
	a = Array{String, 1}(undef, length(pm.authorIDs))
	for i in keys(pm.authorIDs)
		a[pm.authorIDs[i]] = i
	end
	a = a[vcat(nonmembers, members)]
	aIDs = Dict{String, Int64}()
	for i in 1:length(a)
		aIDs[a[i]] = i
	end
	return PubMat(m, aIDs, pm.paperIDs) 
end

"""
"""
function rewire_communities(pubmat::PubMat,
														communities::Dict{Int, Array{Int,1}},
														prewire::Float64, file="comm_rewire")
	pm = deepcopy(pubmat)
	for i in keys(communities)
		members = communities[i]
		#println("i: ", i, ", community size: ", length(members))
		length(members) < 2 && continue
		pm = rewire_community(pm, members, prewire)
	end
	filename = string(file, "-", prewire)
	write_scimat(string(filename, "-pubmat.mat"), pm)
	res = Dict()
	pm = selectauthors(pm, 3)
	res[:puma] = pm
	rcm = collaborationmatrix(pm)
	res[:coma] = rcm
	fn = string(filename, "-colmat.mat")
	write_scimat(fn, rcm)
	rcg = collaborationgraph(rcm)
	res[:coga] = rcg
	run(`./comm_detection.sh $(fn)`)
	f = open(string(filename, "-colmat.Q"))
	Q = readlines(f)
	close(f)
	Q = parse(Float64, Q[1])
	return Q, res
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

"""
    read_louvain_tree(file)

Read the results of louvain community detection from `file`.
"""
function read_louvain_tree(file::String)
	f = open(file)
	lines = readlines(f)
	close(f)
	top = Dict()
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
    combine_sparsematrices(m1, m2)

Combine the two matrices given as argument.
"""
function combine_sparsematrices(m1::SparseMatrixCSC{Float64,Int64},
																m2::SparseMatrixCSC{Float64,Int64})
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
    degreedegreecor(pm)

Calculate the join degree distribution of a bipartite network.
"""
function degreedegreecor(pm::PubMat)
	da = Int.(papernumbers(pm))
	dp = Int.(authornumbers(pm))
	p, A, val = findnz(pm.mat)
	ddc = spzeros(Int, Int(maximum(dp)), Int(maximum(da)))
	for i in 1:length(val)
		ddc[dp[p[i]], da[A[i]]] += 1
	end
	return ddc
end

function fillnodes(distr::Array{Int,1})
	n = sum(distr)
	degrees = Array{Int, 1}(undef, n)
	stubs = Array{Int, 1}(undef, n)
	i = 1
	for j in 1:length(distr)
		for d in 1:distr[j]
			degrees[i] = j
			stubs[i] = j
			i += 1
		end
	end
	i = randperm(n)
	return degrees[i], stubs[i]
end

"""
    generaterndbipartite(bjd)

Generate a random bipartite graph on the basis of the bivariable joint
degree distribution. It follows [this
paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5972839/)
"""
function generaterndbipartite(bjd::SparseMatrixCSC{Int64, Int64})
	pd = Int.(sum(bjd, dims=2) ./ (1:size(bjd, 1)))
	pd = reshape(pd, size(pd, 1))
	Ad = Int.(sum(bjd, dims=1)' ./ (1:size(bjd, 2)))
	Ad = reshape(Ad, size(Ad, 1))
	np = sum(pd)
	nA = sum(Ad)
	ps_d, ps_s = fillnodes(pd)
	As_d, As_s = fillnodes(Ad)
	pin, Ain, nedges = findnz(bjd)
	m = spzeros(length(ps_d), length(As_d))
	#println("O1: ", length(nedges))
	counter = 0
	while length(nedges) > 0 && counter < 100
		#println("I1: ", length(nedges), " p: ", pin, " A: ", Ain, " e: ",
					 #nedges)
		i = randperm(length(pin))
		pin = pin[i]
		Ain = Ain[i]
		nedges = nedges[i]
		for i in 1:length(pin)
			p = findfirst((k) -> (ps_d[k] == pin[i] && ps_s[k] > 0), 1:np)
			A = findfirst((k) -> (As_d[k] == Ain[i] && As_s[k] > 0), 1:nA)
			ps_s[p] -= 1
			As_s[A] -= 1
			nedges[i] -= 1
			m[p,A] += 1.0
		end
		i = nedges .> 0
		pin = pin[i]
		Ain = Ain[i]
		nedges = nedges[i]
		#println("I2: ", length(nedges), " p: ", pin, " A: ", Ain, " e: ",
					 #nedges)
		counter += 1
	end
	#println("O2: ", length(nedges))
	#return ps_s, ps_d, As_s, As_d, m
	return m
end

"""
    simplifybipartite(m)

Simplify a randomly generated bipartite graph represented by `m`.
Simplification means to resolve multiply edges between nodes. The
algorithm follows [this
paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5972839/)
"""
function simplifybipartite!(m::SparseMatrixCSC{Float64, Int64})
	p, A = size(m)
	ps = collect(1:p)
	As = collect(1:A)
	n_medges = sum(values(m) .> 1)
	n_medges == 0 && return nothing
	udegrees = Int.(sum(m, dims=1))
	udegrees = reshape(udegrees, size(udegrees, 2))
	ldegrees = Int.(sum(m, dims=2))
	ldegrees = reshape(ldegrees, size(ldegrees, 1))
	rind, cind, vals = findnz(m)
	i_medges = findall((y) -> y > 1, vals)
	shuffle!(i_medges)
	for i in i_medges
		u = cind[i]
		l = rind[i]
		ld = ldegrees[l]
		lds = findall((k) -> ld == k, ldegrees)
		lds = lds[lds .!= l]
		breakout = false
		if length(lds) > 0
			shuffle!(lds)
			for v in lds
				ws = m[v, :]
				ws = As[ws .> 0.0]
				shuffle!(ws)
				for w in ws
					if m[l,w] == 0.0 && m[v,u] == 0.0
						m[l,u] -= 1.0
						m[v,w] -= 1.0
						m[l,w] = 1.0
						m[v,u] = 1.0
						breakout = true
						break
					end
				end
				breakout && break
			end
		end
		breakout && continue
		ud = udegrees[u]
		uds = findall((k) -> ud == k, udegrees)
		uds = uds[uds .!= u]
		breakout = false
		if length(uds) > 0
			shuffle!(uds)
			for w in uds
				vs = m[:, w]
				vs = ps[vs .> 0.0]
				shuffle!(vs)
				for v in vs
					if m[l,w] == 0.0 && m[v,u] == 0.0
						m[l,u] -= 1.0
						m[v,w] -= 1.0
						m[l,w] = 1.0
						m[v,u] = 1.0
						breakout = true
						break
					end
				end
				breakout && break
			end
		end
		breakout && continue
		ws = m[:,u]
		ws = ps[ws .<= 0.0]
		shuffle!(ws)
		wps = m[l,:]
		wps = As[wps .<= 0.0]
		shuffle!(wps)
		for w in ws, wp in wps
			ups = m[w,:]
			ups = As[ups .> 0.0]
			shuffle!(ups)
			vps = m[:,wp]
			vps = ps[vps .> 0.0]
			shuffle!(vps)
			for up in ups, vp in vps
				if m[vp,up] <= 0.0 
					m[l,u] -= 1.0
					m[w,up] -= 1.0
					m[vp,wp] -= 1.0
					m[w,u] = 1.0
					m[l,wp] = 1.0
					m[vp,up] = 1.0
					breakout = true
					break
				end
			end
			breakout && break
		end
		breakout && continue
		println("simplification is unsuccessful")
	end
	return nothing
end

function checkdegreedegree(bjd::SparseMatrixCSC{Int64, Int64})
	pd = Int.(sum(bjd, dims=2) ./ (1:size(bjd, 1)))
	pd = reshape(pd, size(pd, 1))
	Ad = Int.(sum(bjd, dims=1)' ./ (1:size(bjd, 2)))
	Ad = reshape(Ad, size(Ad, 1))
	return pd, Ad
end
