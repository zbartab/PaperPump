
# Functions used to simulate publication networks and cartels

using StatsBase, Distributions

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

Default values are those derived for the MTMT dataset.
"""
function saturatedexpectedpapers(A, p=100000, gamma=2.5, k_sat=10,
																 k_cut=450, alpha=1.0)
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
    generate_publicationmatrix(k_rand, maxpapers=100_000)

Generate a random publication network.

Publication networks are bipartitite networks where authors are one
kind of nodes, while papers are the other kind. An author is connected
to a paper if the author (co)writes that paper.

# Arguments

- `k_rand`: is the number of papers each author produces.
- `maxpapers`: the maximal number of papers the authors can produce.
"""
function generate_publicationmatrix(k_rand, maxpapers=100_000)
	A = length(k_rand)
	p = maxpapers
	if p < maximum(k_rand)
		p = Int(maximum(k_rand))
	end
	pubnet = spzeros(Int, p, A)
	for i in 1:A
		papers = sample(1:p, Int(k_rand[i]), replace=false)
		for j in papers
			pubnet[j, i] = 1
		end
	end
	an = authornumbers(pubnet)
	pubnet = pubnet[an .> 0,:]
	p, A = size(pubnet)
	auIDs = Dict{String,Int64}()
	for i in 1:A
		auIDs["A$(i)"] = i
	end
	paIDs = Dict{String,Int64}()
	for i in 1:p
		paIDs["p$(i)"] = i
	end
	return PubMat(pubnet, auIDs, paIDs)
end
function generate_publicationmatrix(mat::SparseMatrixCSC{Int64,Int64})
	p, A = size(mat)
	auIDs = Dict{String,Int64}()
	for i in 1:A
		auIDs["A$(i)"] = i
	end
	paIDs = Dict{String,Int64}()
	for i in 1:p
		paIDs["p$(i)"] = i
	end
	return PubMat(mat, auIDs, paIDs)
end

"""
    grow_pubnet(k, p_new, alfa)

Grow a full publication network.
"""
function grow_pubnet(k::Array{Int64,1}, p_new::Float64, alfa::Float64)
	pn = Dict()
	pm = grow_publicationmatrix(k, p_new, alfa)
	pn[:puma] = generate_publicationmatrix(pm)
	pn[:puma] = selectauthors(pn[:puma], 3)
	pn[:coma] = collaborationmatrix(pn[:puma])
	pn[:coga] = collaborationgraph(pn[:coma])
	return pn
end

"""
    grow_publicationmatrix(k, p_new, alfa)

Grow a publication matrix similarly to the Barabási-Albert model. `k`
gives the number of papers each author should have, `p_new` is the
probability to start a new paper and `alfa` is the exponent to shape the
preferential attachement.

An author start a new paper with probability `p_new`. With probability
`1-p_new` it joins a paper already started by someone else. The
propbability to join a given paper is proportional to the number of
authors on that particular paper. This later proportional propbability
is shaped by `alfa`.
"""
function grow_publicationmatrix(k::Array{Int64,1}, p_new::Float64,
																alfa::Float64)
	maxpapers = sum(k)
	A = length(k)
	pm = spzeros(Int, maxpapers, A+3)
	papers = zeros(Int, maxpapers)
	pk = zeros(Float64, maxpapers)
	# set the initial authors' node
	pm[1,1] = 1.0
	pm[2,2] = 1.0
	pm[3,3] = 1.0
	npapers = 3 # counter of papers
	papers[1:npapers] .= 1
	pk[1:npapers] .= update_pk(papers[1:npapers], alfa)
	# add new authors
	for j in 4:A
		for p in 1:k[j]
			# step over the number of papers author j has
			if rand() < p_new
				# start new paper
				npapers += 1
				pm[npapers, j] = 1
				papers[npapers] = 1
			else
				# join existing paper
				i = 0
				ii = 0
				while true
					ii += 1
					i = papertojoin(pk)
					(pm[i,j] == 0 || ii > npapers) && break
				end
				if ii > npapers
					npapers += 1
					pm[npapers, j] = 1
					papers[npapers] = 1
				else
					pm[i,j] = 1
					papers[i] += 1
				end
			end
			pk[1:npapers] .= update_pk(papers[1:npapers], alfa)
		end
	end
	maxpapers = maximum(rowvals(pm))
	pm = pm[1:maxpapers, 4:end]
	return pm
end

"""
    update_pk(pprs, alfa)

Update the probability distribution of choosing an existing paper to
join. `pk` contains the distribution, `papers` has the number of authors
for each paper.
"""
function update_pk(pprs::Array{Int64,1}, alfa::Float64)
	pk = pprs .^ alfa
	return pk ./ sum(pk)
end

"""
    papertojoin(pk)

Return the index of already exisiting paper to join. The probability
distribution of choosing a paper is given in `pk`.
"""
function papertojoin(pk::Array{Float64,})
	l_pk = length(pk)
	p = rand()
	sump = 0.0
	for i in 1:l_pk
		sump += pk[i]
		sump > p && return i
	end
	#return l_pk
end

"""
    rnd_collaborationnetwork(A, gamma, k_sat, k_cut, maxpapers, p)

Create a random collaboration network.
"""
function rnd_collaborationnetwork(A; gamma=2.5, k_sat=10, k_cut=450,
																	maxpapers=100_000, p=100_000)
	k_rand = saturatedexpectedpapers(A, p, gamma, k_sat, k_cut)
	pubmat = generate_publicationmatrix(k_rand, maxpapers)
	colmat = collaborationmatrix(pubmat)
	colnet = collaborationgraph(colmat)
	return colnet
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
function sim_collaborationnetwork(A, cutoff=0.5, gamma=2.5, k_sat=5,
																		maxpapers=100_000)
	cn = rnd_collaborationnetwork(A, gamma, k_sat, maxpapers, 100000)
	return describecartels(cn, cutoff)
end

"""
    run_sims(A, n_iter, cutoff, gamma, k_sat, maxpapers)

Run a series of simulation to collect data on cartels over `n_iter`
random collaboration network.
"""
function run_sims(A; n_iter=10, cutoff=0.5, gamma=2.5, k_sat=5,
									maxpapers=100_000)
	res = Dict()
	for i in 1:n_iter
		res[i] = sim_collaborationnetwork(A, cutoff, gamma, k_sat, maxpapers)
	end
	return res
end

"""
    addcartel!(pubmat, cartel)

Create publication cartel in the publication matrix. A cartel is a set
of authors who take part each others publications.
"""
function addcartel!(pubmat::SparseMatrixCSC{Int64, Int64},
										cartel::Array{Int,1}, collprob=1.0)
	#println(typeof(pubmat))
	cartpub = pubmat[:, cartel]
	i = sum(cartpub, dims=2) .> 0
	i = reshape(i, size(pubmat, 1))
	j = collect(1:size(pubmat, 1))[i]
	for jj in j, c in cartel
		if pubmat[jj,c] <= 0 && rand() < collprob
			pubmat[jj,c] = 1
		end
	end
	return nothing
end
function addcartel!(pubmat::PubMat, cartel::Array{Int,1}, collprob=1.0)
	addcartel!(pubmat.mat, cartel, collprob)
	return nothing
end

"""
    formcartel(pm, cartel, collprob)

Form a publication cartel by rewiring links so that for each link added
between authors part of the cartel another link is removed to preserve the
degrees distribution.
"""
function formcartel!(pubmat::PubMat, cartel::Array{Int,1}, collprob=1.0)
	cartpub = pubmat.mat[:, cartel]
	i = sum(cartpub, dims=2) .> 0
	i = reshape(i, size(pubmat.mat, 1))
	j = collect(1:size(pubmat.mat, 1))[i]
	for jj in j, c in cartel
		if pubmat.mat[jj,c] == 0 && rand() < collprob
			pnzs = findall((k) -> k == 1, pubmat.mat[:,c])
			pnzs = setdiff(pnzs, jj)
			if length(pnzs) > 0
				shuffle!(pnzs)
				for pnz in pnzs
					if pubmat.mat[pnz, c] == 1
						breakout = false
						cnzs = findall((c) -> c == 1, pubmat.mat[jj,:])
						cnzs = setdiff(cnzs, cartel)
						if length(cnzs) > 0
							shuffle!(cnzs)
							for cnz in cnzs
								if pubmat.mat[pnz, cnz] == 0 && pubmat.mat[jj,cnz] == 1
									pubmat.mat[jj,c] = 1
									pubmat.mat[jj,cnz] = 0
									pubmat.mat[pnz,cnz] = 1
									pubmat.mat[pnz,c] = 0
									breakout = true
									break
								end
							end
						end
						breakout && break
					end
				end
			end
		end
	end
	return nothing
end

"""
    rndpubnet(k, pratio, commsizes, cartsizes, p_carts)

Create a random publication matrix in which the authors' paper
distribution is sampled from `k`, the number of papers is around
`pratio` times the author numbers. The generated network builds up from
communities whose size is given by `commsizes`.
"""
function rndpubnet3(k, pratio, commsizes, p_rewire=0.01,
										cartsizes=nothing, p_carts=0.2, α=5, β=1)
	local pm
	first = true
	csizes = []
	cartels = Array{Array{Int, 1}, 1}()
	for cs in commsizes
		k0 = sample(k, cs, replace=false)
		p = Int(round(cs*pratio))
		pm0 = generate_publicationmatrix(k0, p)
		push!(csizes, size(pm0))
		if first
			pm = pm0
			first = false
		else
			m = combine_sparsematrices(pm.mat, pm0.mat)
			pm = generate_publicationmatrix(m)
		end
	end
	pmrw = rewire(pm, p_rewire)
	pm = deepcopy(pmrw)
	if !isnothing(cartsizes)
		pmc = deepcopy(pmrw)
		sampledist = Beta(α,β)
		na = size(pmc, 2)
		ncartmembs = Int(round(p_carts * na))
		local scartsizes
		while true
			scartsizes = cumsum(cartsizes)
			if scartsizes[end] > ncartmembs
				break
			else
				cartsizes = vcat(cartsizes, cartsizes)
			end
		end
		ncarts = findlast(scartsizes .<= ncartmembs)
		cartmembs = sort(sample(1:na, Int(scartsizes[ncarts]), replace=false))
		fi = 1
		for i in 1:ncarts
			li = Int(fi+cartsizes[i])-1
			cart = cartmembs[fi:li]
			addcartel!(pmc.mat, cart, rand(sampledist))
			push!(cartels, cart)
			fi = li + 1
		end
		addrndlinks!(pm.mat, sum(pmc.mat) - sum(pm.mat))
		return pmrw, pm, pmc, cartels
	else
		return pmrw, pm, pm, cartels
	end
end
function rndpubnet2(k, pratio, commsizes, p_rewire=0.01,
									 cartsizes=nothing, p_carts=0.2, α=5, β=1)
	local pm
	first = true
	csizes = []
	cartels = Array{Array{Int, 1}, 1}()
	for cs in commsizes
		k0 = sample(k, cs, replace=false)
		p = Int(round(cs*pratio))
		pm0 = generate_publicationmatrix(k0, p)
		push!(csizes, size(pm0))
		if first
			pm = pm0
			first = false
		else
			m = combine_sparsematrices(pm.mat, pm0.mat)
			pm = generate_publicationmatrix(m)
		end
	end
	pmrw = rewire(pm, p_rewire)
	pm = deepcopy(pmrw)
	if !isnothing(cartsizes)
		pmc = deepcopy(pmrw)
		sampledist = Beta(α,β)
		minc = minimum(cartsizes)
		membersadded = 0
		while membersadded < size(pmc.mat,2) * p_carts
			sr = 1
			sc = 1
			for cs in csizes
				rand() > p_carts && continue
				mw = pm.mat[sr:(sr+cs[1]-1),sc:(sc+cs[2]-1)]
				mwc = pmc.mat[sr:(sr+cs[1]-1),sc:(sc+cs[2]-1)]
				local c
				na = size(mwc, 2)
				if minc <= na
					while true
						c = sample(cartsizes, 1)
						c = c[1]
						c <= na && break
					end
					cart = sample(1:na, Int(c), replace=false)
					addcartel!(mwc, cart, rand(sampledist))
					membersadded += c
					push!(cartels, cart .+ sc .- 1)
				end
				addrndlinks!(mw, sum(mwc)-sum(mw))
				pm.mat[sr:(sr+cs[1]-1),sc:(sc+cs[2]-1)] = mw
				pmc.mat[sr:(sr+cs[1])-1,sc:(sc+cs[2]-1)] = mwc
				sr += cs[1]
				sc += cs[2]
			end
		end
		return pmrw, pm, pmc, cartels
	else
		return pmrw, pm, pm, cartels
	end
end
function rndpubnet(k, pratio, commsizes, p_rewire=0.01,
									 cartsizes=nothing, p_carts=0.2, α=5, β=1)
	local pm
	first = true
	csizes = []
	cartels = Array{Array{Int, 1}, 1}()
	for cs in commsizes
		k0 = sample(k, cs, replace=false)
		p = Int(round(cs*pratio))
		pm0 = generate_publicationmatrix(k0, p)
		push!(csizes, size(pm0))
		if first
			pm = pm0
			first = false
		else
			m = combine_sparsematrices(pm.mat, pm0.mat)
			pm = generate_publicationmatrix(m)
		end
	end
	#println(csizes)
	pmrw = rewire(pm, p_rewire)
	pm = deepcopy(pmrw)
	if !isnothing(cartsizes)
		pmc = deepcopy(pmrw)
		sampledist = Beta(α,β)
		minc = minimum(cartsizes)
		sr = 1
		sc = 1
		for cs in csizes
			mw = pm.mat[sr:(sr+cs[1]-1),sc:(sc+cs[2]-1)]
			mwc = pmc.mat[sr:(sr+cs[1]-1),sc:(sc+cs[2]-1)]
			#println(size(mwc))
			local c
			na = size(mwc, 2)
			if minc <= na
				membersadded = 0
				while membersadded < na * p_carts
					while true
						c = sample(cartsizes, 1)
						c = c[1]
						c <= na && break
					end
					cart = sample(1:na, Int(c), replace=false)
					#addcartel!(mwc, cart, rand(sampledist))
					addcartel!(mwc, cart, rand(sampledist))
					membersadded += c
					push!(cartels, cart .+ sc .- 1)
				end
			end
			addrndlinks!(mw, sum(mwc)-sum(mw))
			#println(size(mw), ", ", size(mwc))
			pm.mat[sr:(sr+cs[1]-1),sc:(sc+cs[2]-1)] = mw
			pmc.mat[sr:(sr+cs[1])-1,sc:(sc+cs[2]-1)] = mwc
			sr += cs[1]
			sc += cs[2]
		end
		return pmrw, pm, pmc, cartels
	else
		return pmrw, pm, pm, cartels
	end
end
function rndpubnet2(k, pratio, commsizes, cartsizes, p_carts=0.2)
	sampledist = Beta(5,1)
	minc = minimum(cartsizes)
	first = true
	local pm
	local pmc
	for cs in commsizes
		k0 = sample(k, cs, replace=false)
		p = Int(round(cs*pratio))
		pm0 = generate_publicationmatrix(k0, p)
		local c
		na = size(pm0.mat, 2)
		if minc <= na
			membersadded = 0
			while membersadded < na * p_carts
				while true
					c = sample(cartsizes, 1)
					c = c[1]
					c <= na && break
				end
				cart = sample(1:na, Int(c), replace=false)
				pm0, pmc0 = pubmatpair(pm0, cart, rand(sampledist))
				#pm0, pmc0 = pubmatpair(pm0, cart, 1.0)
				membersadded += c
			end
		else
			pmc0 = deepcopy(pm0)
		end
		if first
			pm = pm0
			pmc = pmc0
			first = false
		else
			m = combine_sparsematrices(pm.mat, pm0.mat)
			pm = generate_publicationmatrix(m)
			m = combine_sparsematrices(pmc.mat, pmc0.mat)
			pmc = generate_publicationmatrix(m)
		end
	end
	r = Dict()
	r[:pm] = pm
	r[:pmc] = pmc
	return r
end

"""
    pubmatpair(pm, cart, probcart)

Duplicate `pm` add `cart` cartel to the duplum and then add random link
to `pm` to match the total link number to the publication matrix with
cartel
"""
function pubmatpair(pm::PubMat, cart::Array{Int,}, probcart::Float64)
	pmc = deepcopy(pm)
	n = sum(pm.mat)
	addcartel!(pmc, cart, probcart)
	nc = sum(pmc.mat)
	addrndlinks!(pm, nc-n)
	return pm, pmc
end

"""
    rndpubmat(ideal, auIDs, communities)

Generate a random publication network with similar characteristics as
`ideal`. The net is generated separatelly for all communities listed in
`communities`. `auIDs` is used to identify authors.
"""
function rndpubmat(ideal::PubMat, auIDs::Array{String,1},
									 communities::Dict{Int,Array{Int,1}})
	commIDs = keys(communities)
	first = true
	local m::SparseMatrixCSC{Int,Int}
	for cid in commIDs
		c = communities[cid]
		mid = auIDs[c]
		s = selectauthors(ideal, getauthorindex(ideal, mid))
		d = degreedegreecor(s)
		r = generaterndbipartite(d)
		if sum(r .< 0) > 0
			println("# members: ", length(c), " in communities ", cid,
							" ##########")
			error("Negative values!")
		end
		if first
			m = copy(r)
			first = false
		else
			m = combine_sparsematrices(m, r)
		end
	end
	counter = 0
	while sum(m .> 1) > 0 && counter < 10
		simplifybipartite!(m)
		counter += 1
	end
	counter == 10 && error("Simplification of generated publication matrix failed!")
	return generate_publicationmatrix(m)
end


"""
    rewire!(mat, niter)

Randomly rewire in-place the bipartite graph represented by the `mat`
adjacency matrix, `niter` times.
"""
function rewire!(mat::SparseMatrixCSC{Int64, Int64}, niter=100)
	countwiring = 0
	while countwiring < niter
		rowind, colind, vals = findnz(mat)
		n_entries = length(rowind)
		ri = randperm(n_entries)
		wi = 2
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
function rewire(pubmat::PubMat, niter::Int=100)
	matcp = deepcopy(pubmat)
	rewire!(matcp.mat, niter)
	return matcp
end
function rewire(pubmat::PubMat, piter::Float64=0.1)
	matcp = deepcopy(pubmat)
	nedges = nnz(matcp.mat)
	niter = Int(round(nedges * piter))
	niter == 0 && return pubmat
	rewire!(matcp.mat, niter)
	return matcp
end

"""
    rewire_community(pm, members, piter)

Randomly rewire the edges between members of a subgroup of publication
matrix `pm`. Only `piter` portion of the edges rewired.
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
    rewire_communities(pm, communities, prewire, file)

Randomly rewire the communities of publication matrix `pm`. Only
`prewire` portion of edges are rewired.
"""
function rewire_communities(pubmat::PubMat,
														communities::Dict{Int, Array{Int,1}},
														prewire::Float64, file="comm_rewire")
	pm = deepcopy(pubmat)
	for i in keys(communities)
		members = communities[i]
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
	m = spzeros(Int, length(ps_d), length(As_d))
	counter = 0
	while length(nedges) > 0 #&& counter < 100
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
			m[p,A] += 1
		end
		i = nedges .> 0
		pin = pin[i]
		Ain = Ain[i]
		nedges = nedges[i]
		counter += 1
	end
	counter = 0
	while sum(m .> 1) > 0 && counter < 25
		simplifybipartite!(m)
		counter += 1
	end
	return m
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
    simplifybipartite!(m)

Simplify a randomly generated bipartite graph represented by `m`.
Simplification means to resolve multiply edges between nodes. The
algorithm follows [this
paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5972839/)
"""
function simplifybipartite!(m::SparseMatrixCSC{Int64, Int64})
	p, A = size(m)
	ps = collect(1:p)
	As = collect(1:A)
	n_medges = sum(m .> 1)
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
		m[l,u] == 0 && continue
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
						m[l,u] -= 1
						m[v,w] -= 1
						m[l,w] = 1
						m[v,u] = 1
						if sum(m .< 0) > 0
							println("m[l,u] = ", m[l,u], ", m[v,w] = ", m[v,w])
							error("first")
						end
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
						m[l,u] -= 1
						m[v,w] -= 1
						m[l,w] = 1
						m[v,u] = 1
						if sum(m .< 0) > 0
							println("m[l,u] = ", m[l,u], ", m[v,w] = ", m[v,w])
							error("second")
						end
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
				if m[vp,up] <= 0 
					m[l,u] -= 1
					m[w,up] -= 1
					m[vp,wp] -= 1
					m[w,u] = 1
					m[l,wp] = 1
					m[vp,up] = 1
						if sum(m .< 0) > 0
							println("m[l,u] = ", m[l,u], ", m[w,up] = ", m[w,up],
											", m[vp,wp] = ", m[vp, wp])
							error("third")
						end
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

"""
    reshuffle(pm, nauthors)

Take a publication matrix `pm` and shuffle the publication records of
`nauthors` to randomise the matrix.
"""
function reshuffle(pm::SparseMatrixCSC{Int,Int}, nauthors::Int)
	p, A = size(pm)
	As = collect(1:A)
	np = papernumbers(pm)
	nps = reverse(sort(np))
	nps = nps[nauthors]
	rA = As[np .>= nps]
	println(length(rA))
	for a in rA
		pm[:,a] = shuffle(pm[:,a])
	end
	np = authornumbers(pm)
	pm = pm[np .> 0,:]
	dropzeros!(pm)
	return pm
end
function reshuffle(pm::PubMat, nauthors::Int)
	mat = copy(pm.mat)
	mat = reshuffle(mat, nauthors)
	return PubMat(mat, pm.authorIDs, pm.paperIDs)
end

"""
    generate_rndnet(k, p, cartel, pc)

Generate two random publication networks. One is a fully random network,
based on the distribution of papers published by the authors, `k`, and the
number of possible papers, `p`. The other network is the same as the
previous one except that authors listed in `cartel` form a cartel, i.e.
they invite each others to be coauthors with probability `pc`.
"""
function generate_rndnet(k, p, cartel::Array{Int64,1}, pc;
												 filename="proba", doplot=true, papercutoff=0)
	pnnc = Dict()
	pnwc = Dict()
	pnnc[:puma] = generate_publicationmatrix(k, p)
	papercutoff > 0 && (pnnc[:puma] = selectauthors(pnnc[:puma], papercutoff))
	pnnc[:coma] = collaborationmatrix(pnnc[:puma])
	pnnc[:coga] = collaborationgraph(pnnc[:coma])
	pnwc[:puma] = deepcopy(pnnc[:puma])
	addcartel!(pnwc[:puma], cartel, pc)
	papercutoff > 0 && (pnwc[:puma] = selectauthors(pnwc[:puma], papercutoff))
	pnwc[:coma] = collaborationmatrix(pnwc[:puma])
	pnwc[:coga] = collaborationgraph(pnwc[:coma])
	lx, ly = spring_layout(pnwc[:coga], C=20)
	if doplot
		graphplot(pnnc[:coga], lx, ly, filename=filename)
		graphplot(pnwc[:coga], lx, ly, filename=string(filename, "-cartel"))
	end
	return (pnnc, pnwc)
end
function generate_rndnet(k, p, cartel::Array{Array{Int64,1},1}, pc;
												 filename="proba", doplot=true, papercutoff=0)
	pnnc = Dict()
	pnwc = Dict()
	pnnc[:puma] = generate_publicationmatrix(k, p)
	papercutoff > 0 && (pnnc[:puma] = selectauthors(pnnc[:puma], papercutoff))
	pnnc[:coma] = collaborationmatrix(pnnc[:puma])
	pnnc[:coga] = collaborationgraph(pnnc[:coma])
	pnwc[:puma] = deepcopy(pnnc[:puma])
	for c in cartel
		addcartel!(pnwc[:puma], c, pc)
	end
	papercutoff > 0 && (pnwc[:puma] = selectauthors(pnwc[:puma], papercutoff))
	pnwc[:coma] = collaborationmatrix(pnwc[:puma])
	pnwc[:coga] = collaborationgraph(pnwc[:coma])
	lx, ly = spring_layout(pnwc[:coga], C=20)
	if doplot
		graphplot(pnnc[:coga], lx, ly, filename=filename)
		graphplot(pnwc[:coga], lx, ly, filename=string(filename, "-cartel"))
	end
	return (pnnc, pnwc)
end

"""
    generate_cartels(nauthors, cartel_sizes)

Generate an array of arrays containing members of cartels randomly drawn
from `nauthors` number of authors. The sizes of the generated cartels
are given by `cartel_sizes`, a `Dict`.
"""
function generate_cartels(nauthors, cartel_sizes)
	cartels = Array{Int64,1}[]
	nmembers = 0
	for k in keys(cartel_sizes)
		nmembers += k*cartel_sizes[k]
	end
	members = sample(1:nauthors, nmembers, replace=false)
	i = 1
	for k in keys(cartel_sizes)
		for s in 1:cartel_sizes[k]
			j = i + k
			push!(cartels, members[i:(j-1)])
			i = j
		end
	end
	return cartels
end

"""
    prop_cartels(nmembers, sample_cartel)

Make a cartel size distribution based on the distribution given in
`sample_cartel`. The generated cartel size distribution containse no
more than `nmembers` members.
"""
function prop_cartels(nmembers, sample_cartel)
	css = Int[]
	for k in keys(sample_cartel)
		for i in 1:sample_cartel[k]
			push!(css, k)
		end
	end
	shuffle!(css)
	cs = cumsum(css)
	css = css[cs .<= nmembers]
	return histogram(css)
end

"""
    addrndlink!(pm)

Add a random link between a paper and an author in publication matrix
`pm`.
"""
function addrndlink!(mat::SparseMatrixCSC{Int64, Int64})
	p, A = size(mat)
	ip = rand(1:p)
	iA = rand(1:A)
	while mat[ip,iA] == 1
		ip = rand(1:p)
		iA = rand(1:A)
	end
	mat[ip, iA] = 1
	return nothing
end
function addrndlink!(pm::PubMat)
	addrndlink!(pm.mat)
	return nothing
end

"""
"""
function addrndlinks!(pm, nlinks::Int)
	if typeof(pm) == PubMat
		nz = sum(pm.mat .== 0)
	else
		nz = sum(pm .== 0)
	end
	nlinks = minimum([nlinks, nz])
	for i in 1:nlinks
		addrndlink!(pm)
	end
	return nothing
end
