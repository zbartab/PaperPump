
# Functions used to simulate publication networks and cartels

using StatsBase

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
				pm[npapers, j] = 1.0
				papers[npapers] = 1
			else
				# join existing paper
				i = 0
				ii = 0
				while true
					ii += 1
					i = papertojoin(pk)
					(pm[i,j] == 0.0 || ii > npapers) && break
				end
				if ii > npapers
					npapers += 1
					pm[npapers, j] = 1.0
					papers[npapers] = 1
				else
					pm[i,j] = 1.0
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
function sim_collaborationnetwork(A, cutoff=0.4, gamma=2.5, k_sat=5,
																		maxpapers=100_000)
	cn = rnd_collaborationnetwork(A, gamma, k_sat, maxpapers, 100000)
	return describecartels(cn, cutoff)
end

"""
    run_sims(A, n_iter, cutoff, gamma, k_sat, maxpapers)

Run a series of simulation to collect data on cartels over `n_iter`
random collaboration network.
"""
function run_sims(A; n_iter=10, cutoff=0.4, gamma=2.5, k_sat=5,
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
function addcartel!(pubmat::PubMat, cartel::Array{Int,1}, collprob=1.0)
	cartpub = pubmat.mat[:, cartel]
	i = sum(cartpub, dims=2) .> 0
	i = reshape(i, size(pubmat.mat, 1))
	j = collect(1:size(pubmat.mat, 1))[i]
	for jj in j, c in cartel
		if pubmat.mat[jj,c] <= 0 && rand() < collprob
			pubmat.mat[jj,c] = 1
		end
	end
	#pubmat.mat[i, cartel] .= 1
	return nothing
end

"""
    showpubmat(pubmat)

Visualise the collaboration graph created from the publication matrix
`pubmat`.
"""
function showpubmat(pubmat, cutoff=0.4)
	coma = collaborationmatrix(pubmat)
	coga = collaborationgraph(coma)
	d = describecartels(coga, cutoff)
	for k in keys(d)
		println(k, ": ", d[k])
	end
	graphplot(coga)
	hist(Weights(coga), 0:0.025:1)
	return nothing
end

"""
"""
function rndpubnet(k,pratio,commsizes)
	cs = commsizes[1]
	k0 = sample(k, cs, replace=false)
	p = Int(round(cs*pratio))
	pm = generate_publicationmatrix(k0, p)
	for cs in commsizes[2:end]
		k0 = sample(k, cs, replace=false)
		p = Int(round(cs*pratio))
		pm0 = generate_publicationmatrix(k0, p)
		m = combine_sparsematrices(pm.mat, pm0.mat)
		pm = generate_publicationmatrix(m)
	end
	return pm
end
