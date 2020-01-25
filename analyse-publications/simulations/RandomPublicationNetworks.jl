
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
	pubnet = spzeros(p, A)
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
    showpubmat)pubmat)

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
