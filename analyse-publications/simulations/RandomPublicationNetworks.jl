
# Functions used to simulate publication networks and cartels

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

