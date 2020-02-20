
using StatsBase, Statistics

# function to analyses publication networks

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
    meanpapers(pm, group)

Calculate the average number of papers, and their standard deviation,
authored by people in `group`.
"""
function meanpapers(pm, group, f::Function=papernumbers)
	pn = f(pm)
	(mean(pn[group]), std(pn[group]))
end

"""
    rndnonmembers(people, group, n)

Select `n` random members from `people` who are not part of `group`.
"""
function rndnonmembers(people, group, n=length(group))
	d = setdiff(people, group)
	nonmembers = sample(d, n, replace=false)
	return nonmembers
end

"""
    papersdist(pm, group, nrep)

Generate a distribution of average number of papers by people not in
`group`.
"""
function papersdist(pm, group, nrep::Int=1000, f::Function=papernumbers)
	A = collect(1:length(pm.authorIDs))
	mpapers = []
	for i in 1:nrep
		nongroup = rndnonmembers(A, group)
		m, s = meanpapers(pm, nongroup, f)
		push!(mpapers, m)
	end
	return mpapers
end

"""
    weightedpapers(pm)

Calculate for each author in publication matrix `pm` the number of
papers weighted by the inverse of the author numbers. 
"""
function weightedpapers(pm::PubMat)
	na = authornumbers(pm)
	row, col, val = findnz(pm.mat)
	for i in 1:length(row)
		val[i] /= na[row[i]]
	end
	mb = sparse(row, col, val, size(pm.mat)...)
	return papernumbers(mb)
end

"""
    eCCDF(x)

Calculate the empirical complementary cummulative distribution function
of `x`.
"""
function eCCDF(x::Array{Float64, 1})
	hx = histogram(x)
	xs = sort(collect(keys(hx)))
	y = Int[]
	for k in xs
		push!(y, hx[k])
	end
	y = 1 .- (cumsum(y) ./ sum(y))
	y = vcat(1, y[1:(end-1)])
	return xs, y
end
function eCCDF(x::Array{Int, 1})
	hx = histogram(x)
	xs = sort(collect(keys(hx)))
	y = Int[]
	for k in xs
		push!(y, hx[k])
	end
	y = 1 .- (cumsum(y) ./ sum(y))
	y = vcat(1, y[1:(end-1)])
	return xs, y
end
function eCCDF2(x::Array{Int, 1})
	xs = sort(x)
	y = collect(length(xs):-1:1) ./ length(xs)
	#y = vcat(1, y[1:(end-1)])
	return xs, y
end
function eCCDF2(x::Array{Float64, 1})
	xs = sort(x)
	y = collect(length(xs):-1:1) ./ length(xs)
	#y = vcat(1, y[1:(end-1)])
	return xs, y
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
    pubnetstats(pn)

Calculate several statistics for the authors of publication net `pn`.
"""
function pubnetstats(pn)
	pubstats = Dict()
	pubstats[:npapers] = papernumbers(pn[:puma])
	pubstats[:wpapers] = weightedpapers(pn[:puma])
	pubstats[:weights] = Weights(pn[:coga])
	pubstats[:degree] = degree(pn[:coga])
	return pubstats
end

"""
    rewire(pn, nrand)

Randomly rewire, `nrand` times, the publication network `pn` (a
collection of publication and collaboration matrices and a collaboration
graph).
"""
function rewire(pn::Dict{Any,Any}, nrand::Int=1000)
	rpm = rewire(pn[:puma], nrand)
	rcm = collaborationmatrix(rpm)
	rcg = collaborationgraph(rcm)
	return Dict(:puma => rpm, :coma => rcm, :coga => rcg)
end
