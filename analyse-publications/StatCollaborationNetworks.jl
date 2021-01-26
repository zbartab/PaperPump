
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
	v = Array{Float64,1}(undef, length(row))
	for i in 1:length(row)
		v[i] = val[i] / na[row[i]]
	end
	mb = sparse(row, col, v, size(pm.mat)...)
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
function eCCDF3(x::Array{Int, 1}, normfun::Function=mean)
	xx = x ./ normfun(x)
	xs = sort(xx)
	y = collect(length(xs):-1:1) ./ length(xs)
	#y = vcat(1, y[1:(end-1)])
	return xs, y
end
function eCCDF3(x::Array{Float64, 1}, normfun::Function=mean)
	xx = x ./ normfun(x)
	xs = sort(xx)
	y = collect(length(xs):-1:1) ./ length(xs)
	#y = vcat(1, y[1:(end-1)])
	return xs, y
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
function describecartels(colnet, cutoff=0.5)
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
"""
function cartels(colgraph::MetaGraph{Int64, Float64}, cutoff::Float64)
	cartsn = subnet(colgraph, cutoff)
	ca = connected_components(cartsn)
	caIDs = Array{String, 1}[]
	for c in ca
		push!(caIDs, getauthorIDs(cartsn, c))
	end
	return caIDs
end

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
function papernumbers(mat::SparseMatrixCSC{Int64,Int64})
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
function authornumbers(mat::SparseMatrixCSC{Int64,Int64})
	no_authors = sum(mat, dims=2)
	no_authors = reshape(no_authors, size(no_authors, 1))
	return no_authors
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
function degreesdegrees(pm::PubMat)
	da = Int.(papernumbers(pm))
	dp = Int.(authornumbers(pm))
	p, A, val = findnz(pm.mat)
	dhp = Int[]
	dhA = Int[]
	for i in 1:length(val)
		push!(dhp, dp[p[i]])
		push!(dhA, da[A[i]])
	end
	return dhp, dhA
end

# Productivity

"""
    groupproductivity(pm, group)

Return the per capita group productivity of `group` members in publication
matrix `pm`. Group productivity is the total number of papers the
authors in group wrote.
"""
function groupproductivity(pm::PubMat, group::Array{Int,1})
	gpm = selectauthors(pm, group)
	np = size(gpm)
	np = np[1]
	return np / length(group)
end
function groupproductivity(pm::PubMat, group::Array{String,1})
	inds = getauthorindex(pm, group)
	gp = groupproductivity(pm, inds)
	#gp < 2 && length(inds) == 2 && error("$(gp), $(inds)")
	return gp
end

"""
    convertgpQ(pn, groups)

Convert the productivity values of groups, given in `groups` to
percentile values generated from random groups.
"""
function convertgpQ(pn::PubNet, groups::Array{Array{String,1} ,1})
	gp = groupsproductivity(pn, groups)
	l = length(gp.groupsize)
	groupprod = Float64[]
	npapers = Float64[]
	wpapers = Float64[]
	for i in 1:l
		push!(groupprod, getQ(gp.groupprod[i],
													pn.prodsampleQ[:g][!,Symbol(gp.groupsize[i])]))
		push!(npapers, getQ(gp.npapers[i],
													pn.prodsampleQ[:n][!,Symbol(gp.groupsize[i])]))
		push!(wpapers, getQ(gp.wpapers[i],
													pn.prodsampleQ[:w][!,Symbol(gp.groupsize[i])]))
	end
	return (groupsize=gp.groupsize, groupprod=groupprod, npapers=npapers,
					wpapers=wpapers)
end

"""
    getQ(v, quants)

		Return the probability value for a given observed value `v`. The
		probability value is obtained as largest probability in QUANTS for
		which it is true that `v` is greater than the given quantile.
"""
function getQ(v, quants)
	i = sum(quants .<= v)
	if i == 0
		return 0.0
	else
		return QUANTS[i]
	end
end

"""
    QQvals(obsQ, theoP)

Return the proportion of values faling between probabilities in theoP.
"""
function QQvals(obsQ, theoP)
	n = countmap(obsQ)
	z = zeros((length(n),2))
	i = 1
	for k in sort(collect(keys(n)))
		z[i,1] = k
		z[i,2] = n[k]
		i += 1
	end
	ztot = sum(z[:,2])
	z[:,2] = cumsum(z[:,2])
	z[:,2] ./= ztot
end

"""
    groupsproductivity(pn, groups)

Return group size, group productivity, mean number of papers and
mean weighted papers of individuals in each group in `groups`.
"""
function groupsproductivity(pn::PubNet,
														groups::Array{Array{String,1}, 1})
	gsize = Int[]
	m_npapers = Float64[]
	m_wpapers = Float64[]
	m_gpapers = Float64[]
	for g in groups
		membersi = getauthorindex(pn.puma, g)
		push!(gsize, length(membersi))
		push!(m_npapers, mean(pn.npapers[membersi]))
		push!(m_wpapers, mean(pn.wpapers[membersi]))
		push!(m_gpapers, groupproductivity(pn.puma, membersi))
	end
	return (groupsize=gsize, groupprod=m_gpapers, npapers=m_npapers,
					wpapers=m_wpapers)
end
function groupsproductivity(pn::PubNet,
														groups::Dict{Int64, Array{String,1}})
	gsize = Int[]
	m_npapers = Float64[]
	m_wpapers = Float64[]
	m_gpapers = Float64[]
	for k in keys(groups)
		membersi = getauthorindex(pn.puma, groups[k])
		push!(gsize, length(membersi))
		push!(m_npapers, mean(pn.npapers[membersi]))
		push!(m_wpapers, mean(pn.wpapers[membersi]))
		push!(m_gpapers, groupproductivity(pn.puma, membersi))
	end
	return (groupsize=gsize, groupprod=m_gpapers, npapers=m_npapers,
					wpapers=m_wpapers)
end

"""
    samplingproductivityraw(pn, groups, nrep)

Calculate the productivity for random samples drawn from individuals not
members of groups. Return raw walues.
"""
function samplingproductivityraw(pm::PubMat, cm::ColMat, npapers::Array{Int,1},
															wpapers::Array{Float64,},
															groups::Array{Array{String,1},1}, nrep::Int=10)
	members = collect(Iterators.flatten(groups))
	nonmembers = setdiff(collect(keys(cm.authorIDs)), members)
	groupsizes = histogram(map(length, groups))
	gnull = Dict{Int, Array{Float64,1}}()
	nnull = Dict{Int, Array{Float64,1}}()
	wnull = Dict{Int, Array{Float64,1}}()
	for s in keys(groupsizes)
		gnull[s] = Float64[]
		nnull[s] = Float64[]
		wnull[s] = Float64[]
		for i in 1:nrep
			nonc = sample(nonmembers, s, replace=false)
			#println(nonc)
			mi = getauthorindex(pm, nonc)
			#println(mi)
			push!(gnull[s], groupproductivity(pm, nonc))
			push!(nnull[s], mean(npapers[mi]))
			push!(wnull[s], mean(wpapers[mi]))
		end
	end
	return gnull, nnull, wnull
end
function samplingproductivityraw(pn::PubNet,
														 groups::Array{Array{String,1},1}, nrep=10)
	return samplingproductivityraw(pn.puma, pn.coma, pn.npapers, pn.wpapers,
															groups, nrep)
end
function samplingproductivityraw(filename::String, pn::PubNet,
														 groups::Array{Array{String,1},1}, nrep=10)
	g, n, w = samplingproductivityraw(pn.puma, pn.coma, pn.npapers, pn.wpapers,
															groups, nrep)
	prodres = vcat(DataFrame(g), DataFrame(n), DataFrame(w))
	prodres[!,:measure] = vcat(repeat(["group"], Int(nrow(prodres)/3)),
														repeat(["npaper"], Int(nrow(prodres)/3)),
													 repeat(["wpaper"], Int(nrow(prodres)/3)))
	CSV.write(filename, prodres)
end

"""
    relativeproductivity(pn, groupprod)

Compare the productivity of groups to randomly chosen authors not in
groups.
"""
function relativeproductivity(pn::PubNet, groupprod, quantval=0.5)
	q50 = QUANTS .== quantval
	gs = Float64[]
	ns = Float64[]
	ws = Float64[]
	for i in 1:length(groupprod.groupsize)
		c = groupprod.groupsize[i]
		gnc = pn.prodsampleQ[:g][q50, Symbol(c)][1]
		nnc = pn.prodsampleQ[:n][q50, Symbol(c)][1]
		wnc = pn.prodsampleQ[:w][q50, Symbol(c)][1]
		push!(gs, groupprod.groupprod[i]/gnc)
		push!(ns, groupprod.npapers[i]/nnc)
		push!(ws, groupprod.wpapers[i]/wnc)
	end
	return (groupsize=groupprod.groupsize, groupprod=gs, npapers=ns,
					wpapers=ws)
end

# Local network measures of groups

"""
    coherence(pn, group)

Calculate the coherence of group `group` in publication network `pn`.
Coherence calculation follows Wachs & Kertész 2019 SciRep 9:10818.
"""
function coherence(pn::PubNet, gi::Array{Int,1})
	gsize = length(gi)
	W = Float64[]
	for i in 1:(gsize-1)
		for j in (i+1):gsize
			has_edge(pn.coga, gi[i], gi[j]) && push!(W,
																							 get_prop(pn.coga, gi[i], gi[j],
																												:weight))
		end
	end
	if length(W) < 1
		return 0.0
	else
		am = mean(W)
		gm = geomean(W)
		return gm/am
	end
end
function coherence(pn::PubNet, group::Array{String,1})
	gi = getauthorindex(pn.coma, group)
	return coherence(pn, gi)
end

"""
    exclusivity(pn, group)

Calculate the exclusivity of group `group` in publication network `pn`.
Exclusivity calculation follows Wachs & Kertész 2019 SciRep 9:10818.
"""
function exclusivity(pn::PubNet, gi::Array{Int, 1})
	gsize = length(gi)
	Win = Float64[]
	for i in 1:(gsize-1)
		for j in (i+1):gsize
			has_edge(pn.coga, gi[i], gi[j]) && push!(Win,
																							 get_prop(pn.coga, gi[i], gi[j],
																												:weight))
		end
	end
	Wout = Float64[]
	for i in gi
		N = neighbors(pn.coga, i)
		for j in N
			!in(j, gi) && push!(Wout, get_prop(pn.coga, i, j, :weight))
		end
	end
	return sum(Win)/(sum(Win)+sum(Wout))
end
function exclusivity(pn::PubNet, group::Array{String, 1})
	gi = getauthorindex(pn.coma, group)
	return exclusivity(pn, gi)
end

"""
    groupstrength(pn, group)

Calculate the average strength of links within `group` in publication
network `pn`.
"""
function groupstrength(pn::PubNet, group::Array{String, 1})
	gi = getauthorindex(pn.coma, group)
	gsize = length(gi)
	Win = Float64[]
	for i in 1:(gsize-1)
		for j in (i+1):gsize
			has_edge(pn.coga, gi[i], gi[j]) && push!(Win,
																							 get_prop(pn.coga, gi[i], gi[j],
																												:weight))
		end
	end
	return sum(Win)/gsize
end

"""
    density(m)

Return the density of links in a publication matrix.
"""
function density(m::SparseMatrixCSC{Int,Int})
	return sum(m)/prod(size(m))
end
function density(pm::PubMat)
	return density(pm.mat)
end
function density(pn::PubNet)
	return density(pn.puma.mat)
end
