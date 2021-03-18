
# script to run publication network simulations to investigate the foot
# prints of cartels

include("../analyse-publications/PaperPump.jl")
include("../analyse-publications/CommunitiesNetworks.jl")
matplotlib.use("PDF")

using LinearAlgebra

# the footprint of cartels

# function definitions {{{1
function overlap(a::Array{Int, 1}, b::Array{Int, 1}) #{{{2
	return length(intersect(a,b))/length(union(a,b))
end

function AUROC(A, B, C) # {{{2
	TP = length(intersect(A, B))
	FP = length(setdiff(B, A))
	FN = length(setdiff(A, B))
	TN = length(setdiff(C, union(A, B)))
	TPR = TP/(TP+FN)
	FPR = FP/(FP+TN)
	recall = TP/(TP + FN)
	precision = TP/(TP + FP)
	Fvalue = (2*recall*precision)/(recall+precision)
	return TPR, FPR, Fvalue
end

function analysemembership(pn::PubNet, cartelsadded::Array{Array{Int, 1},1}, # {{{2
													 cond::Array{Float64, 1})
	cartelsfound = map(x -> getauthorindex(pn.coma, x), pn.cartels)
	allC = collect(1:nv(pn.coga))
	ocomm, f = findcommunities(pn.coga)
	ocomm = ocomm[length.(ocomm) .> 1]
	ocW = map(x -> median(Weights(pn.coga[x])), ocomm)
	ocCo = map(x -> coherence(pn, x), ocomm)
	ocEx = map(x -> exclusivity(pn, x), ocomm)
	ocomm = ocomm[(ocW .>= cond[3]) .& (ocCo .>= cond[1]) .& (ocEx .>= cond[2])]
	cfA = Array{Float64,1}()
	ocA = Array{Float64,1}()
	for ca in cartelsadded
		cfJ = 0.0
		for cf in cartelsfound
			J = overlap(ca, cf)
			cfJ < J && (cfJ = J)
		end
		push!(cfA, cfJ)
		ocJ = 0.0
		for oc in ocomm
			J = overlap(ca, oc)
			ocJ < J && (ocJ = J)
		end
		push!(ocA, ocJ)
	end
	return [length(unique(flat(cartelsadded))),
					length(unique(flat(cartelsfound))),
					length(unique(flat(ocomm))),
					mean(cfA), AUROC(flat(cartelsadded), flat(cartelsfound), allC)...,
					mean(ocA), AUROC(flat(cartelsadded), flat(ocomm), allC)...]
end

function titratecutoff(caD, rcD, cutoffs) #{{{2
	nnets = length(caD[:net])
	resM = zeros(Float64, (nnets*length(cutoffs), 7))
	j = 0
	for co in cutoffs
		for i in 1:nnets
			net = caD[:net][i]
			net.cartels = cartels(net.coga, co)
			r = analysemembership(net, rcD[:comm][i], [0.0, 0.0, 0.0])
			j += 1
			resM[j,:] .= flat([co, r[[1,2,4,5,6,7]]])
		end
	end
	return resM
end

function titratequantile(crD, caD, rcD, Qs) #{{{2
	nnets = length(caD[:net])
	resM = zeros(Float64, (nnets*length(Qs), 7))
	j = 0
	for Q in Qs
		for i in 1:nnets
			qCo = quantile(crD[:Co][i], Q)
			qEx = quantile(crD[:Ex][i], Q)
			r = analysemembership(caD[:net][i], rcD[:comm][i], [qCo, qEx, 0.0])
			j += 1
			resM[j,:] .= flat([Q, r[[1,2,8,9,10,11]]])
		end
	end
	return resM
end

""" {{{2
    flat(A)

Return the flattened version of the argument. Auxialary function.
"""
function flat(A)
	return collect(Iterators.Flatten(A))
end

""" {{{2
    simulatepubnets(nreps)

Simulate random publication networks `nreps` times. Return
	- the simulated networks, no cartels added
	- the simulated networks with as much random links added as in the net
	where cartels added
	- the simulated networks with cartels added
	- the simulated networks with cartels added randomly rewired
	- the cartels
"""
function simulatepubnets(ncomm=500, pratio=6.0, nreps=10, pcarts=0.01, α=5, β=1)
	#carts = collect(2:5)
	carts = expectedpapers(10000, 20, 2.5, 0)
	carts .+= 1
	rnetsrw = []
	rnetsrl = []
	rnetsca = []
	rnetscr = [] # network with cartels rewired randomly
	randcarts = []
	for i in 1:nreps
		cs = expectedpapers(ncomm, 5000, 2.5, 10) # community sizes
		#cs = sample(commsizes, 100, replace=false)
		k = saturatedexpectedpapers(Int(sum(cs)), 1000) # number of papers
		rpmrw, rpm, rpmc, rc = rndpubnet3(k, pratio, Int.(cs), 0.1, carts,
																		 pcarts, α, β)
		rpmrw3 = selectauthors(rpmrw, 3)
		rpm3 = selectauthors(rpm, 3)
		rpmc3 = selectauthors(rpmc, 3)
		rpmcr = rewire(rpmc3, 0.9)
		rpmcr3 = selectauthors(rpmcr, 3)
		rcmrw = collaborationmatrix(rpmrw3)
		rcm = collaborationmatrix(rpm3)
		rcmc = collaborationmatrix(rpmc3)
		rcmcr = collaborationmatrix(rpmcr3)
		push!(rnetsrw, PubNet(rpmrw, rcmrw, 0.5, "",
													"../analyse-publications/comm_detection.sh"))
		push!(rnetsrl, PubNet(rpm, rcm, 0.5, "",
												"../analyse-publications/comm_detection.sh"))
		push!(rnetsca, PubNet(rpmc, rcmc, 0.5, "",
													"../analyse-publications/comm_detection.sh"))
		push!(rnetscr, PubNet(rpmcr, rcmcr, 0.5, "",
													"../analyse-publications/comm_detection.sh"))
		caID = map(x -> getauthorIDs(rnetsca[i].puma, x), rc)
		push!(randcarts, map(x -> getauthorindex(rnetsca[i].coma, x), caID))
	end
	return rnetsrw, rnetsrl, rnetsca, rnetscr, randcarts
end

""" {{{2
    calcmeasures(pn, cliques)

Calculates several measures for the groups listed in `cliques`. Return
  - the median weights within group
	- size of the group
	- group productivity (number of distinc papers the group produced
	divided by group size)
	- average number of papers a group member produced
"""
function calcmeasures(pn::PubNet, cliques::Array{Array{Int64,1},1})
	cliW = Float64[]
	cliN = Float64[]
	cliP = Float64[]
	cliPN = Float64[]
	for c in cliques
		length(c) <= 1 && continue
		g = pn.coga[c]
		ww = Weights(g)
		#println(g)
		#println(ww)
		if length(ww) < 1
			push!(cliW, 0.0)
		else
			push!(cliW, median(ww))
		end
		push!(cliN, length(c))
		ids = map(x -> get_prop(g, x, :id), 1:nv(g))
		ai = getauthorindex(pn.puma, ids)
		push!(cliP, groupproductivity(pn.puma, ai))
		push!(cliPN, mean(pn.npapers[ai]))
	end
	return cliW, cliN, cliP, cliPN
end
function calcmeasures(pn::PubNet, cliques::Array{Array{String,1},1})
	println("#### 1")
	clis = map(x -> getauthorindex(pn.coma, x), cliques)
	println("#### 2")
	return calcmeasures(pn, clis)
end
function calcmeasuresindict(pns, cliques)
	cW = Array{Array{Float64, 1}, 1}()
	cN = Array{Array{Float64, 1}, 1}()
	cP = Array{Array{Float64, 1}, 1}()
	cPN = Array{Array{Float64, 1}, 1}()
	for i in keys(pns)
		W, N, P, PN = calcmeasures(pns[i], cliques[i])
		push!(cW, W)
		push!(cN, N)
		push!(cP, P)
		push!(cPN, PN)
	end
	return cW, cN, cP, cPN
end

""" {{{2
    calctopologyindict(pns, cliques)

Return local topology measures (coherence, exclusivity) for groups in
`cliques`.
"""
function calctopologyindict(pns, cliques)
	coh = Array{Array{Float64, 1}, 1}()
	exc = Array{Array{Float64, 1}, 1}()
	for i in keys(pns)
		c = map(x -> coherence(pns[i], x), cliques[i])
		e = map(x -> exclusivity(pns[i], x), cliques[i])
		push!(coh, c)
		push!(exc, e)
	end
	return coh, exc
end


""" {{{2
    findcommunitiesindict(dict)

Find the communities in networks stored in `dict`.
"""
function findcommunitiesindict(dict)
	commdict = Dict()
	for i in keys(dict)
		commdict[i], f = findcommunities(dict[i].coga)
	end
	return commdict
end

""" {{{2
    combineruns(pubnets, propoerty)

Return a combination of the given property for the nets in pubnets.
"""
function combineruns(pubnets, property)
	drw = getproperty(pubnets[1], property)
	for i in 2:length(pubnets)
		drw = vcat(drw, getproperty(pubnets[i], property))
	end
	drw = drw[drw .> 0]
	return drw
end

function calcsimmeasures(rnetsrw, rnetsca, rnetscr, randcarts) #{{{3
	# find communities {{{2
	commrw = map(x -> findcommunities(x.coga)[1], rnetsrw)
	commca = map(x -> findcommunities(x.coga)[1], rnetsca)
	commcr = map(x -> findcommunities(x.coga)[1], rnetscr)
	## remove communities of size 1 {{{3
	for i in 1:length(commrw)
		commrw[i] = commrw[i][length.(commrw[i]) .> 1]
		commca[i] = commca[i][length.(commca[i]) .> 1]
		commcr[i] = commcr[i][length.(commcr[i]) .> 1]
		randcarts[i] = randcarts[i][length.(randcarts[i]) .> 1]
	end
	# calculate measures {{{2
	## generated random publication network, "no-cartels", and its communities
	rwW, rwN, rwP, rwPN = calcmeasuresindict(rnetsrw, commrw)
	## cartels added to "no-cartels" to create "cartels-added", and its communities
	caW, caN, caP, caPN = calcmeasuresindict(rnetsca, commca)
	## "cartels-added" randomly rewired, "cartels-rewired", and its communities
	crW, crN, crP, crPN = calcmeasuresindict(rnetscr, commcr)
	## "cartels-added" network and the cartels
	rcW, rcN, rcP, rcPN = calcmeasuresindict(rnetsca, randcarts)

	# calculate topology measures {{{2
	rwCo, rwEx = calctopologyindict(rnetsrw, commrw)
	caCo, caEx = calctopologyindict(rnetsca, commca)
	crCo, crEx = calctopologyindict(rnetscr, commcr)
	rcCo, rcEx = calctopologyindict(rnetsca, randcarts)
	return Dict(:net => rnetsrw, :comm => commrw, :Co => rwCo, :Ex =>
							rwEx, :W => rwW, :N => rwN, :P => rwP, :PN => rwPN),
	Dict(:net => rnetsca, :comm => commca, :Co => caCo, :Ex => caEx, :W =>
			 caW, :N => caN, :P => caP, :PN => caPN),
	Dict(:net => rnetscr, :comm => commcr, :Co => crCo, :Ex => crEx, :W =>
			 crW, :N => crN, :P => crP, :PN => crPN),
	Dict(:net => rnetsca, :comm => randcarts, :Co => rcCo, :Ex => rcEx, :W
			 => rcW, :N => rcN, :P => rcP, :PN => rcPN)
end

function classify(caD, rcD; condCo = 0.9, condEx = 0.2, condW = 0.4) #{{{2
	### cartel conditions; communities fulfilling these conditions might be
	#cartels
	TPR = Array{Float64, 1}() # not found
	FPR = Array{Float64, 1}()
	Fvalue = Array{Float64, 1}()
	#for i in 1:length(caD[:Co])
	for i in keys(caD[:comm])
		ica = (caD[:Co][i] .> condCo) .& (caD[:Ex][i] .> condEx) .&
					(caD[:W][i] .> condW)
		### cartels found in communities
		commcai = caD[:comm][i]
		commrci = rcD[:comm][i]
		#println(sum(ica), " ", length(ica), " ", length(commcai))
		caca = commcai[ica] # fulfilling cartel conditions
		TN = length(setdiff(flat(commcai), flat(commrci)))
		TP = length(intersect(flat(caca), flat(commrci)))
		FP = length(setdiff(flat(caca), flat(commrci)))
		FN = length(setdiff(flat(commrci), flat(caca)))
		recall = TP/(TP + FN)
		precision = TP/(TP + FP)
		push!(TPR, TP/(TP + FN))
		push!(FPR, FP/(FP + TN))
		push!(Fvalue, (2*recall*precision)/(recall+precision))
	end
	return mean(TPR), mean(FPR), mean(Fvalue)
end

function bruteforceclassify(caD, rcD; nreps=1000) #{{{2
	TPR = Array{Float64, 1}()
	FPR = Array{Float64, 1}()
	Fva = Array{Float64, 1}()
	cCo = Array{Float64, 1}()
	cEx = Array{Float64, 1}()
	cW = Array{Float64, 1}()
	for ww in 1:nreps
		c = rand()
		e = rand()
		ww = rand()
		t, f, fv = classify(caD, rcD, condCo=c, condEx = e, condW = ww)
		push!(TPR, t)
		push!(FPR, f)
		push!(Fva, fv)
		push!(cCo, c)
		push!(cEx, e)
		push!(cW, ww)
	end
	return hcat(TPR, FPR, Fva, cCo, cEx, cW)
end

## plotting functions {{{3
""" {{{2
   plotres()

Plot the distributional results for randomly generated publication nets.
"""
function plotres(property::Symbol, xlab::String="", plotrndlink=false,
								 newplot=true, plotfun::Function=loglog)
	newplot && figure()
	for k in keys(rnetsrw)
		plotfun(eCCDF2(getproperty(rnetsrw[k], property))..., "--",
						alpha=0.25, color="magenta", ds="steps")
	end
	for k in keys(rnetsca)
		plotfun(eCCDF2(getproperty(rnetsca[k], property))..., "--",
						alpha=0.25, color="green", ds="steps")
	end
	if plotrndlink
		for k in keys(rnetscr)
			plotfun(eCCDF2(getproperty(rnetscr[k], property))..., "--",
							alpha=0.25, color="orange", ds="steps")
		end
	end
	drw = combineruns(rnetsrw, property)
	drl = combineruns(rnetscr, property)
	dca = combineruns(rnetsca, property)
	plotfun(eCCDF2(drw)..., "-", ds="steps", color="magenta",
					label="no-cartels")
	plotfun(eCCDF2(dca)..., "-", ds="steps", color="green",
					label="cartels added")
	plotrndlink && plotfun(eCCDF2(drl)..., "-", ds="steps", color="orange",
												 label="rewired")
	PyPlot.legend()
	PyPlot.ylabel("eCCDF")
	PyPlot.xlabel(xlab)
	PyPlot.grid()
	PyPlot.tight_layout()
end

""" {{{2
    plhist()

Plot the frequency distribution of community measures.
"""
function plhist(rwW, caW, crW, rcW, myxlabel; hbins=-0.05:0.1:1.05)
	figure()
	subplot(2,2,1)
	PyPlot.hist(rwW, bins=hbins)
	title("no-cartels")
	ylabel("frequency")
	subplot(2,2,2)
	PyPlot.hist(caW, bins=hbins)
	title("cartels-added")
	subplot(2,2,3)
	PyPlot.hist(crW, bins=hbins)
	title("cartels-rewired")
	xlabel(myxlabel); ylabel("frequency")
	subplot(2,2,4)
	PyPlot.hist(rcW, bins=hbins)
	title("cartels")
	xlabel(myxlabel)
	tight_layout()
end

function pl2Dhist(x1, x2, y1, y2, xlab, ylab; mybins=0:0.025:1) #{{{2
	figure(figsize=(9,4))
	subplot(1,2,1)
	PyPlot.hist2D(x1, y1, bins=mybins)
	title("cartels-added")
	xlabel(xlab)
	ylabel(ylab)
	colorbar()
	subplot(1,2,2)
	PyPlot.hist2D(x2, y2, bins=mybins)
	title("cartels")
	xlabel(xlab)
	colorbar()
	tight_layout()
end

# generate random nets and calculate simulation measures {{{1
Random.seed!(101)
rnetsrw, rnetsrl, rnetsca, rnetscr, randcarts = simulatepubnets(100, 6.0, 2, 0.1, 1, 1)
rwD, caD, crD, rcD = calcsimmeasures(rnetsrw, rnetsca, rnetscr, randcarts)

# plot distributions {{{1
plotres(:degrees, "degrees", true)
savefig("../paperfigs/cartel_footprint-degree.pdf")

plotres(:strengthes, "strength", true)
savefig("../paperfigs/cartel_footprint-strength.pdf")

plotres(:weights, "weight", true, true, semilogy)
savefig("../paperfigs/cartel_footprint-weight.pdf")

plotres(:npapers, "number of papers")
savefig("../paperfigs/cartel_footprint-npapers.pdf")

plotres(:wpapers, "weigthed number of papers")
savefig("../paperfigs/cartel_footprint-wpapers.pdf")

plotres(:nauthors, "number of authors", false, true, semilogy)
savefig("../paperfigs/cartel_footprint-nauthors.pdf")


# plot 2D distributions {{{1

hcr = fit(Histogram, (flat(crD[:Co]), flat(crD[:Ex])), (0:0.05:1, 0:0.05:1))
hca = fit(Histogram, (flat(caD[:Co]), flat(caD[:Ex])), (0:0.05:1, 0:0.05:1))

hcr = LinearAlgebra.normalize(hcr, mode=:probability)
hca = LinearAlgebra.normalize(hca, mode=:probability)

coQ = quantile(flat(crD[:Co]), 0.9)
exQ = quantile(flat(crD[:Ex]), 0.9)

myvmax = 0.025
figure(figsize=(12,5))
subplot(1,2,1)
PyPlot.pcolormesh(hcr.edges[1], hcr.edges[2], hcr.weights', vmax=myvmax)
xlabel("coherence")
ylabel("exclusivity")
title("degree preserving rewiring")
colorbar()
tight_layout()
subplot(1,2,2)
PyPlot.pcolormesh(hca.edges[1], hca.edges[2], hca.weights', vmax=myvmax)
colorbar()
PyPlot.scatter(flat(rcD[:Co]), flat(rcD[:Ex]), alpha=0.25, s=1,
							 color="red", label="cartels")
PyPlot.plot([0.95, 0.95, 1], [1, 0.2, 0.2], "--", color="white",
						label="cartels zone")
legend()
xlabel("coherence")
ylabel("exclusivity")
title("cartels added")
tight_layout()

savefig("../paperfigs/cartel_footprint-2Dhist.pdf")

# efficiency of Kertesz's algorithm {{{1

Qs = 0.6:0.05:1
rqM = titratequantile(crD, caD, rcD, Qs)
rqM = DataFrame(rqM)

aqM = aggregate(rqM, 1, mean)
sqM = aggregate(rqM, 1, std)

cutoffs = 0:0.05:0.95
rcM = titratecutoff(caD, rcD, cutoffs)
rcM = DataFrame(rcM)

acM = aggregate(rcM, 1, mean)
scM = aggregate(rcM, 1, std)

figure(figsize=(10,4))
subplot(1,2,1)
PyPlot.errorbar(Qs, aqM.x5_mean, sqM.x5_std, ms=3, fmt="o", linewidth=1,
								label="true positive rate")
PyPlot.errorbar(Qs, aqM.x6_mean, sqM.x6_std, ms=3, fmt="o", linewidth=1,
								label="false positive rate")
PyPlot.plot(Qs, aqM.x7_mean, linewidth=1, linestyle="--", label="F-measure")
ylim(-0.05,1.05)
legend()
title("Wasch & Kertész 2019 algorithm")
xlabel("quantiles")
ylabel("rates")
tight_layout()
subplot(1,2,2)
PyPlot.errorbar(cutoffs, acM.x5_mean, scM.x5_std, ms=3, fmt="o",
								linewidth=1, label="true positive rate")
PyPlot.errorbar(cutoffs, acM.x6_mean, scM.x6_std, ms=3, fmt="o",
								linewidth=1, label="false positive rate")
PyPlot.plot(cutoffs, acM.x7_mean, linewidth=1,  linestyle="--", 
						label="F-measure")
ylim(-0.05,1.05)
legend()
title("link weights")
xlabel("cut-offs")
#ylabel("rates")
tight_layout()

savefig("../paperfigs/cartel_footprint-AUROC.pdf")

# check robustness link weights {{{1

## probability of sharing an article {{{2

function classefficiencyalpha(alphas, nreps=20, cutoff=0.35) #{{{3
	M = zeros(Float64, (nreps*length(alphas), 12))
	j = 0
	for α in alphas
		rnetsrwa, rnetsrla, rnetscaa, rnetscra, randcartsa = simulatepubnets(50, 10.0, nreps, α, 1)
		rwDa, caDa, crDa, rcDa = calcsimmeasures(rnetsrwa, rnetscaa, rnetscra, randcartsa)
		for i in keys(caDa[:net])
			net = caDa[:net][i]
			net.cartels = cartels(net.coga, cutoff)
			mm = analysemembership(net, rcDa[:comm][i], [0.0,0.0,0.0])
			j += 1
			M[j,:] .= vcat(α, mm)
			#M[j,:] .= mm
		end
	end
	return M[:,[1,2,3,5,6,7,8]]
end

### plot {{{3

function mybeta(x, α, β)
	y = x .^ (α -1) .* (1 .- x) .^ (β -1)
	y ./= sum(y)
	return y
end

alphas = 1:0.5:15

alphaM = classefficiencyalpha(alphas)

alphaM = DataFrame(alphaM, [:alpha, :cartels, :foundcarts,
																	:overlap, :TPR, :TNR, :Fscore])
alphaMa = aggregate(alphaM, 1, mean)
alphaMs = aggregate(alphaM, 1, std)

figalpha = figure()
PyPlot.errorbar(alphas, alphaMa.TPR_mean, alphaMs.TPR_std, ms=3,
								fmt="o", linestyle="-", linewidth=1, alpha=0.5,
							 label="true positive rate")
PyPlot.errorbar(alphas, alphaMa.TNR_mean, alphaMs.TNR_std, ms=3,
								fmt="o", linestyle="-", linewidth=1, alpha=0.5,
								label="false positive rate")
PyPlot.errorbar(alphas, alphaMa.Fscore_mean, alphaMs.Fscore_std, ms=3,
								fmt="o", linestyle="--", linewidth=1, alpha=0.5,
								label="F-measure")
legend()
xlabel("α parameter of Beta distribution (β=1)")
ylabel(L"rates $\pm$ SD")
tight_layout()
ax2 = figalpha.add_axes([0.15, 0.2, 0.3, 0.3])
x = 0:0.01:1
y = mybeta(x, 2,1)
m = sum(x .* y)
ax2.plot(x, y, label=L"α =  2, $\bar{p}$ = %$(round(m, digits=3))")
y = mybeta(x, 4,1)
m = sum(x .* y)
ax2.plot(x, y, label=L"α =  4, $\bar{p}$ = %$(round(m, digits=3))")
y = mybeta(x, 6,1)
m = sum(x .* y)
ax2.plot(x, y, label=L"α =  6, $\bar{p}$ = %$(round(m, digits=3))")
y = mybeta(x, 12,1)
m = sum(x .* y)
ax2.plot(x, y, label=L"α = 12, $\bar{p}$ = %$(round(m, digits=3))")
ax2.legend(fontsize="xx-small")
xmin, xmax, ymin, ymax = ax2.axis(false)
#ax2.axvline(xmin, ymin, ymax, color="black", linewidth=1)
ax2.axvline(xmin, ymin, 1, color="black", linewidth=1)
ax2.axhline(ymin, xmin, xmax, color="black", linewidth=1)
ax2.set_title("distribution of sharing probability", fontsize="x-small")

savefig("../paperfigs/cartel_footprint-alpha.pdf")


## network density {{{2

function classefficiencyp(ps, nreps=20, cutoff=0.35) #{{{3
	M = zeros(Float64, (nreps*length(ps), 12))
	j = 0
	for p in ps
		rnetsrwa, rnetsrla, rnetscaa, rnetscra, randcartsa = simulatepubnets(50, p, nreps, 6, 1)
		rwDa, caDa, crDa, rcDa = calcsimmeasures(rnetsrwa, rnetscaa, rnetscra, randcartsa)
		for i in keys(caDa[:net])
			net = caDa[:net][i]
			net.cartels = cartels(net.coga, cutoff)
			mm = analysemembership(net, rcDa[:comm][i], [0.0,0.0,0.0])
			j += 1
			M[j,:] .= vcat(p, mm)
		end
	end
	return M[:,[1,2,3,5,6,7,8]]
end

### plot {{{3

ps = 1:2:21

pM = classefficiencyp(ps)

pM = DataFrame(pM, [:p, :cartels, :foundcarts,
																	:overlap, :TPR, :TNR, :Fscore])
pMa = aggregate(pM, 1, mean)
pMs = aggregate(pM, 1, std)

figp = figure()
PyPlot.errorbar(ps, pMa.TPR_mean, pMs.TPR_std, ms=3,
								fmt="o", linestyle="-", linewidth=1, alpha=0.5,
							 label="true positive rate")
PyPlot.errorbar(ps, pMa.TNR_mean, pMs.TNR_std, ms=3,
								fmt="o", linestyle="-", linewidth=1, alpha=0.5,
								label="false positive rate")
PyPlot.errorbar(ps, pMa.Fscore_mean, pMs.Fscore_std, ms=3,
								fmt="o", linestyle="--", linewidth=1, alpha=0.5,
								label="F-measure")
legend()
xlabel(L"network density, $p$")
ylabel(L"rates $\pm$ SD")
tight_layout()

savefig("../paperfigs/cartel_footprint-p.pdf")

exit()

# number of cartels as the function of cutoffs {{{1

function cartelsizes(pn, cutoffs=0.2:0.01:0.8) #{{{2
	ncarts = Array{Array{Int, 1},1}()
	for c in cutoffs
		carts = cartels(pn.coga, c)
		push!(ncarts, length.(carts))
	end
	return ncarts
end

## number of cartel members {{{2
cutoffs = 0.2:0.01:0.8

pnMTMT = PubNet("../analyse-publications/MTMT/MTMTpubmat.mat")
pndblp = PubNet("../analyse-publications/dblp/dblppubmat.mat")

ncMTMT = cartelsizes(pnMTMT, cutoffs)
ncdblp = cartelsizes(pndblp, cutoffs)

res = zeros(Float64, (96, 9))
j = 0
for pcarts in [0.01, 0.05, 0.1]
	for α in [1,3,5,10]
		for D in [1.0, 2.0, 4.0, 10.0]
			local r
			rnetsrw, rnetsrl, rnetsca, rnetscr, randcarts = simulatepubnets(100, D, 2, pcarts, α, 1)
			for i in 1:length(rnetsca)
				r = analysemembership(rnetsca[i], randcarts[i], [0.9, 0.9, 0.0])
			end
			j += 1
			res[j, : ] .= flat([pcarts, α, D, r[[1,2,4,5,6,7]]])
			ncs = []
			for net in rnetsca
				c = cartelsizes(net, cutoffs)
				push!(ncs, c)
			end
			ncsa = []
			for net in rnetsrl
				c = cartelsizes(net, cutoffs)
				push!(ncsa, c)
			end
			figure(figsize=(9, 9))
			subplot(2,2,1)
			plotres(:strengthes, "strength", true, false)
			subplot(2,2,2)
			plotres(:weights, "weight", true, false, semilogy)
			subplot(2,2,3)
			plotres(:degrees, "degrees", true, false)
			subplot(2,2,4)
			PyPlot.plot(cutoffs, log10.(sum.(ncMTMT) ./ nv(pnMTMT.coga)))
			PyPlot.plot(cutoffs, log10.(sum.(ncdblp) ./ nv(pndblp.coga)))
			for i in 1:length(ncs)
				PyPlot.plot(cutoffs, log10.(sum.(ncs[i]) ./ nv(rnetsca[i].coga)), "--", alpha=0.25)
				PyPlot.plot(cutoffs, log10.(sum.(ncsa[i]) ./ nv(rnetsrl[i].coga)), alpha=0.25)
			end
			PyPlot.title("D = $D, α = $α")
			tight_layout()
			savefig("../work/cartsizes-D$(D)-alpha$(α)-pcarts$(pcarts).pdf")
		end
	end
end

resdf = DataFrame(res, [:pcarts, :alpha, :D, :cartadded, :cartfound,
												:overlap, :TPR, :FPR, :Fval])

aggresdf = aggregate(resdf, [1,2,3], mean)

savefig("../work/screen.pdf")

PyPlot.plot(cutoffs, maximum.(ncMTMT), label="MTMT")
PyPlot.plot(cutoffs, maximum.(ncdblp), label="dblp")
PyPlot.plot(cutoffs, maximum.(ncarts), label="cartels added")
legend()
tight_layout()

PyPlot.plot(cutoffs, log10.(sum.(ncMTMT)), label="MTMT")
PyPlot.plot(cutoffs, log10.(sum.(ncdblp)), label="dblp")
PyPlot.plot(cutoffs, log10.(sum.(ncarts)), label="cartels added")
legend()
tight_layout()

PyPlot.plot(cutoffs, mean.(ncMTMT), label="MTMT")
PyPlot.plot(cutoffs, mean.(ncdblp), label="dblp")
PyPlot.plot(cutoffs, mean.(ncarts), label="cartels added")
legend()
tight_layout()


PyPlot.plot(cutoffs, log10.(mean.(ncMTMT)), label="MTMT")
PyPlot.plot(cutoffs, log10.(mean.(ncdblp)), label="dblp")
PyPlot.plot(cutoffs, log10.(mean.(ncarts)), label="cartels added")
legend()
tight_layout()


# sandbox {{{1

## eCCDF {{{2

figure(figsize=(11,6))
subplot(2,3,1)
PyPlot.loglog(eCCDF3(rnetsrw[2].strengthes)..., ds="steps", alpha=0.5, label="no cartels")
#PyPlot.loglog(eCCDF3(rnetsrl[2].strengthes)..., ds="steps", alpha=0.5, label="random links")
PyPlot.loglog(eCCDF3(rnetsca[2].strengthes)..., ds="steps", alpha=0.5, label="cartels added")
PyPlot.loglog(eCCDF3(pnMTMT.strengthes)..., ds="steps", alpha=0.5, label="MTMT")
PyPlot.loglog(eCCDF3(pndblp.strengthes)..., ds="steps", alpha=0.5, label="dblp")
xlabel("strengths")
ylabel("eCCDF")
legend(fontsize="small")
tight_layout()
subplot(2,3,2)
PyPlot.loglog(eCCDF3(rnetsrw[2].degrees)..., ds="steps", alpha=0.5, label="no cartels")
#PyPlot.loglog(eCCDF3(rnetsrl[2].degrees)..., ds="steps", alpha=0.5, label="random links")
PyPlot.loglog(eCCDF3(rnetsca[2].degrees)..., ds="steps", alpha=0.5, label="cartels added")
PyPlot.loglog(eCCDF3(pnMTMT.degrees)..., ds="steps", alpha=0.5, label="MTMT")
PyPlot.loglog(eCCDF3(pndblp.degrees)..., ds="steps", alpha=0.5, label="dblp")
xlabel("degrees")
ylabel("eCCDF")
legend(fontsize="small")
tight_layout()
subplot(2,3,3)
PyPlot.loglog(eCCDF3(rnetsrw[2].npapers)..., ds="steps", alpha=0.5, label="no cartels")
#PyPlot.loglog(eCCDF3(rnetsrl[2].npapers)..., ds="steps", alpha=0.5, label="random links")
PyPlot.loglog(eCCDF3(rnetsca[2].npapers)..., ds="steps", alpha=0.5, label="cartels added")
PyPlot.loglog(eCCDF3(pnMTMT.npapers)..., ds="steps", alpha=0.5, label="MTMT")
PyPlot.loglog(eCCDF3(pndblp.npapers)..., ds="steps", alpha=0.5, label="dblp")
xlabel("number of papers")
ylabel("eCCDF")
legend(fontsize="small")
tight_layout()
subplot(2,3,4)
PyPlot.semilogy(eCCDF2(rnetsrw[2].weights)..., ds="steps", alpha=0.5, label="no cartels")
#PyPlot.semilogy(eCCDF2(rnetsrl[2].weights)..., ds="steps", alpha=0.5, label="random links")
PyPlot.semilogy(eCCDF2(rnetsca[2].weights)..., ds="steps", alpha=0.5, label="cartels added")
PyPlot.semilogy(eCCDF2(pnMTMT.weights)..., ds="steps", alpha=0.5, label="MTMT")
PyPlot.semilogy(eCCDF2(pndblp.weights)..., ds="steps", alpha=0.5, label="dblp")
xlabel("weights")
ylabel("eCCDF")
legend(fontsize="small")
tight_layout()
subplot(2,3,5)
PyPlot.loglog(eCCDF3(rnetsrw[2].clustcoefs)..., ds="steps", alpha=0.5, label="no cartels")
#PyPlot.loglog(eCCDF3(rnetsrl[2].clustcoefs)..., ds="steps", alpha=0.5, label="random links")
PyPlot.loglog(eCCDF3(rnetsca[2].clustcoefs)..., ds="steps", alpha=0.5, label="cartels added")
PyPlot.loglog(eCCDF3(pnMTMT.clustcoefs)..., ds="steps", alpha=0.5, label="MTMT")
PyPlot.loglog(eCCDF3(pndblp.clustcoefs)..., ds="steps", alpha=0.5, label="dblp")
xlabel("clustering coefficients")
ylabel("eCCDF")
legend(fontsize="small")
tight_layout()
subplot(2,3,6)
PyPlot.loglog(eCCDF3(rnetsrw[1].nauthors)..., ds="steps", alpha=0.5, label="no cartels")
#PyPlot.loglog(eCCDF3(rnetsrl[1].nauthors)..., ds="steps", alpha=0.5, label="random links")
PyPlot.loglog(eCCDF3(rnetsca[1].nauthors)..., ds="steps", alpha=0.5, label="cartels added")
PyPlot.loglog(eCCDF3(pnMTMT.nauthors)..., ds="steps", alpha=0.5, label="MTMT")
PyPlot.loglog(eCCDF3(pndblp.nauthors)..., ds="steps", alpha=0.5, label="dblp")
xlabel("number of authors")
ylabel("eCCDF")
legend(fontsize="small")
tight_layout()
savefig("../work/screen.pdf")
