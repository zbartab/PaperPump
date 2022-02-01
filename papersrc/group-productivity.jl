
using HypothesisTests, GLM

include("../analyse-publications/PaperPump.jl")
include("../analyse-publications/CommunitiesNetworks.jl")

matplotlib.use("PDF")

function plprod(sbplt, cartdata, variable, rnddata, xlab="")
	cad = cartdata[:car][variable]
	maxx = (mx1 = maximum(cad)) < (mx2 = maximum(rnddata)) ? mx2 : mx1
	meanx = round(median(cad), digits=1)
	sbplt.hist(cad, density=true, label="tight groups, M = $meanx",
							alpha=0.25, bins=0:5:maxx, log=false)
	meanx = round(median(rnddata), digits=1)
	sbplt.hist(rnddata, density=true, label="random groups, M = $meanx",
							alpha=0.25, bins=0:5:maxx, log=false)
	sbplt.legend(fontsize=:small)
	sbplt.set_xlabel(xlab)
	sbplt.set_ylabel("density")
	ins = sbplt.inset_axes([0.65, 0.05, 0.3, 0.65])
	ins.boxplot(cad, positions=[1], showfliers=false, #labels=(""),
							patch_artist=true, boxprops = Dict("facecolor" => "C0",
																								 "alpha" => 0.25))
	ins.boxplot(rnddata, positions=[2], showfliers=false, #labels=(""),
							patch_artist=true, boxprops = Dict("facecolor" => "C1",
																								 "alpha" => 0.25))
	ins.set_xticklabels(labels="")
	ins.tick_params(axis="x", bottom=false)
	ins.tick_params(axis="y", labelsize=:small)
end


function toppers(cartdata, variable, rnddata, Q)
	cad = cartdata[:car][variable]
	return sum(cad .> quantile(rnddata, Q)) / length(cad)
end

function createDF(proddata)
	proddata[:car][:,:grtype] .= "tightgroup"
	nc = ncol(proddata[:raw])
	i = proddata[:raw][:measure] .== "group"
	rgrp = flat(Matrix(proddata[:raw][i, 1:(nc-1)]))
	i = proddata[:raw][:measure] .== "npaper"
	rnpa = flat(Matrix(proddata[:raw][i, 1:(nc-1)]))
	i = proddata[:raw][:measure] .== "wpaper"
	rwpa = flat(Matrix(proddata[:raw][i, 1:(nc-1)]))
	rgrsize = names(proddata[:raw])[1:(nc-1)]
	rgrsize = parse.(Int, rgrsize)
	rgrsize = flat(map(x -> repeat([x], Int(nrow(proddata[:raw])/3)), rgrsize))
	rdf = DataFrame(:groupsize => rgrsize, :groupprod => rgrp, :npapers => rnpa,
									:wpapers => rwpa)
	rdf[:,:grtype] = "random"
	resdf = vcat(proddata[:car], rdf)
	return(resdf)
end

function calcquantile(df::DataFrame, measure::Symbol)
	grtype = df.grtype .== "tightgroup"
	tgdf = df[grtype,:]
	rgdf = df[.!grtype,:]
	GS = unique(rgdf.groupsize)
	mD = Dict{Int, Array{Float64, 1}}()
	for n in GS
		mD[n] = rgdf[rgdf.groupsize .== n, measure]
	end
	Qs = zeros(nrow(tgdf))
	for i in 1:nrow(tgdf)
		m = tgdf[i, measure]
		Qs[i] = sum(m .>= mD[tgdf.groupsize[i]]) ./ length(mD[tgdf.groupsize[i]])
	end
	return Qs
end

pMTMT = Dict()
pMTMT[:raw] = CSV.read("../analyse-publications/MTMT/MTMTpubmat-productivity_raw.csv",
											 DataFrame)
pMTMT[:car] = CSV.read("../analyse-publications/MTMT/MTMTpubmat-cartel_productivity.csv",
											 DataFrame)

pdblp = Dict()
pdblp[:raw] = CSV.read("../analyse-publications/dblp/dblppubmat-productivity_raw.csv",
											 DataFrame)
pdblp[:car] = CSV.read("../analyse-publications/dblp/dblppubmat-cartel_productivity.csv",
											 DataFrame)

m = Matrix(pMTMT[:raw][pMTMT[:raw][:measure] .== "group",
											 1:(size(pMTMT[:raw], 2)-1)])
gMTMT = reshape(m, prod(size(m)))

m = Matrix(pMTMT[:raw][pMTMT[:raw][:measure] .== "npaper",
											 1:(size(pMTMT[:raw], 2)-1)])
nMTMT = reshape(m, prod(size(m)))

m = Matrix(pMTMT[:raw][pMTMT[:raw][:measure] .== "wpaper",
											 1:(size(pMTMT[:raw], 2)-1)])
wMTMT = reshape(m, prod(size(m)))

m = Matrix(pdblp[:raw][pdblp[:raw][:measure] .== "group",
											 1:(size(pdblp[:raw], 2)-1)])
gdblp = reshape(m, prod(size(m)))

m = Matrix(pdblp[:raw][pdblp[:raw][:measure] .== "npaper",
											 1:(size(pdblp[:raw], 2)-1)])
ndblp = reshape(m, prod(size(m)))

m = Matrix(pdblp[:raw][pdblp[:raw][:measure] .== "wpaper",
											 1:(size(pdblp[:raw], 2)-1)])
wdblp = reshape(m, prod(size(m)))

# number of authors in tight groups
sum(pMTMT[:car].groupsize)
sum(pdblp[:car].groupsize)

# number of tight groups
nrow(pMTMT[:car])
nrow(pdblp[:car])

sbplts = PyPlot.subplots(3,2,figsize=(7,7))
plprod(sbplts[2][1], pMTMT, :npapers, nMTMT, "mean number of papers")
plprod(sbplts[2][2], pMTMT, :groupprod, gMTMT, "mean group productivity")
plprod(sbplts[2][3], pMTMT, :wpapers, wMTMT, "mean weighted number of papers")
plprod(sbplts[2][4], pdblp, :npapers, ndblp, "mean number of papers")
plprod(sbplts[2][5], pdblp, :groupprod, gdblp, "mean group productivity")
plprod(sbplts[2][6], pdblp, :wpapers, wdblp, "mean weighted number of papers")
tight_layout()

savefig("../paperfigs/group_productivity.pdf")

dfMTMT = createDF(pMTMT)
dfdblp = createDF(pdblp)

#dfMTMT[:,:s_groupprod] = (dfMTMT.groupprod .- mean(dfMTMT.groupprod)) ./ std(dfMTMT.groupprod)
#dfdblp[:,:s_groupprod] = (dfdblp.groupprod .- mean(dfdblp.groupprod)) ./ std(dfdblp.groupprod)

dfMTMT[:,:s_groupprod] = (dfMTMT.groupprod .- mean(dfMTMT.groupprod))
dfdblp[:,:s_groupprod] = (dfdblp.groupprod .- mean(dfdblp.groupprod))

fmMTMT0 = lm(@formula(npapers ~ grtype * groupprod), dfMTMT)
fmMTMT1 = lm(@formula(wpapers ~ grtype * groupprod), dfMTMT)
fmdblp0 = lm(@formula(npapers ~ grtype * groupprod), dfdblp)
fmdblp1 = lm(@formula(wpapers ~ grtype * groupprod), dfdblp)

i = dfMTMT.groupsize .== 2
fmMTMT20 = lm(@formula(npapers ~ grtype * s_groupprod), dfMTMT[i, :])
fmMTMT21 = lm(@formula(wpapers ~ grtype * s_groupprod), dfMTMT[i, :])
i = dfdblp.groupsize .== 2
fmdblp20 = lm(@formula(npapers ~ grtype * s_groupprod), dfdblp[i, :])
fmdblp21 = lm(@formula(wpapers ~ grtype * s_groupprod), dfdblp[i, :])

figure(figsize=(12,5))
subplot(1,2,1)
newX = DataFrame(:groupprod =>
								 repeat([minimum(dfMTMT.groupprod), maximum(dfMTMT.groupprod)],
												2))
newX[:,:grtype] = vcat(repeat(["tightgroup"], Int(nrow(newX)/2)),
											 repeat(["random"], Int(nrow(newX)/2)))
#close("all")
i = dfMTMT.grtype .== "random"
PyPlot.scatter(dfMTMT.groupprod[i], dfMTMT.npapers[i], alpha=0.1,
							 label="random groups")
i = dfMTMT.grtype .== "tightgroup"
PyPlot.scatter(dfMTMT.groupprod[i], dfMTMT.npapers[i], alpha=0.1,
							 label="tight groups")
i = newX.grtype .== "random"
PyPlot.plot(newX.groupprod[i], predict(fmMTMT0, newX[i,:]),
							 label="random groups, estimates")
i = newX.grtype .== "tightgroup"
PyPlot.plot(newX.groupprod[i], predict(fmMTMT0, newX[i,:]),
						label="tight groups, estimates")
title("MTMT")
legend()
xlabel("average group productivity")
ylabel("average number of papers")
tight_layout()
#savefig("../work/screen.pdf")
subplot(1,2,2)
newX = DataFrame(:groupprod =>
								 repeat([minimum(dfdblp.groupprod), maximum(dfdblp.groupprod)],
												2))
newX[:,:grtype] = vcat(repeat(["tightgroup"], Int(nrow(newX)/2)),
											 repeat(["random"], Int(nrow(newX)/2)))
#close("all")
i = dfdblp.grtype .== "random"
PyPlot.scatter(dfdblp.groupprod[i], dfdblp.npapers[i], alpha=0.1,
							 label="random groups")
i = dfdblp.grtype .== "tightgroup"
PyPlot.scatter(dfdblp.groupprod[i], dfdblp.npapers[i], alpha=0.1,
							 label="tight groups")
i = newX.grtype .== "random"
PyPlot.plot(newX.groupprod[i], predict(fmdblp0, newX[i,:]),
							 label="random groups, estimates")
i = newX.grtype .== "tightgroup"
PyPlot.plot(newX.groupprod[i], predict(fmdblp0, newX[i,:]),
						label="tight groups, estimates")
title("dblp")
legend()
xlabel("average group productivity")
ylabel("average number of papers")
tight_layout()
savefig("../paperfigs/paper_sharing-relations.pdf")

# estimating c

i = dfMTMT.grtype .== "tightgroup"
#i = repeat([true], nrow(dfMTMT))
alfa = dfMTMT.npapers[i] ./ dfMTMT.groupprod[i]
cMTMT = (alfa .- 1) ./ (dfMTMT.groupsize[i] .- 1)
i = dfdblp.grtype .== "tightgroup"
#i = repeat([true], nrow(dfdblp))
alfa = dfdblp.npapers[i] ./ dfdblp.groupprod[i]
cdblp = (alfa .- 1) ./ (dfdblp.groupsize[i] .- 1)

describe(cMTMT)
describe(cdblp)

figure()
PyPlot.hist(cMTMT, alpha=0.25, density=true, bins=0.2:0.05:1, label="MTMT")
PyPlot.hist(cdblp, alpha=0.25, density=true, bins=0.2:0.05:1, label="dblp")
legend()
xlabel("level of authorship sharing")
ylabel("density")
tight_layout()
savefig("../paperfigs/paper_sharing.pdf")

mybins = 0:0.05:1
qnMTMT = calcquantile(dfMTMT, :npapers)
qndblp = calcquantile(dfdblp, :npapers)
close("all")
PyPlot.hist(qnMTMT, density=true, bins=mybins, label="MTMT", alpha=0.25)
PyPlot.hist(qndblp, density=true, bins=mybins, label="dblp", alpha=0.25)
xlabel("quantiles")
ylabel("density")
title("number of papers")
legend()
tight_layout()
savefig("../work/screen.pdf")

qgMTMT = calcquantile(dfMTMT, :groupprod);
qgdblp = calcquantile(dfdblp, :groupprod);
close("all")
PyPlot.hist(qgMTMT, density=true, bins=mybins, label="MTMT", alpha=0.25)
PyPlot.hist(qgdblp, density=true, bins=mybins, label="dblp", alpha=0.25)
xlabel("quantiles")
ylabel("density")
title("group productivity")
legend()
tight_layout()
savefig("../work/screen.pdf")

qwMTMT = calcquantile(dfMTMT, :wpapers);
qwdblp = calcquantile(dfdblp, :wpapers);
close("all")
PyPlot.hist(qwMTMT, density=true, bins=mybins, label="MTMT", alpha=0.25)
PyPlot.hist(qwdblp, density=true, bins=mybins, label="dblp", alpha=0.25)
xlabel("quantiles")
ylabel("density")
title("weighted number of papers")
legend()
tight_layout()
savefig("../work/screen.pdf")

# proportion of tight groups with significantly lower number of papers than
# random groups
sum(qnMTMT .<= 0.05) / length(qnMTMT)
sum(qndblp .<= 0.05) / length(qndblp)

# proportion of tight groups with significantly lower productivity than
# random groups
sum(qgMTMT .<= 0.05) / length(qgMTMT)
sum(qgdblp .<= 0.05) / length(qgdblp)

# proportion of tight groups with significantly lower weighted number of
# papers than random groups
sum(qwMTMT .<= 0.05) / length(qwMTMT)
sum(qwdblp .<= 0.05) / length(qwdblp)


