
# script to run publication network simulations to investigate the foot
# prints of cartels

include("../analyse-publications/PaperPump.jl")

# the footprint of cartels

function simulatepubnets()
	#carts = collect(2:5)
	carts = expectedpapers(1000, 20, 2.5, 0)
	carts .+= 1
	rnetsrw = Dict()
	rnetsrl = Dict()
	rnetsca = Dict()
	rnetscr = Dict() # network with cartels rewired randomly
	for i in 1:10
		cs = expectedpapers(500, 500, 2.5, 10) # community sizes
		#cs = sample(commsizes, 100, replace=false)
		k = saturatedexpectedpapers(Int(sum(cs)), 100) # number of papers
		rpmrw, rpm, rpmc = rndpubnet(k, 6.0, Int.(cs), 0.1, carts, 0.001)
		rpmrw3 = selectauthors(rpmrw, 3)
		rpm3 = selectauthors(rpm, 3)
		rpmc3 = selectauthors(rpmc, 3)
		rpmcr = rewire(rpmc3, 0.9)
		rpmcr3 = selectauthors(rpmcr, 3)
		rcmrw = collaborationmatrix(rpmrw3)
		rcm = collaborationmatrix(rpm3)
		rcmc = collaborationmatrix(rpmc3)
		rcmcr = collaborationmatrix(rpmcr3)
		rnetsrw[i] = PubNet(rpmrw, rcmrw, 0.5, "",
												"../analyse-publications/comm_detection.sh")
		rnetsrl[i] = PubNet(rpm, rcm, 0.5, "",
											"../analyse-publications/comm_detection.sh")
		rnetsca[i] = PubNet(rpmc, rcmc, 0.5, "",
													"../analyse-publications/comm_detection.sh")
		rnetscr[i] = PubNet(rpmcr, rcmcr, 0.5, "",
													"../analyse-publications/comm_detection.sh")
	end
	return rnetsrw, rnetsrl, rnetsca, rnetscr
end

function combineruns(pubnet, property)
	drw = getproperty(pubnet[1], property)
	for i in 2:length(pubnet)
		drw = vcat(drw, getproperty(pubnet[i], property))
	end
	drw = drw[drw .> 0]
	return drw
end

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

Random.seed!(101)

rnetsrw, rnetsrl, rnetsca, rnetscr = simulatepubnets()

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


