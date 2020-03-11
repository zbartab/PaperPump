
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
	for i in 1:10
		cs = expectedpapers(100, 500, 2.5, 10) # community sizes
		#cs = sample(commsizes, 100, replace=false)
		k = saturatedexpectedpapers(Int(sum(cs)), 100) # number of papers
		rpmrw, rpm, rpmc = rndpubnet(k, 6.0, Int.(cs), 0.1, carts, 0.001)
		rpmrw3 = selectauthors(rpmrw, 3)
		rpm3 = selectauthors(rpm, 3)
		rpmc3 = selectauthors(rpmc, 3)
		rcmrw = collaborationmatrix(rpmrw3)
		rcm = collaborationmatrix(rpm3)
		rcmc = collaborationmatrix(rpmc3)
		rnetsrw[i] = PubNet(rpmrw, rcmrw, 0.5, "",
												"../analyse-publications/comm_detection.sh")
		rnetsrl[i] = PubNet(rpm, rcm, 0.5, "",
											"../analyse-publications/comm_detection.sh")
		rnetsca[i] = PubNet(rpmc, rcmc, 0.5, "",
													"../analyse-publications/comm_detection.sh")
	end
	return rnetsrw, rnetsrl, rnetsca
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
						alpha=0.25, color="orange", ds="steps")
	end
	if plotrndlink
		for k in keys(rnetsrl)
			plotfun(eCCDF2(getproperty(rnetsrl[k], property))..., "--",
							alpha=0.25, color="blue", ds="steps")
		end
	end
	for k in keys(rnetsca)
		plotfun(eCCDF2(getproperty(rnetsca[k], property))..., "--",
						alpha=0.25, color="green", ds="steps")
	end
	drw = combineruns(rnetsrw, property)
	drl = combineruns(rnetsrl, property)
	dca = combineruns(rnetsca, property)
	plotfun(eCCDF2(drw)..., "-", ds="steps", color="orange", label="random")
	plotrndlink && plotfun(eCCDF2(drl)..., "-", ds="steps", color="blue",
												 label="link")
	plotfun(eCCDF2(dca)..., "-", ds="steps", color="green",
					label="cartels added")
	legend()
	ylabel("eCCDF")
	xlabel(xlab)
	PyPlot.grid()
	tight_layout()
end

Random.seed!(101)

rnetsrw, rnetsrl, rnetsca = simulatepubnets()

plotres(:degrees, "degrees")
savefig("../paperfigs/cartel_footprint-degree.pdf")

plotres(:strengthes, "strength")
savefig("../paperfigs/cartel_footprint-strength.pdf")

plotres(:weights, "weight", false, true, semilogy)
savefig("../paperfigs/cartel_footprint-weight.pdf")

plotres(:npapers, "number of papers")
savefig("../paperfigs/cartel_footprint-npapers.pdf")

plotres(:wpapers, "weigthed number of papers")
savefig("../paperfigs/cartel_footprint-wpapers.pdf")

plotres(:nauthors, "number of authors", false, true, semilogy)
savefig("../paperfigs/cartel_footprint-nauthors.pdf")


