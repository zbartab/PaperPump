
include("../analyse-publications/PaperPump.jl")

pn = PubNet("../analyse-publications/MTMT/MTMTpubmat.mat", 0.5)
samplingproductivityraw("../analyse-publications/MTMT/MTMTpubmat-productivity_raw.csv",
												pn, pn.cartels, 10000)
gs, gp, np, wp = groupsproductivity(pn, pn.cartels)
cartprod = DataFrame(groupsize = gs, groupprod = gp, npapers = np,
										 wpapers = wp)
CSV.write("../analyse-publications/MTMT/MTMTpubmat-cartel_productivity.csv",
					cartprod)

pn = PubNet("../analyse-publications/dblp/dblppubmat.mat", 0.5)
samplingproductivityraw("../analyse-publications/dblp/dblppubmat-productivity_raw.csv",
												pn, pn.cartels, 10000)
gs, gp, np, wp = groupsproductivity(pn, pn.cartels)
cartprod = DataFrame(groupsize = gs, groupprod = gp, npapers = np,
										 wpapers = wp)
CSV.write("../analyse-publications/dblp/dblppubmat-cartel_productivity.csv",
					cartprod)
