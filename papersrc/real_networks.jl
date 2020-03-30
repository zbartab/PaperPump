
# calculate some characteristics for real publication networks

include("../analyse-publications/PaperPump.jl")
matplotlib.use("PDF")

function partcopy(pn::PubNet)
	propnames = propertynames(pn)
	r = Dict()
	for n in propnames[4:11]
		r[n] = copy(getproperty(pn, n))
	end
	r[:communities] = map(length, values(pn.hierarchy[:lev1]))
	r[:cartels] = map(length, pn.cartels)
	return r
end

function plotres(x, s, xlab, plotfun=loglog, rmzero=true)
	y = x[:real][s]
	rmzero && (y = y[y .> 0])
	plotfun(eCCDF2(y)..., "-", ds="steps")
	for i in 1:20
		y = x[Symbol("rnd$(i)")][s]
		rmzero && (y = y[y .> 0])
		plotfun(eCCDF2(y)..., "--", ds="steps",
						color="orange", alpha=0.1)
	end
	legend(["real network"])
	xlabel(xlab)
	ylabel("eCCDF")
	tight_layout()
end

mtmt = Dict()
pn = PubNet("../analyse-publications/MTMT/MTMTpubmat.mat")
mtmt[:real] = partcopy(pn)
for i in 1:20
	pn = PubNet("../analyse-publications/MTMT/MTMTpubmat-rewired-$(i).mat")
	mtmt[Symbol("rnd$(i)")] = partcopy(pn)
end
pn = nothing

figure();
plotres(mtmt, :degrees, "degree");
savefig("../paperfigs/MTMT-degree.pdf")
figure();
plotres(mtmt, :strengthes, "strength");
savefig("../paperfigs/MTMT-strength.pdf")
figure();
plotres(mtmt, :weights, "weight", semilogy);
savefig("../paperfigs/MTMT-weight.pdf")
figure();
plotres(mtmt, :nauthors, "number of authors");
savefig("../paperfigs/MTMT-nauthors.pdf")
figure();
plotres(mtmt, :npapers, "number of papers");
savefig("../paperfigs/MTMT-npapers.pdf")
figure();
plotres(mtmt, :wpapers, "weighted number of papers");
savefig("../paperfigs/MTMT-wpapers.pdf")

mtmt = Dict()
pn = PubNet("../analyse-publications/dblp/dblppubmat.mat")
mtmt[:real] = partcopy(pn)
for i in 1:20
	pn = PubNet("../analyse-publications/dblp/dblppubmat-rewired-$(i).mat")
	mtmt[Symbol("rnd$(i)")] = partcopy(pn)
end
pn = nothing

figure();
plotres(mtmt, :degrees, "degree");
savefig("../paperfigs/dblp-degree.pdf")
figure();
plotres(mtmt, :strengthes, "strength");
savefig("../paperfigs/dblp-strength.pdf")
figure();
plotres(mtmt, :weights, "weight", semilogy);
savefig("../paperfigs/dblp-weight.pdf")
figure();
plotres(mtmt, :nauthors, "number of authors");
savefig("../paperfigs/dblp-nauthors.pdf")
figure();
plotres(mtmt, :npapers, "number of papers");
savefig("../paperfigs/dblp-npapers.pdf")
figure();
plotres(mtmt, :wpapers, "weighted number of papers");
savefig("../paperfigs/dblp-wpapers.pdf")

