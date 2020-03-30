
include("../analyse-publications/PaperPump.jl")
using RCall

function calctables(pn::PubNet, basename::String)
	R"source('../analyse-publications/fitpw.R')"
	tocalc = [:degrees, :nauthors, :npapers]
	for i in tocalc
		fname = string(basename,"-", String(i), ".txt")
		d = getproperty(pn, i)
		R"heavy.tailed($d, 100, 1000, $fname)"
	end
end

pn = PubNet("../analyse-publications/MTMT/MTMTpubmat.mat")
calctables(pn, "../paperfigs/MTMT-descr")

pn = PubNet("../analyse-publications/dblp/dblppubmat.mat")
calctables(pn, "../paperfigs/dblp-descr")
