#! /home/apa/bin/julia

# produce the publication and collaboration matrices for the dblp
# dataset

include("../CollaborationNetworks.jl")
include("ProcessDBLPrecords.jl")

function main(ARGV::Array{String,1})
	starty = parse(Int, ARGV[1])
	stopy = parse(Int, ARGV[2])
	dblp_file = joinpath(pwd(), "dblp-data", "dblp-articles.txt")
	lines = read_dblprecords(dblp_file)
	parecsdblp = recordpapers(lines, starty, stopy);
	dblppubmat, dblpauids = recs2pubmat(parecsdblp);
	write_spmatrix("dblppubmat.csv", dblppubmat)
	println("publication matrix written")
	return nothing
end

main(ARGS)
