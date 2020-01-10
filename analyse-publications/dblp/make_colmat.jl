#! /home/apa/bin/julia

# produce the publication and collaboration matrices for the dblp
# dataset

include("../CollaborationNetworks.jl")
include("ProcessDBLPrecords.jl")

function main(ARGV::Array{String,})
	filein = "dblppubmat.txt"
	npapers = -Inf
	if length(ARGV) >= 2
		filein = ARGV[1]
		npapers = parse(Float64, ARGV[2])
	elseif length(ARGV) == 1
		filein = ARGV[1]
	end
	dblppubmat = read_scimat(filein)
	println("publication matrix `$(filein)` read")
	dblppubmat = selectauthors(dblppubmat, npapers)
	fileout = replace(filein, "pub" => "col")
	dblpcolmat = collaborationmatrix(dblppubmat);
	write_scimat(fileout, dblpcolmat)
	println("collaboration matrix `$(fileout)` written")
	return nothing
end

main(ARGS)
