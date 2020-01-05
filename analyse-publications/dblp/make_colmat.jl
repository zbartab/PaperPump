#! /home/apa/bin/julia

# produce the publication and collaboration matrices for the dblp
# dataset

include("../CollaborationNetworks.jl")
include("ProcessDBLPrecords.jl")

function main(filein::String)
	filein == "" && (filein = "dblppubmat.csv")
	dblppubmat = read_spmatrix(filein)
	println("publication matrix `$(filein)` read")
	fileout = replace(filein, "pub" => "col")
	dblpcolmat = collaborationmatrix(dblppubmat);
	write_spmatrix("dblpcolmat.csv", dblpcolmat)
	println("collaboration matrix `$(fileout)` written")
	return nothing
end

main(ARGS[1])
