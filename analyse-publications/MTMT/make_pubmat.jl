#! /home/apa/bin/julia

# produce the publication and collaboration matrices for the MTMT
# dataset

include("../CollaborationNetworks.jl")
include("ProcessMTMTrecords.jl")

function main(ARGV::Array{String,1})
	if length(ARGV) == 1
		fileout = ARGV[1]
	else
		fileout = "MTMTpubmat.txt"
	end
	MTMT_dir = joinpath(pwd(), "MTMT-downloads");
	pm = files2MTMTpubmat(MTMT_dir)
	write_scimat(fileout, pm)
	println("publication matrix `$(fileout)` written")
	return nothing
end

main(ARGS)
