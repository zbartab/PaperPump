#! /home/apa/bin/julia

# produce the publication and collaboration matrices for the dblp
# dataset

include("../CollaborationNetworks.jl")
include("ProcessDBLPrecords.jl")

function readargs(ARGV::Array{String,1})
	if length(ARGV) == 3
		niter = parse(Int, pop!(ARGV))
	else
		niter = 1000000
	end
	if length(ARGV) == 2
		nrewire = parse(Int, pop!(ARGV))
	else
		nrewire = 20
	end
	if length(ARGV) == 1
		filein = pop!(ARGV)
	else
		filein = "dblppubmat.csv"
	end
	return filein, nrewire, niter
end

function main(ARGV::Array{String,1})
	filein, nrewire, niter = readargs(ARGV)
	dblppubmat = read_spmatrix(filein)
	println("publication matrix `$(filein)` read")
	filein = replace(filein, r"\.csv$" => "")
	fileout = replace(filein, "pub" => "col")
	for i in 1:nrewire
		re_pubmat = rewire(dblppubmat, niter)
		write_spmatrix(string(filein, "-rewired-", i, ".csv"), re_pubmat)
		re_colmat = collaborationmatrix(re_pubmat)
		write_spmatrix(string(fileout, "-rewired-", i, ".csv"), re_colmat)
	end
	return nothing
end

main(ARGS)
