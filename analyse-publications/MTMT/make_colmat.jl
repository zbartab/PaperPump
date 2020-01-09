#! /home/apa/bin/julia

# produce the publication and collaboration matrices for the dblp
# dataset

include("../CollaborationNetworks.jl")
include("ProcessMTMTrecords.jl")

function main(ARGV::Array{String,})
	filein = "MTMTpubmat.txt"
	npapers = -Inf
	if length(ARGV) >= 2
		filein = ARGV[1]
		npapers = parse(Float64, ARGV[2])
	elseif length(ARGV) == 1
		filein = ARGV[1]
	end
	MTMTpubmat = read_scimat(filein)
	println("publication matrix `$(filein)` read")
	MTMTpubmat = selectauthors(MTMTpubmat, npapers)
	fileout = replace(filein, "pub" => "col")
	MTMTcolmat = collaborationmatrix(MTMTpubmat);
	write_scimat(fileout, MTMTcolmat)
	println("collaboration matrix `$(fileout)` written")
	return nothing
end

main(ARGS)
