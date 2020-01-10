#! /home/apa/bin/julia -p 6

# produce the publication and collaboration matrices for the dblp
# dataset

@everywhere begin
	include("../CollaborationNetworks.jl")
	include("ProcessMTMTrecords.jl")
end

function readargs(ARGV::Array{String,1})
	if length(ARGV) == 4
		npapers = parse(Int, pop!(ARGV))
	else
		npapers = 1000000
	end
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
		filein = "MTMTpubmat.txt"
	end
	return filein, nrewire, niter, npapers
end

function main(ARGV::Array{String,1})
	filein, nrewire, niter, npapers = readargs(ARGV)
	MTMTpubmat = read_scimat(filein)
	println("publication matrix `$(filein)` read")
	filein = replace(filein, r"\....$" => "")
	fileout = replace(filein, "pub" => "col")
	dorewire = function(i)
		re_pubmat = rewire(MTMTpubmat, niter)
		write_scimat(string(filein, "-rewired-", i, ".txt"), re_pubmat)
		re_pubmat = selectauthors(re_pubmat, npapers)
		re_colmat = collaborationmatrix(re_pubmat)
		write_scimat(string(fileout, "-rewired-", i, ".txt"), re_colmat)
		return nothing
	end
	pmap(dorewire, 1:nrewire)
	return nothing
end

main(ARGS)
