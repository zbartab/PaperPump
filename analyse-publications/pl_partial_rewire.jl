#! /usr/bin/julia -p 14

# partially rewire publication matrices
# partial rewiring means that only a given portion of the edges are
# reconnected

@everywhere begin
	include("../CollaborationNetworks.jl")
end

function readargs(ARGV::Array{String,1})
	length(ARGV) != 3 && error("Not enough arguments given!")
	npapers = parse(Int, pop!(ARGV))
	nrewire = parse(Int, pop!(ARGV))
	filein = pop!(ARGV)
	return filein, nrewire, npapers
end

function main(ARGV::Array{String,1})
	filein, nrewire, npapers = readargs(ARGV)
	pubmat = read_scimat(filein)
	println("publication matrix `$(filein)` read")
	filein = replace(filein, r"\....$" => "")
	fileout = replace(filein, "pub" => "col")
	cases = String[]
	for prw in [0.001, 0.005, 0.01, 0.05, 0.1, 0.25], i in 1:nrewire
		push!(cases, string(prw,",",i))
	end
	dorewire = function(s)
		prw, i = split(s, ",")
		prw = parse(Float64, prw)
		i = parse(Int, i)
		re_pubmat = rewire(pubmat, prw)
		write_scimat(string(filein, "-", prw, "-rewired-", i, ".mat"), re_pubmat)
		re_pubmat = selectauthors(re_pubmat, npapers)
		re_colmat = collaborationmatrix(re_pubmat)
		write_scimat(string(fileout, "-", prw, "-rewired-", i, ".mat"), re_colmat)
		return nothing
	end
	pmap(dorewire, cases)
	return nothing
end

main(ARGS)
