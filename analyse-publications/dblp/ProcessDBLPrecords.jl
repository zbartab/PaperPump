
# code to process the dblp records and produce collection of papers

using SparseArrays

"""
    read_dblprecords(file)

Read the records of the dblp database from `file` into an array of strings.
"""
function read_dblprecords(file)
	f = open(file)
	lines = String[]
	for line in eachline(f)
		push!(lines, line)
	end
	close(f)
	return lines
end

"""
    recordpapers(startyear, stopyear)

Record the authors who write the papers in the dblp database, passed as
`lines`. Only papers published between `startyear` and `stopyear`
(inclusive) considered.
"""
function recordpapers(lines::Array{String,1}, startyear::Int,
											stopyear::Int)
	n_records = div(length(lines), 3)
	dblp_dict = Dict{UInt64, Dict{String, Union{Int, Array{UInt64,1}}}}()
	for i in 1:n_records
		j = i * 3
		year = parse(Int, lines[j])
		startyear <= year <= stopyear || continue
		doi = hash(lines[j-2])
		dblp_dict[doi] = Dict("authors" => hash.(split(lines[j-1], "; ")),
													"year" => year)
	end
	return dblp_dict
end

"""
    collectauthors(aurecs)

Return a dict of unique paper IDs.
"""
function collectauthors(aurecs::Dict{UInt64,
																		Dict{String,Union{Int64,
																										 Array{UInt64,1}}}})
	authors = Dict{UInt64,Int}()
	for ps in values(aurecs)
		for p in ps["authors"]
			authors[p] = 1
		end
	end
	i = 1
	for k in keys(authors)
		authors[k] = i
		i += 1
	end
	return authors
end

"""
    collectpapers(aurecs)

Return a dict of unique papers IDs.
"""
function collectpapers(aurecs::Dict{UInt64,
																		 Dict{String,Union{Int64,
																											 Array{UInt64,1}}}})
	papers = Dict{UInt64,Int}()
	i = 1
	for k in keys(aurecs)
		papers[k] = i
		i += 1
	end
	return papers
end

"""
    recs2pubmat(aurecs)

Create a publication matrix from the author records (i..e. it works on
the result of `recordpapers`). It returns (as a tuple) the publication
matrix and a hash table to allow the identification of authors in the
publication matrix.
"""
function recs2pubmat(aurecs::Dict{UInt64,
																	Dict{String,Union{Int64,
																										Array{UInt64,1}}}})
	ps = collectpapers(aurecs)
	as = collectauthors(aurecs)
	pubmat = spzeros(length(ps), length(as))
	for k in keys(aurecs)
		p_i = ps[k]
		for a in aurecs[k]["authors"]
			pubmat[p_i, as[a]] = 1
		end
	end
	return pubmat, as
end
