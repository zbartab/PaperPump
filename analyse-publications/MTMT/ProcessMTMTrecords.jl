
# Handle MTMT records

"""
    read_MTMT(file)

Read the MTMT records of an author from JSON file downloaded from
mtmt.hu.
"""
function read_MTMT(file)
	local r
	try
		r = JSON.parsefile(file)
	catch
		return nothing
	end
	return r["content"]
end

"""
    recordauthors(filelist)

Record the papers written by authors stored in MTMT json files in
`filelist`.
"""
function recordauthors(filelist::Array{String,1})
	authorrecs = Dict{UInt64, Array{UInt64,1}}()
	authorIDs = Dict{UInt64, String}()
	for f in filelist
		(recs = read_MTMT(f)) == nothing && continue
		idstr = replace(f, r".*/(.*)_.*$" => s"\1")
		id = hash(idstr)
		authorIDs[id] = idstr
		!haskey(authorrecs, id) && (authorrecs[id] = Array{UInt64}[])
		for r in recs
			hasranking(r) && push!(authorrecs[id], hash(r["mtid"]))
		end
		length(authorrecs[id]) == 0 && delete!(authorrecs, id)
	end
	return authorrecs, authorIDs
end

"""
    hasranking(cikk)

Returns `true` if MTMT paper record `cikk` has journal ranking field.
"""
function hasranking(cikk)
	return haskey(cikk, "ratings") && haskey(cikk["ratings"][1], "ranking")
end

"""
    collectpapers(aurecs)

Return a dict of unique paper IDs.
"""
function collectpapers(aurecs::Dict{UInt64, Array{UInt64,1}})
	papers = Dict{UInt64,Int}()
	for ps in values(aurecs)
		length(ps) == 0 && continue
		for p in ps
			papers[p] = 1
		end
	end
	i = 1
	for k in keys(papers)
		papers[k] = i
		i += 1
	end
	return papers
end

"""
    collectauthors(aurecs)

Return a dict of unique authors IDs.
"""
function collectauthors(aurecs::Dict{UInt64, Array{UInt64,1}})
	authors = Dict{UInt64,Int}()
	i = 1
	for k in keys(aurecs)
		authors[k] = i
		i += 1
	end
	return authors
end

"""
    recs2pubmat(aurecs)

Create a publication matrix from the author records. It returns (as a
tuple) the publication matrix and a hash table to allow the
identification of authors in the publication matrix.
"""
function recs2pubmat(aurecs::Dict{UInt64, Array{UInt64,1}})
	ps = collectpapers(aurecs)
	as = collectauthors(aurecs)
	pubmat = spzeros(length(ps), length(as))
	for k in keys(aurecs)
		ai = as[k]
		for p in aurecs[k]
			pubmat[ps[p], ai] = 1
		end
	end
	return pubmat, as
end
