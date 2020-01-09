
# Handle MTMT records

using JSON

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
	papersIDs = Dict{UInt64, String}()
	for f in filelist
		(recs = read_MTMT(f)) == nothing && continue
		idstr = replace(f, r".*/(.*)_.*$" => s"\1")
		id = hash(idstr)
		authorIDs[id] = idstr
		!haskey(authorrecs, id) && (authorrecs[id] = Array{UInt64}[])
		for r in recs
			if hasranking(r) 
				p = r["mtid"]
				pid = hash(p)
				papersIDs[pid] = string(p)
				push!(authorrecs[id], pid)
			end
		end
		length(authorrecs[id]) == 0 && delete!(authorrecs, id)
	end
	return authorrecs, authorIDs, papersIDs
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
function recs2pubmat(aurecs::Dict{UInt64, Array{UInt64,1}},
										 auIDs::Dict{UInt64, String},
										 pIDs::Dict{UInt64, String})
	ps = collectpapers(aurecs)
	as = collectauthors(aurecs)
	pubmat = spzeros(length(ps), length(as))
	paperIDs = Array{String,1}(undef, length(ps))
	authorIDs = Array{String,1}(undef, length(as))
	for k in keys(aurecs)
		ai = as[k]
		authorIDs[ai] = auIDs[k]
		for p in aurecs[k]
			p_i = ps[p]
			pubmat[p_i, ai] = 1
			paperIDs[p_i] = pIDs[p]
		end
	end
	return PubMat(pubmat, authorIDs, paperIDs)
end

"""
    files2MTMTpubmat()

Create a publication matrix from the files containing MTMT records.
Return a dict with the matrix and its row and column names. `datadir`
gives the directory where the files reside.
"""
function files2MTMTpubmat(datadir::String)
	MTMT_dir = joinpath(pwd(), datadir);
	fMTMT = readdir(MTMT_dir);
	filter!((x) -> occursin(r"[0-9]\.json$", x), fMTMT);
	fMTMT = joinpath.(MTMT_dir, fMTMT);
	fMTMT = fMTMT[1:100]
	recs, auIDs, pIDs = recordauthors(fMTMT)
	return recs2pubmat(recs, auIDs, pIDs)
end
