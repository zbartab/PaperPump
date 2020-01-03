
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

