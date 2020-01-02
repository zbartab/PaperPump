
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
