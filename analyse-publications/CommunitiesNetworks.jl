
using LightGraphs, MetaGraphs

# community detection algorithm according to @lancichinetti2009


"""
    findcommunities(graph, seed, α = 1.5, β = 1.5)

Returns a cover of `graph`, i.e. a division of the graph into
communities.
"""
function findcommunities(graph; α = 1.5, β = 1.5)
	#println("# 1")
	uncovered = collect(1:nv(graph))
	cover = Array{Array{Int, 1}, 1}()
	commfits = Array{Float64, 1}()
	while length(uncovered) > 0
		seed = rand(uncovered)
		comm, fit = findcommunity(graph, seed, α = α, β = β)
		#comm, f = [1,2], 0.1234
		push!(commfits, fit)
		push!(cover, comm)
		setdiff!(uncovered, comm)
	end
	return cover, commfits
end

"""
    groupneighbours(graph, group)

Return the neighbours of a group, e.g. nodes that are connected to the
nodes in the group, but they are not members of the group.
"""
function groupneighbours(graph, group::Array{Int, 1})
	neighbours = Set{Int}()
	for m in group
		ns = neighbors(graph, m)
		if length(ns) > 0
			for n in ns
				push!(neighbours, n)
			end
		end
	end
	return collect(setdiff(neighbours, group))
end

"""
    groupfitness(graph, group)

Return the fitness of the group calculated according to @wachs2019a
"""
function groupfitness(graph, group::Array{Int, 1},
											neighbours::Array{Int, 1}; α = 1.0, β = 0.0)
	#neighbours = groupneighbours(graph, group)
	strengthout = 0.0
	for n in neighbours
		for m in group
			has_edge(graph, m, n) && (strengthout += LightGraphs.weights(graph)[m, n])
		end
	end
	strengthin = 0.0
	lg = length(group)
	for i in 1:(lg-1)
		for j in (i+1):lg
			has_edge(graph, group[i], group[j]) && 
			(strengthin += LightGraphs.weights(graph)[group[i], group[j]])
		end
	end
	#println(">>> ", strengthin, "\t", strengthout)
	return strengthin/((strengthin + strengthout)^α *
										 length(group)^β)
end

"""
    findcommunity(graph, seed, alpha, beta)

Return a community find by the FLK algorithm [@wachs2019a]. Parameter
`graph` ia the network in which to find community, `seed` is the
starting node, while `alpha` and `beta` are free parameters influencing grouping.
"""
function findcommunity(graph, seed::Int; α=1.0, β=0.0)
	cneighbours = neighbors(graph, seed)
	length(cneighbours) == 0 && return [seed], 0.0
	community = [seed]
	cfitness = groupfitness(graph, community, cneighbours, α=α, β=β)
	while true
		#println("### ", community)
		maxfit = -1.0e10
		newcfit = Array{Float64, 1}()
		newmember = Array{Int, 1}()
		newcn = Array{Array{Int, 1}, 1}()
		for n in cneighbours
			comm = vcat(community, n)
			cn = groupneighbours(graph, comm)
			cf = groupfitness(graph, comm, cn, α=α, β=β)
			#println(n, "\t", round(cf, digits=5), "\t",
							#round(cfitness, digits=5), "\t",
							#round(cf-cfitness, digits=5))
			#println(cf)
			fit = cf -cfitness
			if fit > maxfit
				maxfit = fit
				newcfit = [cf]
				newmember = [n]
				newcn = [cn]
			elseif fit == maxfit
				push!(newcfit, cf)
				push!(newmember, n)
				push!(newcn, cn)
			end
		end
		if maxfit < 0.0
			return (community, cfitness)
		else
			i = rand(1:length(newcfit))
			community = vcat(community, newmember[i])
			cfitness = newcfit[i]
			cneighbours = newcn[i]
			#println(cfitness, "\t", length(community), "\t", newmember[i])
		end
	end
end

function plgraph(g, comm)
	# to run this you need GraphRecipes and Plots
	cols = ones(Int, nv(g))
	cols[comm] .= 2
	graphplot(g, curves=false, names=1:nv(g), nodesize=0.2,
						method=:spring, nodecolor=cols)
end
