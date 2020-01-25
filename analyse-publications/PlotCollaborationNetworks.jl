
# functions to plot the analyses of collaboration networks


using Colors
using PyPlot, Cairo, Compose
using GraphPlot

## the functions

"""
    plothistogram(histdict)

Draw a bar plot on the result of `histogram`.
"""
function plothistogram(histdict::Dict{Any,Any}; xlab="", ylab="",
											 plottitle="", relfreq=false)
	x = []
	y = []
	for k in keys(histdict)
		push!(x, k)
		push!(y, histdict[k])
	end
	relfreq && (y = y ./ sum(y))
	bar(x, y)
	title(plottitle)
	xlabel(xlab)
	ylabel(ylab)
	PyPlot.grid(true)
	tight_layout()
end

"""
    plotloglog(x, plottitle)

Plot the distribution of `x` on a log-log scale.
"""
function plotloglog(x, label="", plottitle="")
	d = loglogbins(x)
	if label == ""
		plot(log10.(d["kn"]), log10.(d["pk"]))
	else
		plot(log10.(d["kn"]), log10.(d["pk"]), label=label)
		legend()
	end
	scatter(log10.(d["kn"]), log10.(d["pk"]))
	title(plottitle)
	xlabel(L"$\log_{10}(k)$")
	ylabel(L"$\log_{10}(p_k)$")
	PyPlot.grid(true)
	tight_layout()
end

"""
    graphplot(g)

Produce a pdf plot of graph `g`
"""
function graphplot(g, loc_x, loc_y; cutoff=0.4, backend=PDF, filename="proba",
									 nodelabel=true)
	backend == PDF && (filename *= ".pdf")
	backend == PNG && (filename *= ".png")
	nodesize = degree(g)
	if nodelabel
		nodelabs = labels(g)
	else
		nodelabs = ""
	end
	#nodefillc = distinguishable_colors(nv(g), colorant"blue")
	#nodefillc = range(colorant"lightsalmon", stop=colorant"darksalmon",
										#length=nv(g))
	nodefillc = colorant"turquoise"
	if ne(g) == 0
		Compose.draw(backend(filename, 50cm, 50cm),
								 gplot(g, loc_x, loc_y, nodelabel=nodelabs,
											 nodefillc=nodefillc))
	else
		W = Weights(g)
		tocutoff = W .> cutoff
		edgewidth = ones(length(W))
		edgewidth[tocutoff] .= 5
		edgecolor = [colorant"lightgrey" for i in 1:length(W)]
		edgecolor[tocutoff] .= colorant"darkgrey"
		if sum(tocutoff) == 0
			maxlinewidth = 0.4
		else
			maxlinewidth = 2
		end
		Compose.draw(backend(filename, 50cm, 50cm),
								 gplot(g, loc_x, loc_y, edgelinewidth=edgewidth,
											 EDGELINEWIDTH=maxlinewidth, nodesize=nodesize,
											 nodelabel=nodelabs,
											 nodefillc=nodefillc, edgestrokec=edgecolor))
	end
	return nothing
end
function graphplot(g; layout::Function=random_layout, keyargs...)
	graphplot(g, layout(g)...; keyargs...)
	return nothing
end

"""
    donothing(x)

Do nothing, just return `x`. It acts as a placeholder.
"""
function donothing(x)
	return x
end

"""
    plotcartelseffects(psnc, pswc, measure, mytitle)

Plot how the authors' measure changes as a result of adding a cartel to
the pubkication network.
"""
function plotcartelseffects(psnc, pswc, measure; f=donothing, mytitle="")
	dx = 0.3
	boxplot(f(psnc[measure])[cartel], patch_artist=true, positions=[1],
					boxprops=Dict("facecolor" => "blue"))
	boxplot(f(pswc[measure])[cartel], patch_artist=true, positions=[2],
					boxprops=Dict("facecolor" => "lightblue"))
	boxplot(f(psnc[measure])[noncartel], patch_artist=true, positions=[3],
					boxprops=Dict("facecolor" => "red"))
	boxplot(f(pswc[measure])[noncartel], patch_artist=true, positions=[4],
					boxprops=Dict("facecolor" => "orange"))
	for i in cartel
		plot([1.0+dx,2.0-dx], [f(psnc[measure])[i], f(pswc[measure])[i]])
	end
	for i in noncartel
		plot([3.0+dx,4.0-dx], [f(psnc[measure])[i], f(pswc[measure])[i]])
	end
	title(mytitle)
	tight_layout()
end
