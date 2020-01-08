
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
											 plottitle="")
	x = []
	y = []
	for k in keys(histdict)
		push!(x, k)
		push!(y, histdict[k])
	end
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
	plot(log10.(d["kn"]), log10.(d["pk"]), label=label)
	scatter(log10.(d["kn"]), log10.(d["pk"]), label=label)
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
function graphplot(g; cutoff=0.4, backend=PDF, filename="proba",
									 nodelabel=true, layout=random_layout)
	backend == PDF && (filename *= ".pdf")
	backend == PNG && (filename *= ".png")
	W = Weights(g)
	tocutoff = W .> cutoff
	edgewidth = ones(length(W))
	edgewidth[tocutoff] .= 5
	edgecolor = [colorant"lightgrey" for i in 1:length(W)]
	edgecolor[tocutoff] .= colorant"black"
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
	if sum(tocutoff) == 0
		maxlinewidth = 0.2
	else
		maxlinewidth = 2
	end
	Compose.draw(backend(filename, 50cm, 50cm),
							 gplot(g, layout=layout, edgelinewidth=edgewidth,
										 EDGELINEWIDTH=maxlinewidth, nodesize=nodesize,
										 nodelabel=nodelabs,
										 nodefillc=nodefillc, edgestrokec=edgecolor))
end

