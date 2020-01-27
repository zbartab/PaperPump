
# script to visualise a sample publication network

include("../analyse-publications/CollaborationNetworks.jl")
include("../analyse-publications/PlotCollaborationNetworks.jl")
include("../analyse-publications/StatCollaborationNetworks.jl")
include("../analyse-publications/simulations/RandomPublicationNetworks.jl")

mat = [0 1 0 0 0;
			 0 1 1 0 0;
			 1 0 0 1 0;
			 0 1 1 0 1;
			 0 0 1 0 1;
			 1 1 0 0 0]
mat = float.(mat)
mat = sparse(mat)
pm = generate_publicationmatrix(mat)

com = collaborationmatrix(pm)
cog = collaborationgraph(com)
#lx, ly = spring_layout(cog, C=20)
lx = [1.0, -0.4292, -0.1319, -1.0, 0.95489]
ly = [-0.6075, 1.0, -1.0, -0.0386, 0.66608]
graphplot(cog, lx, ly, cutoff=2, width=15cm, height=15cm, maxlinewidth=3,
					filename="sample_graph")

addcartel!(pm, [3,4], 1.0)
com = collaborationmatrix(pm)
coga = collaborationgraph(com)

graphplot(cog, lx, ly, cutoff=0.9, width=15cm, height=15cm, maxlinewidth=3,
					filename="sample_graph_cartel")

W = Weights(cog)
Wc = Weights(coga)

#semilogy(eCCDF(W)..., ".-", label="no cartel")
#semilogy(eCCDF(Wc)..., ".-", label="cartel added")
#legend()
#xlim(-0.05, 1.05)
#xlabel("weights")
#ylabel("CCDF")
