
# script to visualise a sample publication network

include("../analyse-publications/PaperPump.jl")

Random.seed!(101)

mat = [0 1 0 0 0 0 0;
			 0 0 0 0 0 1 1;
			 0 0 0 0 1 0 1;
			 0 1 1 0 0 0 0;
			 1 0 0 1 0 1 0;
			 0 1 1 0 0 0 0;
			 0 0 1 0 1 0 0;
			 1 1 0 0 0 0 1]
mat = sparse(mat)
pm = generate_publicationmatrix(mat)

com = collaborationmatrix(pm)
cog = collaborationgraph(com)
lx, ly = spring_layout(cog, C=20)
#lx = [1.0, -0.4292, -0.1319, -1.0, 0.95489]
#ly = [-0.6075, 1.0, -1.0, -0.0386, 0.66608]
graphplot(cog, lx, ly, cutoff=2, width=15cm, height=15cm, maxlinewidth=2,
					filename="../paperfigs/sample_graph")

pma = deepcopy(pm)
addcartel!(pma, [1,4,6], 1.0)
coma = collaborationmatrix(pma)
coga = collaborationgraph(coma)
graphplot(coga, lx, ly, cutoff=0.9, width=15cm, height=15cm, maxlinewidth=2,
					filename="../paperfigs/sample_graph_cartel")

W = Weights(cog)
Wc = Weights(coga)

semilogy(eCCDF2(W)..., "-", label="no cartel", ds="steps")
semilogy(eCCDF2(Wc)..., "--", label="cartel added", ds="steps")
legend()
xlim(-0.05, 1.05)
xlabel("weights")
ylabel("ECCDF")
tight_layout()
savefig("../paperfigs/sample_graph_weights.pdf")

function import_pubmat(pm::PubMat, file::String="../paperfigs/sample_graph.mat")
	f = open(file, "w")
	m = Matrix(pm.mat)
	for i in 1:size(m,2)
		print(f, " & \\node {$(findfirst((x) -> x == i, pm.authorIDs))};")
	end
	println(f, " \\\\ ")
	for i in 1:size(m,1)
		print(f, "\\node {$(findfirst((x) -> x == i, pm.paperIDs))};")
		for j in 1:size(m,2)
			print(f, " & ")
			if m[i,j] == 0
				print(f, "\\node {.};")
			else
				print(f, "\\node {1};")
			end
		end
		println(f, " \\\\ ")
	end
	close(f)
end

import_pubmat(pm)
import_pubmat(pma, "../paperfigs/sample_graph_cartel.mat")


function import_colmat(pm::ColMat,
											 file::String="../paperfigs/sample_graph-col.mat")
	f = open(file, "w")
	m = Matrix(pm.mat)
	for i in 1:size(m,2)
		print(f, " & \\node {$(findfirst((x) -> x == i, pm.authorIDs))};")
	end
	println(f, " \\\\ ")
	for i in 1:size(m,1)
		print(f, "\\node {$(findfirst((x) -> x == i, pm.authorIDs))};")
		for j in 1:size(m,2)
			print(f, " & ")
			if m[i,j] == 0
				print(f, "\\node {.};")
			else
				print(f, "\\node {$(round(m[i,j]*100)/100)};")
			end
		end
		println(f, " \\\\ ")
	end
	close(f)
end

import_colmat(com)
import_colmat(coma, "../paperfigs/sample_graph_cartel-col.mat")
