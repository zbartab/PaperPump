
# script to analyse the simulation

include("../analyse-publications/CollaborationNetworks.jl")
include("../analyse-publications/PlotCollaborationNetworks.jl")
include("../analyse-publications/StatCollaborationNetworks.jl")
include("../analyse-publications/simulations/RandomPublicationNetworks.jl")

# Elso eset
## mindenki egyformán teljesít
## nincs kollaboráció

k = repeat([3], 30)
#k = collect(1:30) .^ 2
#k = collect(1:30)
cartel = collect([1,2,3,29,30])
noncartel = setdiff(1:length(k), cartel)
pnnc, pnwc = generate_rndnet(k, 100000, cartel, 1.0,
														 filename="sa_eq-prod_nocoauth")
psnc =pubnetstats(pnnc)
pswc =pubnetstats(pnwc)

figure(figsize=(8, 3))
xlegend = ("no cartel", "cartel added")
subplot(1,2,1)
for i in cartel
	plot(xlegend, [psnc[:npapers][i], pswc[:npapers][i]])
end
xlim(-0.5, 1.5)
title("cartel members")
ylabel("number of papers")
subplot(1,2,2)
for i in noncartel
	plot(xlegend, [psnc[:npapers][i], pswc[:npapers][i]])
end
#xlim(0.75, 2.25)
xlim(-0.5, 1.5)
title("non-members")
ylabel("number of papers")
tight_layout()

savefig("sa_eq-prod_npapers.pdf")

figure(figsize=(8, 3))
xlegend = ("no cartel", "cartel added")
subplot(1,2,1)
for i in cartel
	plot(xlegend, [psnc[:wpapers][i], pswc[:wpapers][i]])
end
xlim(-0.5, 1.5)
title("cartel members")
ylabel("weighted number of papers")
subplot(1,2,2)
for i in noncartel
	plot(xlegend, [psnc[:wpapers][i], pswc[:wpapers][i]])
end
#xlim(0.75, 2.25)
xlim(-0.5, 1.5)
title("non-members")
ylabel("weighted number of papers")
tight_layout()

savefig("sa_eq-prod_wpapers.pdf")


# Masodik eset
## mindenki egyformán teljesít
## van kollaboráció

Random.seed!(101)
k = repeat([3], 30)
#k = collect(1:30) .^ 2
#k = collect(1:30)
cartel = collect([1,2,3,29,30])
noncartel = setdiff(1:length(k), cartel)
#pnnc, pnwc = generate_rndnet(k, 100, cartel, 1.0,
pnnc, pnwc = generate_rndnet(k, 60, cartel, 1.0,
														 filename="sa_eq-prod_coauth")
psnc =pubnetstats(pnnc)
pswc =pubnetstats(pnwc)

figure(figsize=(8, 3))
xlegend = ("no cartel", "cartel added")
subplot(1,2,1)
for i in cartel
	plot(xlegend, [psnc[:npapers][i], pswc[:npapers][i]])
end
xlim(-0.5, 1.5)
title("cartel members")
ylabel("number of papers")
subplot(1,2,2)
for i in noncartel
	plot(xlegend, [psnc[:npapers][i], pswc[:npapers][i]])
end
#xlim(0.75, 2.25)
xlim(-0.5, 1.5)
title("non-members")
ylabel("number of papers")
tight_layout()

savefig("sa_eq-prod_coauth_npapers.pdf")

figure(figsize=(8, 3))
xlegend = ("no cartel", "cartel added")
subplot(1,2,1)
for i in cartel
	plot(xlegend, [psnc[:wpapers][i], pswc[:wpapers][i]])
end
xlim(-0.5, 1.5)
title("cartel members")
ylabel("weighted number of papers")
subplot(1,2,2)
for i in noncartel
	plot(xlegend, [psnc[:wpapers][i], pswc[:wpapers][i]])
end
#xlim(0.75, 2.25)
xlim(-0.5, 1.5)
title("non-members")
ylabel("weighted number of papers")
tight_layout()

savefig("sa_eq-prod_coauth_wpapers.pdf")

# Harmadik eset
## szerzők teljesítménye különbözik
## van kollaboráció

Random.seed!(101)
#k = repeat([3], 30)
#k = collect(1:30) .^ 2
k = collect(1:30)
cartel = collect([1,2,3,29,30])
noncartel = setdiff(1:length(k), cartel)
pnnc, pnwc = generate_rndnet(k, 500, cartel, 1.0,
														 filename="sa_diff-prod_coauth")
psnc =pubnetstats(pnnc)
pswc =pubnetstats(pnwc)

figure(figsize=(8, 3))
xlegend = ("no cartel", "cartel added")
subplot(1,2,1)
for i in cartel
	plot(xlegend, [psnc[:npapers][i], pswc[:npapers][i]])
end
xlim(-0.5, 1.5)
title("cartel members")
ylabel("number of papers")
subplot(1,2,2)
for i in noncartel
	plot(xlegend, [psnc[:npapers][i], pswc[:npapers][i]])
end
#xlim(0.75, 2.25)
xlim(-0.5, 1.5)
title("non-members")
ylabel("number of papers")
tight_layout()

savefig("sa_diff-prod_coauth_npapers.pdf")

figure(figsize=(8, 3))
xlegend = ("no cartel", "cartel added")
subplot(1,2,1)
for i in cartel
	plot(xlegend, [psnc[:wpapers][i], pswc[:wpapers][i]])
end
xlim(-0.5, 1.5)
title("cartel members")
ylabel("weighted number of papers")
subplot(1,2,2)
for i in noncartel
	plot(xlegend, [psnc[:wpapers][i], pswc[:wpapers][i]])
end
#xlim(0.75, 2.25)
xlim(-0.5, 1.5)
title("non-members")
ylabel("weighted number of papers")
tight_layout()

savefig("sa_diff-prod_coauth_wpapers.pdf")
