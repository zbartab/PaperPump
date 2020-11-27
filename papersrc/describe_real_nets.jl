
using Printf

include("../analyse-publications/PaperPump.jl")
matplotlib.use("PDF")

function printnet(net, variable)
	x = getproperty(pn[net], variable)
	dx = [length(x), mean(x), quantile(x)...]
	print(f,net)
	for i in 1:length(dx)
		if i == 1
			@printf(f, " | %d", dx[i])
		elseif i == 2
			@printf(f, " | %0.3f", dx[i])
		elseif vars[variable][2]
			@printf(f, " | %d", dx[i])
		else
			@printf(f, " | %0.3f", dx[i])
		end
	end
	println(f)
end

#pn = Dict()
#pn["MTMT"] = PubNet("../analyse-publications/MTMT/MTMTpubmat.mat")
#pn["dblp"] = PubNet("../analyse-publications/dblp/dblppubmat.mat")

# variables

prop_names = propertynames(pnMTMT)[4:10]
vars = Dict()
vars[:nauthors] = ("authors", true)
vars[:npapers] = ("papers", true)
vars[:wpapers] = ("weighted papers", false)
vars[:strengthes] = ("strengths", false)
vars[:weights] = ("weights", false)
vars[:degrees] = ("degrees", true)
vars[:clustcoefs] = ("clustering coeff.", false)

f = open("../paperfigs/real_nets_description.txt", "w")

# table header

println(f, "Network | _n_ | _m_ | min | LQ | M | UQ | max ")
println(f, ":------:|:---:|:---:|:---:|:--:|:-:|:--:|:---:")

# table
for v in prop_names
	println(f, "_", vars[v][1], "_ | | | | | | | ")
	printnet("MTMT", v)
	printnet("dblp", v)
end

close(f)

loglog(eCCDF2(pn["MTMT"].degrees)..., "-", ds="steps", color="blue",
			 label="MTMT")
loglog(eCCDF2(pn["dblp"].degrees)..., "-", ds="steps", color="orange",
			 label="dblp")
legend()
xlabel("degree")
ylabel("eCCDF")
tight_layout()
savefig("../paperfigs/real_nets_description.pdf")
