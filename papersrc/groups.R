
# produce a figure for the author groups

library(igraph)

create.network <- function(gr.name, gr.size, n.coauthors, gr.connected=FALSE) {
	#group <- paste(gr.name, LETTERS[1:gr.size], sep="")
	group <- paste(gr.name, 1:gr.size, sep="")
	coauthors <- sapply(group,
											function(a) paste(a, letters[1:n.coauthors], sep=""))
	coauthors <- as.vector(coauthors)
	edges <- data.frame(authors=rep(group, rep(n.coauthors, gr.size)),
											coauthors=coauthors)
	edges <- as.matrix(edges)
	n.edges <- nrow(edges)
	if(gr.connected) {
		i <- length(group)
		p <- (i*(i-1))/2
		ee <- matrix("", ncol=2, nrow=p)
		k <- 1
		for(i in 1:(gr.size-1)) {
			for(j in (i+1):gr.size) {
				ee[k, 1] <- group[i]
				ee[k, 2] <- group[j]
				k <- k+1
			}
		}
		edges <- rbind(edges, ee)
	}
	g <- graph_from_edgelist(edges, directed=FALSE)
	V(g)$size <- 30
	V(g)$color[grepl("[a-c]$", names(V(g)))] <- "salmon"
	V(g)$color[grepl("[^a-c]$", names(V(g)))] <- "lightgreen"
	E(g)$color <- "black"
	E(g)$width <- 4
	E(g)$width[1:n.edges] <- 1
	g
}

G.A <- G.B <- 4 # group size
c.A <- c.B <- 3 # number of coauthors from outside
network.A <- create.network("A", G.A, c.A)
network.B <- create.network("B", G.B, c.B, TRUE)
coords <- layout_(network.B, nicely())
#png(file="groups.png", width=700, height=400)
pdf(file="paperfigs/groups.pdf", width=10, height=6)
op <- par(mfcol=c(1,2))
plot(network.A, layout=coords, main="group A",
		 mark.groups=list(V(network.A)[nchar(names(V(network.A))) == 2]))
plot(network.B, layout=coords, main="group B",
		 mark.groups=list(V(network.B)[nchar(names(V(network.B))) == 2]))
par(op)
dev.off()
