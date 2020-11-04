
# load libraries

library(lattice)

# functions

## setting the data to plot

plot.prod <- function(m.raw, m.cart, measure, ylab, do.plot=TRUE) {
	pd = list()
	for(n in sort(unique(m.cart$groupsize))) {
		m <- sub("prod", "", measure)
		m <- sub("s$", "", m)
		pd[[paste("r", n, sep="")]] = m.raw[m.raw$measure == m,
																				paste("X", n, sep="")]
		pd[[paste("c", n, sep="")]] = m.cart[m.cart$groupsize == n, measure]
	}
	if(do.plot) {
		l <- unlist(pd)
		s <- sub("^([rc])([1-8]0*).*", "\\1\\2", names(l))

		a <- 1:(1.5*length(pd))
		a <- a[a %% 3 != 0]

		plot(1,1, xlim=range(a), ylim=range(l), log="y", type="n", xlab="group size",
				 ylab=ylab, xaxt="n")
		axis(side=1, at=rowMeans(matrix(a, ncol=2, byrow=TRUE)),
				 labels=sort(unique(m.cart$groupsize)))
		for(i in 1:length(a)) {
			if(a[i] %% 3 == 1) {
				#mcol <- rgb(0.5, 0, 0, 0.05)
				mcol <- rgb(1, 0.6, 0.2, 0.05)
				#bcol <- rgb(1, 0.6, 0.2, 1)
				bcol <- "black"
			} else {
				mcol <- rgb(0, 0.5, 0, 0.25)
				bcol <- "black"
			}
			points(jitter(rep(a[i], length(pd[[i]])), amount=0.25), pd[[i]], pch=16, 
						 col=mcol)
			boxplot(pd[[i]], outline=FALSE, log="y", add=TRUE, at=a[i], col=mcol,
							border=bcol, yaxt="n")
		}
	}
	return(invisible(pd))
}

cartel.z.scores <- function(p) {
	l <- length(p)
	ns <- as.numeric(unique(sub("^[rc]", "", names(p))))
	m.z <- list()
	n.c <- list()
	for(n in ns) {
		mr <- p[[paste("r", n, sep="")]]
		mc <- p[[paste("c", n, sep="")]]
		tot <- c(mr, mc)
		gr <- rep(c("r", "c"), c(length(mr), length(mc)))
		stot <- scale(tot)
		m.z[[paste("c", n, sep="")]] <- stot[gr=="c"]
		n.c[[paste("c", n, sep="")]] <- rep(n, length(mc))
	}
	unlist(m.z)
}

cartel.z.diff <- function(p) {
	l <- length(p)
	ns <- as.numeric(unique(sub("^[rc]", "", names(p))))
	m.z <- numeric(length(ns))
	n.c <- numeric(length(ns))
	i <- 1
	for(n in ns) {
		mr <- p[[paste("r", n, sep="")]]
		mc <- p[[paste("c", n, sep="")]]
		tot <- c(mr, mc)
		gr <- rep(c("r", "c"), c(length(mr), length(mc)))
		stot <- scale(tot)
		m.z[i] <- mean(stot[gr=="c"])
		n.c[i] <- length(mc)
		i <- i+1
	}
	weighted.mean(m.z, n.c)
}

cartel.diff <- function(p, Q=0.5) {
	l <- length(p)
	ns <- as.numeric(unique(sub("^[rc]", "", names(p))))
	m.z <- numeric(length(ns))
	n.c <- numeric(length(ns))
	i <- 1
	for(n in ns) {
		mr <- p[[paste("r", n, sep="")]]
		mc <- p[[paste("c", n, sep="")]]
		m.z[i] <- sum(mc < quantile(mr, Q))
		n.c[i] <- length(mc)
		i <- i+1
	}
	sum(m.z)/sum(n.c)
}

cartel.median.diff <- function(p) {
	l <- length(p)
	ns <- as.numeric(unique(sub("^[rc]", "", names(p))))
	m.z <- list()
	n.c <- numeric(length(ns))
	i <- 1
	for(n in ns) {
		mr <- p[[paste("r", n, sep="")]]
		mc <- p[[paste("c", n, sep="")]]
		m.z[[paste("c", n, sep="")]] <- (mc-median(mr))/median(mr)
		n.c[i] <- length(mc)
		i <- i+1
	}
	median(unlist(m.z))
}

cartel.m.scores <- function(p) {
	l <- length(p)
	ns <- as.numeric(unique(sub("^[rc]", "", names(p))))
	m.z <- list()
	for(n in ns) {
		mr <- p[[paste("r", n, sep="")]]
		mc <- p[[paste("c", n, sep="")]]
		tot <- c(mr, mc)
		gr <- rep(c("r", "c"), c(length(mr), length(mc)))
		stot <- (tot-median(tot))/IQR(tot)
		m.z[[paste("c", n, sep="")]] <- stot[gr=="c"]
	}
	unlist(m.z)
}

compare.rnd2cart <- function(p, Qs=c(0,0.05, 0.25, 0.5, 0.75, 0.95, 1)) {
	d <- numeric(length(Qs))
	i <- 1
	for(Q in Qs) {
		d[i] <- cartel.diff(p, Q)
		i <- i+1
	}
	data.frame(Q=Qs, d=d)
}

opar <- par()

# MTMT

## read data

m.raw <- read.csv("../analyse-publications/MTMT/MTMTpubmat-productivity_raw.csv")
m.cart <- read.csv("../analyse-publications/MTMT/MTMTpubmat-cartel_productivity.csv")

## process and plot data


pdf(file="../paperfigs/group_productivity.pdf", width=7, height=10)
layout(matrix(1:8, ncol=2))
par(mar=c(4,4,0,1)+0.1)

prods <- list()
prods[["groupprod"]] <- plot.prod(m.raw, m.cart, "groupprod",
																	"group productivity")
prods[["npapers"]] <- plot.prod(m.raw, m.cart, "npapers",
																"mean number of papers")
prods[["wpapers"]] <- plot.prod(m.raw, m.cart, "wpapers",
																"mean weighted number of papers")

gp <- compare.rnd2cart(prods$groupprod, seq(0,1,0.1))
np <- compare.rnd2cart(prods$npapers, seq(0,1,0.1))
wp <- compare.rnd2cart(prods$wpapers, seq(0,1,0.1))

## plot data

plot(gp, xlim=c(0,1), ylim=c(0,1), type="o", pch=16,
		 xlab="quantiles of random samples",
		 ylab="proportion in subspicious groups")
grid()
abline(a=0, b=1, lty=2)
lines(np, type="o", pch=17, col=2)
lines(wp, type="o", pch=18, col=3)
legend("bottomright", legend=c("group productivity", "number of papers",
															 "weighted number of papers"), col=1:3,
			 pch=16:18, lty=1, bty="n", title="measures")

# dblp

## read data

m.raw <- read.csv("../analyse-publications/dblp/dblppubmat-productivity_raw.csv")
m.cart <- read.csv("../analyse-publications/dblp/dblppubmat-cartel_productivity.csv")

## process and plot data

prods <- list()
prods[["groupprod"]] <- plot.prod(m.raw, m.cart, "groupprod",
																	"group productivity")
prods[["npapers"]] <- plot.prod(m.raw, m.cart, "npapers",
																"mean number of papers")
prods[["wpapers"]] <- plot.prod(m.raw, m.cart, "wpapers",
																"mean weighted number of papers")

gp <- compare.rnd2cart(prods$groupprod, seq(0,1,0.1))
np <- compare.rnd2cart(prods$npapers, seq(0,1,0.1))
wp <- compare.rnd2cart(prods$wpapers, seq(0,1,0.1))

# plot data

plot(gp, xlim=c(0,1), ylim=c(0,1), type="o", pch=16,
		 xlab="quantiles of random samples",
		 ylab="proportion in subspicious groups")
grid()
abline(a=0, b=1, lty=2)
lines(np, type="o", pch=17, col=2)
lines(wp, type="o", pch=18, col=3)
legend("bottomright", legend=c("group productivity", "number of papers",
															 "weighted number of papers"), col=1:3,
			 pch=16:18, lty=1, bty="n", title="measures")

layout(1)
par(opar)
dev.off()

#z.gp = cartel.z.scores(prods$groupprod)
#z.np = cartel.z.scores(prods$npapers)
#z.wp = cartel.z.scores(prods$wpapers)
#
#boxplot(z.gp, z.np, z.wp, outline=FALSE)
#
#m.gp = cartel.m.scores(prods$groupprod)
#m.np = cartel.m.scores(prods$npapers)
#m.wp = cartel.m.scores(prods$wpapers)
#
#boxplot(m.gp, m.np, m.wp, outline=FALSE)
