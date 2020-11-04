
# fit power-law and similar distributions


# load MTMT degree data, saved from julia {{{1

#d = read.csv("MTMT-degree.csv")
#d = d$degrees
#d = d[d>0]

# load code {{{1

source("../analyse-publications/pli-R-v0.0.3-2007-07-25/powerexp.R")
source("../analyse-publications/pli-R-v0.0.3-2007-07-25/discexp.R")
source("../analyse-publications/pli-R-v0.0.3-2007-07-25/discpowerexp.R")
source("../analyse-publications/pli-R-v0.0.3-2007-07-25/zeta.R")
source("../analyse-publications/pli-R-v0.0.3-2007-07-25/disclnorm.R")
source("../analyse-publications/pli-R-v0.0.3-2007-07-25/weibull.R")
source("../analyse-publications/pli-R-v0.0.3-2007-07-25/discweib.R")
source("../analyse-publications/pli-R-v0.0.3-2007-07-25/pareto.R")

# functions {{{1

## functions to fit heavy tailed distributions {{{2

fit.weibull <- function(d, x.min) {
	zd <- try(discweib.fit(d, x.min), silent=TRUE)
	if("try-error" %in% class(zd)) {
		return(NA)
	}
	x <- 0:(100*max(d))
	fd <- ddiscweib(x, zd$shape, zd$scale, zd$threshold)
	i <- !is.na(fd)
	if(sum(i) == 0) {
		return(NA)
	}
	fd <- cumsum(fd[i])
	#fd <- c(0, fd[-length(fd)])
	list(zd=zd, fd=fd, x=x[i])
}

fit.exponential <- function(d, x.min) {
	zd <- discexp.fit(d, x.min)
	x <- 0:(100*max(d))
	fd <- ddiscexp(x, zd$lambda, zd$threshold)
	i <- !is.na(fd)
	fd <- cumsum(fd[i])
	#fd <- c(0, fd[-length(fd)])
	list(zd=zd, fd=fd, x=x[i])
}

fit.powerlaw <- function(d, x.min) {
	zd <- try(zeta.fit(d, x.min), silent=TRUE)
	if("try-error" %in% class(zd)) {
		return(NA)
	}
	x <- 0:(100*max(d))
	fd <- dzeta(x, zd$threshold, zd$exponent)
	i <- !is.na(fd)
	fd <- cumsum(fd[i])
	#fd <- c(0, fd[-length(fd)])
	list(zd=zd, fd=fd, x=x[i])
}

fit.lognormal <- function(d, x.min) {
	zd <- try(fit.lnorm.disc(d, x.min), silent=TRUE)
	if("try-error" %in% class(zd)) {
		return(NA)
	}
	x <- 0:(100*max(d))
	fd <- dlnorm.disc(x, zd$meanlog, zd$sdlog)
	i <- !is.na(fd) & x >= x.min
	fd <- fd[i]/sum(fd[i])
	fd <- cumsum(fd)
	#fd <- c(0, fd[-length(fd)])
	list(zd=zd, fd=fd, x=x[i])
}

fit.powerlawexp <- function(d, x.min) {
	zd <- try(discpowerexp.fit(d, x.min), silent=TRUE)
	if("try-error" %in% class(zd)) {
		return(NA)
	}
	x <- 0:(100*max(d))
	fd <- ddiscpowerexp(x, zd$exponent, zd$rate, zd$threshold)
	i <- !is.na(fd)
	fd <- cumsum(fd[i])
	#fd <- c(0, fd[-length(fd)])
	list(zd=zd, fd=fd, x=x[i])
}

fitter.funs <- list(lognormal=fit.lognormal, weibull=fit.weibull, exponential=fit.exponential,
										powerlaw=fit.powerlaw, powerlawexp=fit.powerlawexp)

## functions to find x.min {{{2

find.x.max <- function(d, n.x.max=100) {
	# find the maximum of `d` for which `length(d[d >= x.max]) > n.x.max`
	t.d <- table(d)
	cs.d <- cumsum(rev(t.d))
	i <- which(cs.d >= n.x.max)
	return(as.numeric(names(cs.d)[i])[1])
}

find.x.min <- function(d, x.min, fitter=fit.powerlaw, x.max=NULL) {
	if(length(x.min) == 1) {x.min = c(0, x.min)}
	if(is.null(x.max)) {x.max = find.x.max(d)}
	if(x.min[2] > x.max) {x.min[2] = x.max}
	u.d <- sort(unique(d))
	xs <- u.d[x.min[1] < u.d & u.d <= x.min[2]]
	Ds <- rep(NA, length(xs))
	mD <- 1.0e12
	mx <- 0
	fitted <- NA
	for(i in 1:length(xs)) {
		print(i)
		r <- fitter(d, xs[i])
		if(is.na(r)) {next}
		rD <- KS.D(d, r)
		Ds[i] <- rD
		if(rD < mD) {
			mD = rD
			mx = xs[i]
			fitted = r
		}
	}
	#list(D=mD, x=mx, fitted=fitted, Ds=Ds, xs=xs)
	list(D=mD, x=mx, fitted=fitted, x.max=x.max)
}

KS.D <- function(d, fit) {
	# `fit` is returned by `fit.` functions
	md <- d[d >= fit$zd$threshold]
	Fd <- ecdf(md)
	ud <- sort(unique(md))
	i <- fit$x %in% ud
	D = abs(Fd(ud) - fit$fd[i])
	return(max(D))
}

## functions to plot the fits and the distributions {{{2

plot.powerlaw <- function(d, zd) {
	dm <- d[d >= zd$threshold]
	ud <- sort(unique(dm))
	Fd <- ecdf(dm)
	plot(ud, 1-Fd(ud), type="s", log="xy")
	lines(ud, 1-pzeta(ud, zd$threshold, zd$exponent), col=2, lty=2)
}

coord.obs <- function(d, x.min) {
	md <- d[d >= x.min]
	Fd <- ecdf(md)
	ud <- sort(unique(md))
	cbind(ud, 1-Fd(ud))
}

coord.distrib <- function(d, fit) {
	md <- d[d >= fit$zd$threshold]
	ud <- sort(unique(md))
	i <- fit$x %in% ud
	cbind(fit$x[i], 1-fit$fd[i])
}

## functions to generate random numbers from heavy tailed distributions

rand.weibull <- function(d, fit) {
	n <- length(d)
	i <- d >= fit$zd$threshold
	rthreshold <- fit$zd$threshold-0.5
	ntail <- sum(i)
	ptail <- ntail/n
	ne <- n * (1/pweibull(rthreshold, fit$zd$shape, fit$zd$scale))
	rtail <- NULL
	while(length(rtail) < n) {
		rt <- rweibull(ne, fit$zd$shape, fit$zd$scale)
		rt <- rt[rt >= rthreshold]
		rtail <- c(rtail, rt)
	}
	rtail <- round(rtail[1:n])
	if(sum(!i) > 0) {
		r <- sample(d[!i], n, replace=TRUE)
		ii <- runif(n) < ptail
		return(c(r[!ii], rtail[ii]))
	} else {
		return(rtail)
	}
}

rand.exponential <- function(d, fit) {
	n <- length(d)
	i <- d >= fit$zd$threshold
	rthreshold <- fit$zd$threshold-0.5
	ntail <- sum(i)
	ptail <- ntail/n
	ne <- n * (1/pexp(rthreshold, fit$zd$lambda))
	rtail <- NULL
	while(length(rtail) < n) {
		rt <- rexp(ne, fit$zd$lambda)
		rt <- rt[rt >= rthreshold]
		rtail <- c(rtail, rt)
	}
	rtail <- round(rtail[1:n])
	if(sum(!i) > 0) {
		r <- sample(d[!i], n, replace=TRUE)
		ii <- runif(n) < ptail
		return(c(r[!ii], rtail[ii]))
	} else {
		return(rtail)
	}
}

rand.lognormal <- function(d, fit) {
	n <- length(d)
	i <- d >= fit$zd$threshold
	rthreshold <- fit$zd$threshold-0.5
	ntail <- sum(i)
	ptail <- ntail/n
	ne <- n * (1/plnorm(rthreshold, fit$zd$meanlog, fit$zd$sdlog))
	rtail <- NULL
	while(length(rtail) < n) {
		rt <- rlnorm(ne, fit$zd$meanlog, fit$zd$sdlog)
		#print(length(rt))
		rt <- rt[rt >= rthreshold]
		rtail <- c(rtail, rt)
		#print(length(rtail))
	}
	rtail <- round(rtail[1:n])
	if(sum(!i) > 0) {
		r <- sample(d[!i], n, replace=TRUE)
		ii <- runif(n) < ptail
		return(c(r[!ii], rtail[ii]))
	} else {
		return(rtail)
	}
}

rand.powerlaw <- function(d, fit) {
	n <- length(d)
	i <- d >= fit$zd$threshold
	ntail <- sum(i)
	ptail <- ntail/n
	rtail <- round(rpareto(n, fit$zd$threshold-0.5, fit$zd$exponent))
	if(sum(!i) > 0) {
		r <- sample(d[!i], n, replace=TRUE)
		ii <- runif(n) < ptail
		return(c(r[!ii], rtail[ii]))
	} else {
		return(rtail)
	}
}

rand.powerlawexp <- function(d, fit) {
	n <- length(d)
	i <- d >= fit$zd$threshold
	ntail <- sum(i)
	ptail <- ntail/n
	rtail <- round(rpowerexp(n, fit$zd$threshold-0.5, fit$zd$exponent,
													 fit$zd$rate))
	if(sum(!i) > 0) {
		r <- sample(d[!i], n, replace=TRUE)
		ii <- runif(n) < ptail
		return(c(r[!ii], rtail[ii]))
	} else {
		return(rtail)
	}
}

rgen.funs <- list(powerlaw=rand.powerlaw, lognormal=rand.lognormal,
									exponential=rand.exponential, weibull=rand.weibull,
									powerlawexp=rand.powerlawexp)

## functions to calculate goodness-of-fit {{{2

goodness.of.fit <- function(d, fit, x.max, nrep=10, distribution="powerlaw") {
	x.mins <- round(fit$zd$threshold * c(0.5, 1.5))
	Dobs <- KS.D(d, fit)
	Ds <- rep(NA, nrep)
	for(i in 1:nrep) {
		cat("GOF (", distribution, "): ", i, "\n", sep="")
		r <- rgen.funs[[distribution]](d, fit)
		if(sum(is.na(r)) > 0) {
			Ds[i] <- NA
		} else {
			#cat("GOF generating random deviates finished.\n")
			f <- find.x.min(r, x.mins, fitter.funs[[distribution]], x.max)
			Ds[i] <- f$D
		}
	}
	Ds <- Ds[!is.na(Ds)]
	if(length(Ds) > 0) {
		return(list(p=sum(Ds > Dobs)/length(Ds), Dobs=Dobs, Drand=Ds))
	} else {
		return(list(p=NA, Dobs=Dobs, Drand=rep(NA, length(d))))
	}
}

## the main work horse function

my.format <- function(n, digits) {
	if(is.na(n)) {
		return(NA)
	} else {
		formatC(n, digits=digits, format="f")
	}
}

param.str <- function(distribution, fitted) {
	return(switch(distribution,
				 powerlaw = paste("$\\\\gamma = ", my.format(fitted$exponent,3),"$",
													sep=""),
				 powerlawexp = paste("$\\\\gamma = ", my.format(fitted$exponent,3),
														 "$, $\\\\lambda = ",
														 my.format(fitted$rate,3),"$", sep=""),
				 lognormal = paste("log(mean) = ", my.format(fitted$meanlog,3),
													 ", log(SD) = ", my.format(fitted$sdlog,3), sep=""),
				 weibull = paste("$k = ", my.format(fitted$shape,3), 
												 "$, $\\\\lambda = ",
												 my.format(fitted$scale, 3),"$", sep=""),
				 exponential = paste("$\\\\lambda = ", my.format(fitted$lambda, 3),
														 "$", sep="")
				 ))
}

heavy.tailed <- function(d, x.min, nrep=10, filename=NULL, variable="") {
	if(is.null(filename)) stop("ERROR: filename must be given!")
	#sink(filename)
	dist.names <- names(fitter.funs)
	cat(" **Ditribution** | _p_ | _D_~obs~ | _x_~min~ | _n_ | **Parameters**\n", file=filename)
	cat(":------------|:----:|:----:|:-----:|:----:|:-------------\n",
			file=filename, append=TRUE)
	for(dn in dist.names) {
		cat("Distribution:", dn, '\n')
		r1 <- find.x.min(d, x.min, fitter.funs[[dn]])
		if(length(r1$fitted)==1) {
			cat("x.min was not found\n")
			cat(dn, "|", my.format(NA, 3), "|", my.format(NA,3), "|",
					NA, "|", NA, "|", NA,
					"\n", sep=" ", file=filename, append=TRUE)
			cat("Results written\n")
		} else {
			cat("x.min found\n")
			print(r1$fitted$zd)
			r2 <- goodness.of.fit(d, r1$fitted, r1$x.max, nrep, dn)
			cat("goodness of fit calculated\n")
			cat(dn, "|", my.format(r2$p, 3), "|", my.format(r2$Dobs,3), "|",
					r1$x, "|", sum(d >= r1$x), "|", param.str(dn, r1$fitted$zd),
					"\n", sep=" ", file=filename, append=TRUE)
			cat("Results written\n")
		}
	}
	#sink()
}

# fit power law with exponential cut-off {{{1

sandbox <- function() {

d.sat <- find.x.min(d, 100)

ypl <- fit.powerlaw(d, d.sat$x)
yexp <- fit.powerlawexp(d, d.sat$x)
yln <- fit.lognormal(d, d.sat$x)
yexpo <- fit.exponential(d, d.sat$x)
yweib <- fit.weibull(d, d.sat$x)

c.obs <- coord.obs(d, 1)
c.pl <- coord.distrib(d, ypl)
c.exp <- coord.distrib(d, yexp)
c.ln <- coord.distrib(d, yln)
c.expo <- coord.distrib(d, yexpo)
c.weib <- coord.distrib(d, yweib)

plot(c.obs, log="xy")
lines(c.pl, col=2)
lines(c.exp, col=3)
lines(c.ln, col=4)
lines(c.expo, col=5)
lines(c.weib, col=6)
legend("bottomleft", legend=c("observations", "powerlaw",
															"powerlaw-exponential cutoff",
															"lognormal", "exponential", "weibull"),
			 col=1:6, lty=1, bty="n")

# fit heavy tailed distributions to data {{{1

yexp1 <- fit.powerlawexp(d, 1)

c.obs <- coord.obs(d, 1)
c.ln <- coord.distrib(d, yexp1)

plot(c.obs, log="xy")
lines(c.ln, col=4)

yln1 <- fit.lognormal(d, 1)

c.obs <- coord.obs(d, 1)
c.ln <- coord.distrib(d, yln1)

plot(c.obs, log="xy")
lines(c.ln, col=4)

yexpo1 <- fit.exponential(d, 1)

c.obs <- coord.obs(d, 1)
c.ln <- coord.distrib(d, yexpo1)

plot(c.obs, log="xy")
lines(c.ln, col=4)

yweib1 <- fit.weibull(d, 1)

c.obs <- coord.obs(d, 1)
c.ln <- coord.distrib(d, yweib1)

plot(c.obs, log="xy")
lines(c.ln, col=4)


ypl.1 <- fit.powerlaw(d, 1)
yexp.1 <- fit.powerlawexp(d, 1)
yln.1 <- fit.lognormal(d, 1)
yexpo.1 <- fit.exponential(d, 1)
yweib.1 <- fit.weibull(d, 1)

c.obs <- coord.obs(d, 1)
c.pl <- coord.distrib(d, ypl.1)
c.exp <- coord.distrib(d, yexp.1)
c.ln <- coord.distrib(d, yln.1)
c.expo <- coord.distrib(d, yexpo.1)
c.weib <- coord.distrib(d, yweib.1)

plot(c.obs, log="xy")
lines(c.pl, col=2)
lines(c.exp, col=3)
lines(c.ln, col=4)
lines(c.expo, col=5)
lines(c.weib, col=6)
legend("bottomleft", legend=c("observations", "powerlaw",
															"powerlaw-exponential cutoff",
															"lognormal", "exponential", "weibull"),
			 col=1:6, lty=1, bty="n")


mpl <- find.x.min(d, 100)

c.obs <- coord.obs(d, mpl$x)
c.pl <- coord.distrib(d, mpl$fitted)
plot(c.obs, log="xy")
lines(c.pl, col=2)

mplexp <- find.x.min(d, 100, fit.powerlawexp)

c.obs <- coord.obs(d, mplexp$x)
c.pl <- coord.distrib(d, mplexp$fitted)
plot(c.obs, log="xy")
lines(c.pl, col=2)

mln <- find.x.min(d, 100, fit.lognormal)

c.obs <- coord.obs(d, mln$x)
c.pl <- coord.distrib(d, mln$fitted)
plot(c.obs, log="xy")
lines(c.pl, col=2)

mweib <- find.x.min(d, 100, fit.weibull)

c.obs <- coord.obs(d, mweib$x)
c.pl <- coord.distrib(d, mweib$fitted)
plot(c.obs, log="xy")
lines(c.pl, col=2)

mexp <- find.x.min(d, 100, fit.exponential)

c.obs <- coord.obs(d, mexp$x)
c.exp <- coord.distrib(d, mexp$fitted)
plot(c.obs, log="xy")
lines(c.exp, col=2)



source("pli-R-v0.0.3-2007-07-25/power-law-test.R")

power.powerexp.lrt(mpl$fitted, mplexp$fitted)

d10 = d[d >= 10]
F10 = ecdf(d10)
zd10 = zeta.fit(d10, 10)
u10 = sort(unique(d10))

plot(u10, 1-F10(u10), type="s", log="xy")
lines(u10, 1-pzeta(u10, 10, zd10$exponent), type="s", col=2)

D = abs(F10(u10) - pzeta(u10, 10, zd10$exponent))
plot(D, type="l")


lnd = fit.lnorm.disc(d, 10)
u1 = sort(unique(d))
F1 = ecdf(d)
plot(u1, 1-F1(u1), type="s", log="xy")
lines(u1, 1-plnorm(u1, lnd$meanlog, lnd$sdlog), type="s", col=2)
}
