
# figure to illustrate the weighted publication production

## functions
weighted.pubs.A <- function(c.co, a=3, G=4) {
	a/(1+c.co)
}
weighted.pubs.B <- function(c.co, a=3, G=4) {
	(G*a)/(G+c.co)
}
w.p.A <- function(co, a=3, G=4, b=0.2) {
	a*b + (1-b)*(a/(1+co))
}
w.p.B <- function(co, a=3, G=4, b=0.2) {
	a*b + (1-b)*((G*a)/(G+co))
}

## plots

#png(file="weighted_production.png", width=700, height=600)
pdf(file="paperfigs/weighted_production.pdf", width=10, height=8)
op <- par(mfcol=c(1, 2))
cc <- seq(0,10,0.01)
plot(cc, weighted.pubs.A(cc), type="l", bty="l",
		 xlab=expression(paste("number of collaborators, ", italic(c))),
		 ylab=expression(paste("weighted publication performance, ", italic(w)[i])))
lines(cc, weighted.pubs.B(cc), col=2)
lines(cc, w.p.A(cc), lty=2)
lines(cc, w.p.B(cc), lty=2, col=2)
legend("topright", legend=c("author A", "author A with bonus",
														"author B", "author B with bonus"), lty=c(1,2,1,2),
														col=c(1,1,2,2), bty="n")
mtext("a", side=3, adj=0.1)

ww0 <- weighted.pubs.B(cc)/weighted.pubs.A(cc)
ww02 <- w.p.B(cc)/w.p.A(cc)
ww05 <- w.p.B(cc, b=0.5)/w.p.A(cc, b=0.5)
plot(cc, ww0, type="l", bty="l",
		 xlab=expression(paste("number of collaborators, ", italic(c))),
		 ylab=expression(paste("relative weighted performance, ",
													 italic(w)[B]/italic(w)[A])))
lines(cc, ww02, lty=2)
lines(cc, ww05, lty=2)
legend("topleft", legend=c("no bonus", "first authorship bonus"),
			 lty=c(1,2), bty="n")
mtext(expression(italic(b) == "0.0"), side=4, at=ww0[length(ww0)])
mtext(expression(italic(b) == "0.2"), side=4, at=ww02[length(ww02)])
mtext(expression(italic(b) == "0.5"), side=4, at=ww05[length(ww05)])
mtext("b", side=3, adj=0.1)
par(op)
dev.off()

