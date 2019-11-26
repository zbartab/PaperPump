
# functions to process MTMT publication records json

library("rjson")

read.MTMT <- function(file) {
	mtmt.json <- fromJSON(file=file)
	mtmt.json$content
}

p.rank <- function(cikk) {
	r <- cikk$ratings[[1]]$ranking
	if(is.null(r)) r <- NA
	r
}

p.n.authors <- function(cikk) {
	length(cikk$authorships)

}

p.authors <- function(cikk) {
	a <- cikk$authorships
	givenName <- sapply(a, function(aa) aa$givenName)
	familyName <- sapply(a, function(aa) aa$familyName)
	paste(givenName, familyName, sep="-")
}

p.author.rank <- function(cikk, id) {
	a <- cikk$authorships
	r <- sapply(a, function(l) {
							if(length(l$author) > 1) {l$author$mtid} else {NA}})
	o <- which(r == id)
	o[1]
}

p.author.label <- function(cikk, id) {
	a <- cikk$authorships
	r <- sapply(a, function(l) {
							if(length(l$author) > 1) {l$author$mtid} else {NA}})
	o <- which(r == id)
	i <- o[1]
	a[[o[1]]]$label
}

p.n.citations <- function(cikk) {
	cikk$independentCitationCount
}

pub.measures <- function(cikkek, au.id) {
	#cikkek <- read.MTMT(file)
	id <- sapply(cikkek, function(p) p$mtid)
	types <- sapply(cikkek, function(p) p$type$label)
	year <- sapply(cikkek, function(p) p$publishedYear)
	ranks <- sapply(cikkek, p.rank)
	citations <- sapply(cikkek, p.n.citations)
	n.authors <- sapply(cikkek, p.n.authors)
	r.authors <- sapply(cikkek, p.author.rank, au.id)
	d <- data.frame(id=id, types=types, year=year, ranks=ranks,
									citations=citations, n.authors=n.authors,
									r.author=r.authors, first.au=r.authors == 1,
									last.au=r.authors == n.authors)
	d[!is.na(d$ranks),]
}

create.fileID <- function(staff.record) {
	# create a hopefully unique file name to store a staff member's
	# publication records
	n <- staff.record[1]
	n <- gsub("\\<(.)[^ ]+", "\\1", n)
	n <- gsub(" ", "", n)
	fileID <- paste(n, staff.record[3], staff.record[2], sep="-")
	fileID
}

create.nameID <- function(staff.record) {
	# create an ID from the name and the department of a staff member
	n <- staff.record[1]
	n <- gsub("\\<(..)[^ ]+", "\\1", n)
	n <- gsub(" ", "", n)
	d <- sub("^(..).*", "\\1", staff.record[3])
	ID <- paste(n, d, sep="-")
	ID
}

download.MTMTrecords <- function(mtmt.id, outfile) {
	mtmt.url <- "https://m2.mtmt.hu/api/publication?cond=published%3Beq%3Btrue&cond=core%3Beq%3Btrue&cond=authors.mtid%3Beq%3B<-mtmtid->&ty_on=1&ty_on_check=1&st_on=1&st_on_check=1&url_on=1&url_on_check=1&cite_type=2&sort=publishedYear%2Cdesc&sort=firstAuthor%2Casc&size=5000&page=1&format=json"
	mtmt.url <- sub("<-mtmtid->", mtmt.id, mtmt.url)
	#outfile <- paste(outfileID, "json", sep=".")
	download.file(mtmt.url, outfile, method="libcurl")
}

calc.bib.measures <- function(cikkek) {
	# total number of independent citations:
	cit <- sum(cikkek$citations)
	# total number of coauthors:
	coauth <- sum(cikkek$n.authors)
	# number of citations corrected for the number of coauthors:
	cit.coauth <- sum(cikkek$citations/cikkek$n.authors)
	# number of papers:
	n.papers <- nrow(cikkek)
	# number of papers corrected for the number of coauthors:
	n.papers.corr <- sum(1/cikkek$n.authors)
	# number of ranked papers:
	n.WoS <- sum(!is.na(cikkek$ranks))
	# number of ranked papers corrected for the number of coauthors:
	n.WoS.corr <- sum(1/cikkek$n.authors[!is.na(cikkek$ranks)])
	# number of D1 papers
	n.D1 <- sum(cikkek$ranks == "D1", na.rm=TRUE)
	n.D1.corr <- sum(1/cikkek$n.authors[cikkek$ranks == "D1"], na.rm=TRUE)
	# proportion of D1 papers to the number of ranked papers:
	p.D1 <- n.D1/n.WoS
	c(n.citations=cit,
		n.coauthors=coauth,
		avg.citations=cit/n.papers, # average number of citations per paper
		n.citations.corr=cit.coauth,
		n.papers=n.papers,
		n.papers.corr=n.papers.corr,
		n.WoS=n.WoS,
		n.WoS.corr=n.WoS.corr,
		n.D1=n.D1,
		n.D1.corr=n.D1.corr,
		p.D1=p.D1)
}

create.directed.assoc.mat <- function() {
	# create a directed association matrix based on common papers
	staff <- names(staff.measures)
	n.staff <- length(staff)
	assoc.mat <- matrix(0, ncol=n.staff, nrow=n.staff)
	staffID <- apply(BOI.staff, 1, create.nameID)
	rownames(assoc.mat) <- staff
	colnames(assoc.mat) <- staff
	for(i in 1:(n.staff-1)) {
		for(j in (i+1):n.staff) {
			m1.ps <- staff.measures[[staff[i]]]$id
			m2.ps <- staff.measures[[staff[j]]]$id
			assoc.mat[staff[i], staff[j]] <- length(intersect(m1.ps,
																												m2.ps))/length(m2.ps)
			assoc.mat[staff[j], staff[i]] <- length(intersect(m2.ps,
																												m1.ps))/length(m1.ps)
		}
	}
	rownames(assoc.mat) <- staffID
	colnames(assoc.mat) <- staffID
	assoc.mat
}

create.assoc.mat <- function(staff.measures) {
	# create an undirected association matrix based on common papers
	staff <- names(staff.measures)
	n.staff <- length(staff)
	assoc.mat <- matrix(0, ncol=n.staff, nrow=n.staff)
	#i <- BOI.staff$mtmt.id %in% staff
	#staffID <- apply(BOI.staff[i,], 1, create.nameID)
	rownames(assoc.mat) <- staff
	colnames(assoc.mat) <- staff
	for(i in 1:(n.staff-1)) {
		for(j in (i+1):n.staff) {
			m1.ps <- staff.measures[[staff[i]]]$id
			m2.ps <- staff.measures[[staff[j]]]$id
			assoc.mat[staff[i], staff[j]] <- 
				length(intersect(m1.ps, m2.ps))/length(union(m1.ps, m2.ps))
		}
	}
	#rownames(assoc.mat) <- staffID
	#colnames(assoc.mat) <- staffID
	assoc.mat
}

comp.authors <- function(a1, a2) {
	# compare two author label and try to decide if they belong to the
	# same author
	if(a1 == a2) {
		return(TRUE)
	} else {
		# the same names but in different order
		aa <- sort(unlist(strsplit(a1, " ")))
		l.aa <- length(aa)
		ab <- sort(unlist(strsplit(a2, " ")))
		l.ab <- length(ab)
		if(l.aa == l.ab) {
			if(sum(aa == ab) == l.aa) {
				return(TRUE)
			}
		}
	}
	return(FALSE)
}

hash.id <- function(x, base=36, n.digits=5) {
	# convert `x` from base ten to base `base`
	digits <- c(0:9, letters[1:(base-10)])
	n <- floor(log(x, base))
	res <- ""
	for(i in n:0) {
		j <- x %/% base^i
		res <- paste(res, digits[j+1], sep="")
		x <- x %% base^i
	}
	paste(paste(rep(0, n.digits-nchar(res)), collapse=""), res, sep="")
}
