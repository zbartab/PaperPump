


parsed.list <- "parsed-dblp-list.Rdata"

dblp.file <- scan(file="dblp-articles.txt", what=character(), sep="\n")
n.records <- length(dblp.file)/3
#n.records <- 10
dblp.ls <- vector(mode = "list", length = n.records)
i <- 1:n.records
i <- (i*3) - 2
n <- sub("^doi\\.org/", "", dblp.file[i])
names(dblp.ls) <- n
for(i in 1:n.records) {
	j <- i*3
	dblp.ls[[i]] <- list(authors=unlist(strsplit(dblp.file[j-1], "; ")),
											 year=as.numeric(dblp.file[j]))
	#if(i %% 10000 == 0) {print(i)}
}

save(dblp.ls, file=parsed.list)
