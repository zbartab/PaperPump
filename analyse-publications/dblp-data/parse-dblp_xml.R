
library(XML)


parsed.object <- "parsed-dblp.Rdata"
parsed.list <- "parsed-dblp-list.Rdata"
d.parsed <- xmlParse("dblp-articles.xml")
save(d.parsed, file=parsed.object)
dblp.ls <- xmlToList(d.parsed)
save(dblp.ls, file=parsed.list)
