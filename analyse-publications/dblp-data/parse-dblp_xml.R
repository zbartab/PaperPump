
library(XML)

parsed.object <- "parsed-dblp.Rdata"
d <- xmlParse("dblp-2019-12-01.xml")
dblp.ls <- xmlToList(d)
save(dblp.ls, file=parsed.object)
