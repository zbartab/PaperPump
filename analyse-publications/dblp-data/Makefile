
# prepare the publication database from dblp xml file

parsed-dblp-list.Rdata: parse-dblp_xml.R dblp-articles.xml
	#Rscript $< > /dev/null

dblp-articles.xml: dblp-nl.xml sort-articles-out.pl
	./sort-articles-out.pl $< > $@

dblp-nl.xml: dblp-2019-12-01.xml
	sed -e 's/></>\n</g' $< > $@
