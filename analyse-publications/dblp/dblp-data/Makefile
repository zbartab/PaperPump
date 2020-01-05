
# prepare the publication database from dblp xml file

dblp-articles.txt: dblp-nl.xml sort-articles-out.pl
	./sort-articles-out.pl $< > $@

dblp-nl.xml: dblp-2019-12-01.xml
	sed -e 's/></>\n</g' $< > $@

dblp-2019-12-01.xml:
	wget https://dblp.org/xml/release/dblp-2019-12-01.xml.gz .
	wget https://dblp.org/xml/release/dblp-2019-11-22.dtd .
	gunzip	dblp-2019-12-01.xml.gz
