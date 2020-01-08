
# this file documents the production of publication and collaboration
# matrices from the dbpl dataset

.PHONY: all
all: dblppubmat.csv dblpcolmat-full.csv rewired

dblppubmat.csv: make_pubmat.jl dblp-data/dblp-articles.txt
	./make_pubmat.jl 2001 2010

dblpcolmat-full.csv: dblppubmat-full.csv make_colmat.jl 
	./make_colmat.jl $<

rewired: dblppubmat.csv rewire_pubmat.jl
	./rewire_pubmat.jl $< 20 1000000
	touch $@

dblp-data/dblp-articles.txt:
	$(MAKE) --directory=dblp-data
