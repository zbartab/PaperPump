
# this file documents the production of publication and collaboration
# matrices from the dbpl dataset

PUBMAT = dblppubmat-full.csv
COLMAT = dblpcolmat-full.csv
DBLPDATA = dbpl/dbpl-articles.txt

.PHONY: all
all: $(PUBMAT) $(COLMAT) rewired

$(PUBMAT): make_pubmat.jl $(DBLPDATA)
	./make_pubmat.jl 0 1000000 $@ 

$(COLMAT): $(PUBMAT) make_colmat.jl 
	./make_colmat.jl $< 3

rewired: $(PUBMAT) rewire_pubmat.jl
	./rewire_pubmat.jl $< 20 1000000
	touch $@

$(DBLPDATA):
	$(MAKE) --directory=dblp-data
