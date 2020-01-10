
# this file documents the production of publication and collaboration
# matrices from the dbpl dataset

PUBMAT = dblppubmat-full.csv
COLMAT = dblpcolmat-full.csv
DBLPDATA = dbpl/dbpl-articles.txt
MK_PUBMAT = make_pubmat.jl
MK_COLMAT = make_colmat.jl
REWIRE = rewire_pubmat.jl

.PHONY: all
all: $(PUBMAT) $(COLMAT) rewired

$(PUBMAT): $(MK_PUBMAT) $(DBLPDATA)
	./$(MK_PUBMAT) 0 1000000 $@ 

$(COLMAT): $(PUBMAT) $(MK_COLMAT) 
	./$(MK_COLMAT) $< 3

rewired: $(PUBMAT) $(REWIRE)
	./$(REWIRE) $< 2 10
	touch $@

$(DBLPDATA):
	$(MAKE) --directory=dblp-data
