
# this file documents the production of publication and collaboration
# matrices from the dbpl dataset

PUBMAT = dblppubmat.txt
COLMAT = dblpcolmat.txt
DBLPDATA = dbpl/dbpl-articles.txt
MK_PUBMAT = make_pubmat.jl
MK_COLMAT = make_colmat.jl
REWIRE = parallel_rewire_pubmat.jl
NPAPERS = 3

.PHONY: all
all: $(PUBMAT) $(COLMAT) rewired

$(PUBMAT): $(MK_PUBMAT) $(DBLPDATA)
	./$(MK_PUBMAT) 0 1000000 $@ 

$(COLMAT): $(PUBMAT) $(MK_COLMAT) 
	./$(MK_COLMAT) $< $(NPAPERS)

rewired: $(PUBMAT) $(REWIRE)
	./$(REWIRE) $< 20 1000000 $(NPAPERS)
	touch $@

$(DBLPDATA):
	$(MAKE) --directory=dblp-data
