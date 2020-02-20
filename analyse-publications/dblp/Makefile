
# this file documents the production of publication and collaboration
# matrices from the dbpl dataset

PUBMAT = dblppubmat.mat
COLMAT = dblpcolmat.mat
DBLPDATA = dbpl/dbpl-articles.txt
MK_PUBMAT = make_pubmat.jl
MK_COLMAT = make_colmat.jl
REWIRE = parallel_rewire_pubmat.jl
PARTREWIRE = ../pl_partial_rewire.jl
NPAPERS = 3
TREES = $(patsubst %.mat,%.tree,$(wildcard *colmat*.mat))

.PHONY: all
all: $(PUBMAT) $(COLMAT) $(TREES) rewired partrewired

#$(PUBMAT): $(MK_PUBMAT) $(DBLPDATA)
$(PUBMAT): $(MK_PUBMAT) dblpdatadone
	./$(MK_PUBMAT) 0 1000000 $@ 

$(COLMAT): $(PUBMAT) $(MK_COLMAT) 
	./$(MK_COLMAT) $< $(NPAPERS)

rewired: $(PUBMAT) $(REWIRE)
	./$(REWIRE) $< 20 1000000 $(NPAPERS)
	touch $@

partrewired: $(PUBMAT) $(PARTREWIRE)
	$(PARTREWIRE) $< 20 $(NPAPERS)
	touch $@

#$(DBLPDATA): dblpdatadone
dblpdatadone:
	$(MAKE) --directory=dblp-data
	touch dblpdatadone



%.tree: %.mat
		../comm_detection.sh $<

