
# this file documents the production of publication and collaboration
# matrices from the dbpl dataset

PUBMAT = MTMTpubmat.txt
COLMAT = MTMTcolmat.txt
MK_PUBMAT = make_pubmat.jl
MK_COLMAT = make_colmat.jl
REWIRE = rewire_pubmat.jl

.PHONY: all
all: $(PUBMAT) $(COLMAT) rewired

$(PUBMAT): $(MK_PUBMAT)
	./$(MK_PUBMAT) $@ 

$(COLMAT): $(PUBMAT) $(MK_COLMAT)
	./$(MK_COLMAT) $< 3

rewired: $(PUBMAT) $(REWIRE)
	./$(REWIRE) $< 20 1000000
	touch $@

