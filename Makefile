
# contains recipies to assemble to PaperPump paper

SRC = PaperPump.md

FPP = filepp -kc \& -mp !

DFIG = paperfigs
DSCR = papersrc

PDFS=$(SRC:.md=.pdf)
HTML=$(SRC:.md=.html)
TEX=$(SRC:.md=.tex)
DOCX=$(SRC:.md=.docx)

PNGFIGS=$(DFIG)/*.png

html: $(HTML)
pdf: $(PDFS)
docx: $(DOCX)

$(DOCX): $(SRC) groups.png weighted_production.png
	$(FPP) $< | pandoc -c ~/lib/markdown/pandoc.css --mathml -N --standalone \
		--self-contained --filter pandoc-citeproc -o $@

$(PDFS): $(SRC)
	pandoc -o $@ $<

#$(HTML): $(SRC) $(DFIG)/groups.png $(DFIG)/weighted_production.png 
$(HTML): $(SRC) $(PNGFIGS)
	t_sample_graphs_tex
	$(FPP) -DEXT=png $< | pandoc -c ~/lib/markdown/pandoc.css --mathml \
		-N --standalone --toc --self-contained --filter pandoc-citeproc -o $@

$(DFIG)/%.pdf: $(DSCR)/%.R
	Rscript $<

.PRECIOUS: $(DFIG)/%.pdf

$(DFIG)/%.png: $(DFIG)/%.pdf
	convert -density 100 $< -density 100 $@

# produce sample_graph figures

t_sample_graphs_julia: $(DSCR)/sample_publication_network.jl
	cd $(DSCR) && julia $(<F)
	touch $@

t_sample_graphs_tex: $(DSCR)/sample_publication_network.tex \
	t_sample_graphs_julia
	$(FPP) -I paperfigs/ $< > work/$(<F)
	cd work && xetex $(<F) && mv $(<F:.tex=.pdf) ../$(DFIG)
	touch $@
