
# contains recipies to assemble PaperPump Appendix

SRC = PaperPump_Appendix.md

FPP = filepp -kc \& -mp !

DFIG = paperfigs
DSRC = papersrc
DEPSFIRST = $(wildcard $(DFIG)/real_nets_description.*)
DEPSSECOND = $(wildcard $(DFIG)/*-descr-degrees.txt)

PDFS=$(SRC:.md=.pdf)
HTML=$(SRC:.md=.html)
TEX=$(SRC:.md=.tex)
DOCX=$(SRC:.md=.docx)

PDFFIGS=$(wildcard $(DFIG)/*.pdf)
PNGFIGS=$(PDFFIGS:.pdf=.png)

pdf: $(PDFS)
html: $(HTML)
docx: $(DOCX)

$(DOCX): $(SRC) $(DEPSFIRST) $(DEPSSECOND)
	$(FPP) $< | pandoc -c ~/lib/markdown/pandoc.css --mathml -N --standalone \
		--self-contained --filter pandoc-citeproc -o $@

$(PDFS): $(SRC) $(DEPSFIRST) $(DEPSSECOND)
	$(FPP) -D"FIGURE(a)"="" -DEXT=pdf -Dpagebreak="\\newpage" $< | pandoc -N --standalone \
		--pdf-engine=xelatex --self-contained --filter pandoc-citeproc -o $@ 

$(HTML): $(SRC) $(PNGFIGS) $(DEPSFIRST) $(DEPSSECOND)
	$(FPP) -D"FIGURE(a)"="Figure a." -DEXT=png $< | pandoc \
		-c ~/lib/markdown/pandoc.css --mathml -N --standalone --toc \
		--self-contained --filter pandoc-citeproc -o $@

$(DFIG)/%.pdf: $(DSRC)/%.R
	Rscript $<

.PRECIOUS: $(DFIG)/%.pdf

$(DFIG)/%.png: $(DFIG)/%.pdf
	convert -density 100 $< -density 100 $@

$(DFIG)/real_nets_description.pdf: $(DSRC)/describe_real_nets.jl
	cd $(DSRC) && julia $(<F)

$(DFIG)/real_nets_description.txt: $(DSRC)/describe_real_nets.jl
	cd $(DSRC) && julia $(<F)

$(DFIG)/MTMT-descr-degrees.txt: $(DSRC)/describe_distributions.jl
	cd $(DSRC) && julia $(<F)

$(DFIG)/dblp-descr-degrees.txt: $(DSRC)/describe_distributions.jl
	cd $(DSRC) && julia $(<F)

