
# contains recipies to assemble to PaperPump paper

SRC = PaperPump.md

FPP = filepp -kc \& -mp !

DFIG = paperfigs
DSCR = papersrc

PDFS=$(SRC:.md=.pdf)
HTML=$(SRC:.md=.html)
TEX=$(SRC:.md=.tex)
DOCX=$(SRC:.md=.docx)

PDFFIGS=$(wildcard $(DFIG)/*.pdf)
PNGFIGS=$(PDFFIGS:.pdf=.png)

pdf: $(PDFS)
html: $(HTML)
docx: $(DOCX)

$(DOCX): $(SRC) groups.png weighted_production.png
	$(FPP) $< | pandoc -c ~/lib/markdown/pandoc.css --mathml -N --standalone \
		--self-contained --filter pandoc-citeproc -o $@

$(PDFS): $(SRC) t_sample_graphs_tex t_simulation_analyses_tex
	$(FPP) -D"FIGURE(a)"="" -DEXT=pdf $< | pandoc -N --standalone \
		--self-contained --filter pandoc-citeproc -o $@ 

$(HTML): $(SRC) $(PNGFIGS) t_sample_graphs_tex t_simulation_analyses_tex
	$(FPP) -D"FIGURE(a)"="Figure a." -DEXT=png $< | pandoc \
		-c ~/lib/markdown/pandoc.css --mathml -N --standalone --toc \
		--self-contained --filter pandoc-citeproc -o $@

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
	cd work && xetex $(<F)
	pdftk work/$(<F:.tex=.pdf) burst output $(DFIG)/$(<F:.tex=-%02d.pdf)
	touch $@

t_simulation_analyses: $(DSCR)/simulation_analyses.jl
	cd $(DSCR) && julia $(<F)
	touch $@

t_simulation_analyses_tex: $(DSCR)/simulation_analyses.tex \
	t_simulation_analyses
	$(FPP) -I paperfigs/ $< > work/$(<F)
	cd work && xetex $(<F)
	pdftk work/$(<F:.tex=.pdf) burst output $(DFIG)/$(<F:.tex=-%02d.pdf)
	touch $@

