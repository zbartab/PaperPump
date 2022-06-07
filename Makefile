
# contains recipies to assemble to PaperPump paper

SRC = PP-MS_metrics.md

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

edit:
	gvim $(SRC)

view:
	zathura $(PDFS) &

#cp -a $< ~/Dropbox/Draft/work/

2cloud: $(SRC)
	if [ $< -nt ~/Dropbox/Draft/work/$(<F) ]; then \
		cp -a $< ~/Dropbox/Draft/work/ ; \
	else \
		echo Dropbox file is newer, do not copy!;\
	fi

patch:
	diff -u $(SRC) ~/Dropbox/Draft/work/$(SRC) > \
		work/$(SRC).patch || exit 0
	patch < work/$(SRC).patch
	touch $(SRC)

$(DOCX): $(SRC) $(PDFS)
	$(FPP) -DEXT=png $< | pandoc -c ~/lib/markdown/pandoc.css --mathml \
		-N --standalone --self-contained --filter pandoc-citeproc -o $@

$(PDFS): $(SRC) t_sample_graphs_tex t_simulation_analyses_tex \
	t_cartel_footprint_tex t_real_networks_tex t_group_productivity_fig
	$(FPP) -D"FIGURE(a)"="" -DEXT=pdf $< | pandoc -N --standalone \
		--pdf-engine=xelatex --self-contained --filter pandoc-citeproc -o $@ 

$(HTML): $(SRC) $(PNGFIGS) t_sample_graphs_tex t_simulation_analyses_tex \
	t_cartel_footprint_tex
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

t_cartel_footprint: $(DSCR)/cartel_footprint.jl
	cd $(DSCR) && julia $(<F)
	touch $@

t_cartel_footprint_tex: $(DSCR)/cartel_footprint.tex \
	t_cartel_footprint
	$(FPP) -I paperfigs/ $< > work/$(<F)
	cd work && xetex $(<F)
	pdftk work/$(<F:.tex=.pdf) burst output $(DFIG)/$(<F:.tex=-%02d.pdf)
	touch $@

t_real_networks: $(DSCR)/real_networks.jl
	cd $(DSCR) && julia $(<F)
	touch $@

t_real_networks_tex: $(DSCR)/real_networks.tex t_real_networks
	$(FPP) -I paperfigs/ $< > work/$(<F)
	cd work && xetex $(<F)
	mv work/$(<F:.tex=.pdf) $(DFIG)/
	touch $@

t_group_productivity: $(DSCR)/cartel_productivity.jl
	cd $(DSCR) && julia $(<F)
	touch $@

t_group_productivity_fig: $(DSCR)/group-productivity.jl t_group_productivity
	cd $(DSCR) && julia $(<F)
	touch $@

