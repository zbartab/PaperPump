
# contains recipies to assemble to PaperPump paper

SRC = PaperPump.md

PDFS=$(SRC:.md=.pdf)
HTML=$(SRC:.md=.html)
TEX=$(SRC:.md=.tex)
DOCX=$(SRC:.md=.docx)

html: $(HTML)
pdf: $(PDFS)
docx: $(DOCX)

$(DOCX): $(SRC) groups.png weighted_production.png
	pandoc -c ~/lib/markdown/pandoc.css --mathml -N --standalone \
		--self-contained --filter pandoc-citeproc -o $@ $<

$(PDFS): $(SRC)
	pandoc -o $@ $<

$(HTML): $(SRC) groups.png weighted_production.png
	pandoc -c ~/lib/markdown/pandoc.css --mathml -N --standalone --toc \
		--self-contained --filter pandoc-citeproc -o $@ $<

%.png: %.R
	Rscript $<


