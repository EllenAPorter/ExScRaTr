DOCUMENT = talk

PDF_PREVIEWER = evince
PDFLATEX_FLAGS =  -file-line-error -interaction=nonstopmode -halt-on-error

INCLUDED_FILES = \
	$(wildcard images/*.png) \
	$(wildcard diagrams/*.pdf)

default: $(DOCUMENT).pdf

clean:
	rm -rf auto/ *.aux *.log *.nav *.out *.snm *~ *.toc

immaculate: clean
	rm -f $(DOCUMENT).pdf

preview: $(DOCUMENT).pdf
	$(PDF_PREVIEWER) $<

$(DOCUMENT).pdf: $(DOCUMENT).tex $(INCLUDED_FILES)
	pdflatex $(PDFLATEX_FLAGS) $<
	pdflatex $(PDFLATEX_FLAGS) $<
