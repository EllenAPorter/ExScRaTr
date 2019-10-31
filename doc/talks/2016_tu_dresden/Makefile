DOCUMENT = talk

PDF_PREVIEWER = evince
PDFLATEX_FLAGS =  -file-line-error -interaction=nonstopmode -halt-on-error


INCLUDED_FILES = \
	$(wildcard diagrams/*.pdf) \
	$(wildcard image/*.png) \


default: $(DOCUMENT).pdf

clean:
	rm -rf auto/ *.aux *.dvi *.log *.nav *.out *.snm *~ *.toc $(INCLUDED_FILES)

immaculate: clean
	rm -f $(DOCUMENT).pdf

preview: $(DOCUMENT).pdf
	$(PDF_PREVIEWER) $<

$(DOCUMENT).pdf: $(DOCUMENT).tex $(INCLUDED_FILES)
	pdflatex $(PDFLATEX_FLAGS) $<
	pdflatex $(PDFLATEX_FLAGS) $<

diagrams/%.png: diagrams/%.dia
	dia -t png -e $@ $<

diagrams/%.eps: diagrams/%.dia
	dia -t eps -e $@ $<

diagrams/%.pdf: diagrams/%.eps
	epstopdf $<
