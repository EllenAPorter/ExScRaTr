DOCNAME = main

PRINT = lpr
PRINTOPTS =

INCLUDED_FILES = \
	drawings/*.pdf \
	$(wildcard sec*.tex)

PDF_PREVIEW = evince
PDFLATEX_FLAGS =  -file-line-error -interaction=nonstopmode -halt-on-error

.PHONY: default
default: $(DOCNAME).pdf

.PHONY: clean
clean:
	rm -rf auto/ \
	    *.aux *.bbl *.blg *.dlog *.dvi *.lof *.log *.lot *.prv *.toc \
	    _region_.tex $(GENFILES) *~

.PHONY: immaculate
immaculate: clean
	rm -f $(DOCNAME).pdf

.PHONY: preview
preview: $(DOCNAME).pdf
	$(PDF_PREVIEW) $<

.PHONY: print
print: $(DOCNAME).pdf
	$(PRINT) $(PRINTOPTS) $<

$(DOCNAME).pdf: $(DOCNAME).tex $(INCLUDED_FILES) bibliography.bib
	pdflatex $<
	bibtex $(subst .tex,,$<)
	pdflatex $(PDFLATEX_FLAGS) $<
	pdflatex $(PDFLATEX_FLAGS) $<

