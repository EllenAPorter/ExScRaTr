DOCNAME = main

PRINT = lpr
PRINTOPTS =

INCLUDED_FILES = \
	drawings/CnC.pdf \
	drawings/Lights1.pdf \
	drawings/Lights2.pdf \
	drawings/Trace1.pdf \
	drawings/Trace2.pdf \
	drawings/Trace3.pdf \
	drawings/Trace4.pdf \
	$(wildcard sec*.tex)

PDF_PREVIEW = evince
PDFLATEX_FLAGS =  -file-line-error -interaction=nonstopmode -halt-on-error

.PHONY: default
default: $(DOCNAME).pdf

.PHONY: clean
clean:
	rm -rf auto/ \
	    *.toc *.lot *.lof *.aux *.bbl *.blg *.dlog *.dvi *.log *.prv \
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

