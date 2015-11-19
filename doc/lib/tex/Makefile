DOCNAME = porter_gi_2016

PRINT = lpr
PRINTOPTS =

PDF_PREVIEW = evince
PDFLATEX_FLAGS =  -file-line-error -interaction=nonstopmode -halt-on-error

default: $(DOCNAME).pdf

clean:
	rm -rf auto/ \
	    *.aux *.bbl *.blg *.dlog *.dvi *.log *.prv \
	    _region_.tex $(GENFILES) *~

immaculate: clean
	# if we decide not to keep PDF's under git
	# rm -f $(DOCNAME).pdf

preview: $(DOCNAME).pdf
	$(PDF_PREVIEW) $<

print: $(DOCNAME).pdf
	$(PRINT) $(PRINTOPTS) $<

$(DOCNAME).pdf: $(DOCNAME).tex
	pdflatex $<
	bibtex $(subst .tex,,$<)
	pdflatex $(PDFLATEX_FLAGS) $<
	pdflatex $(PDFLATEX_FLAGS) $<
