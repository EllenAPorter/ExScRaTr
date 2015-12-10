DOCNAME = ascr

PRINT = lpr
PRINTOPTS =

PDF_PREVIEW = evince
PDFLATEX_FLAGS =  -file-line-error -interaction=nonstopmode -halt-on-error

default: $(DOCNAME).pdf

clean:
	rm -f *.bbl *.blg *.dvi *.aux *.log *.dlog $(GENFILES) *~

immaculate: clean
	rm -f $(DOCNAME).pdf

preview: $(DOCNAME).pdf
	$(PDF_PREVIEW) $<

print: $(DOCNAME).pdf
	$(PRINT) $(PRINTOPTS) $<

$(DOCNAME).pdf: $(DOCNAME).tex
# if there's a biblography, enable these:
#	pdflatex $<
#	bibtex $(subst .tex,,$<)
	pdflatex $(PDFLATEX_FLAGS) $<
	pdflatex $(PDFLATEX_FLAGS) $<
