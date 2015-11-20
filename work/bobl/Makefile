PS_FILES = $(wildcard *.ps)

PDF_FILES = $(subst .ps,.pdf,$(PS_FILES))

PS2PDF_ARGS = -dEPSCrop

default: $(PDF_FILES)

%.pdf: %.ps
	ps2pdf $(PS2PDF_ARGS) $<
