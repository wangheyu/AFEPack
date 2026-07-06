TEX = xelatex
BIB = bibtex
MAIN = afepack-manual

.PHONY: all clean distclean

all: $(MAIN).pdf

$(MAIN).pdf: $(MAIN).tex references.bib
	$(TEX) -interaction=nonstopmode $(MAIN).tex
	$(BIB) $(MAIN)
	$(TEX) -interaction=nonstopmode $(MAIN).tex
	$(TEX) -interaction=nonstopmode $(MAIN).tex

clean:
	rm -f $(MAIN).aux $(MAIN).bbl $(MAIN).blg $(MAIN).log \
	      $(MAIN).out $(MAIN).toc $(MAIN).lof $(MAIN).lot

distclean: clean
	rm -f $(MAIN).pdf
