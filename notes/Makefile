LATEX:=pdflatex
LATEXOPT:=--shell-escape
NONSTOP:=--interaction=nonstopmode

LATEXMK:=latexmk
LATEXMKOPT:=-pdf

MAIN:=note
SOURCES:=$(MAIN).tex Makefile

default: show

show: $(MAIN).pdf
	vupdf $<

$(MAIN).pdf: $(SOURCES)
	$(LATEXMK) $(LATEXMKOPT) -pdflatex="$(LATEX) $(LATEXOPT) $(NONSTOP) %O %S" $(MAIN)

clean:
	$(LATEXMK) -C $(MAIN)

.PHONY: show clean
