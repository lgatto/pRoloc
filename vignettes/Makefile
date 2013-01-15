LATEXFILES = *.aux\
	*.bbl\
	*.blg\
	*.ilg\
	*.log\
	*.nlo\
	*.nls\
	*.toc\
	*.aux\
	Rplots.pdf\
	*.dvi\
	*.map\
	*.out\
	*.figlist\
	*.dep\
	*.dpth\
	*.makefile\
	pRoloc-tutorial.tex

all: 
	make pRoloc-tutorial

setvars:
ifeq (${R_HOME},)
R_HOME= $(shell R RHOME)
endif

pRoloc-tutorial:
	cp pRoloc-tutorial.Rnw.knitr pRoloc-tutorial.Rnw
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit2pdf('pRoloc-tutorial.Rnw');"
	bibtex pRoloc-tutorial
	"$(R_HOME)/bin/R" CMD pdflatex pRoloc-tutorial.tex
	"$(R_HOME)/bin/R" CMD pdflatex pRoloc-tutorial.tex

.PHONY: clean allclean 

clean:	
	rm -f $(LATEXFILES)
	rm -f .Rhistory
	rm -f *~
	rm -f */*~
	cp pRoloc-tutorial.Rnw.stub pRoloc-tutorial.Rnw

allclean: clean
	rm -f  pRoloc-tutorial.pdf
	rm -rf figure/
	rm -rf cache/
