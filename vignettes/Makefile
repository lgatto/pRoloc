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
	pRoloc-imputation.tex\
	pRoloc-ml.tex\
	pRoloc-tutorial.tex

all: 
	make tut

setvars:
ifeq (${R_HOME},)
R_HOME= $(shell R RHOME)
endif

tut:
	"$(R_HOME)/bin/R" CMD pdflatex pRoloc-tutorial.tex
	bibtex pRoloc-tutorial
	"$(R_HOME)/bin/R" CMD pdflatex pRoloc-tutorial.tex
	"$(R_HOME)/bin/R" CMD pdflatex pRoloc-tutorial.tex

## manual building
tut2: 
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


allclean: clean
	rm -f  pRoloc-tutorial.pdf
	rm -rf figure/
	rm -rf cache/
