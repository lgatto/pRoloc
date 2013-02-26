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
	cp pRoloc-tutorial.Rnw.knitr pRoloc-tutorial.Rnw
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit2pdf('pRoloc-tutorial.Rnw');"
	bibtex pRoloc-tutorial
	"$(R_HOME)/bin/R" CMD pdflatex pRoloc-tutorial.tex
	"$(R_HOME)/bin/R" CMD pdflatex pRoloc-tutorial.tex

imp:
	cp pRoloc-imputation.Rnw.knitr pRoloc-imputation.Rnw
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit2pdf('pRoloc-imputation.Rnw');"
	bibtex pRoloc-imputation
	"$(R_HOME)/bin/R" CMD pdflatex pRoloc-imputation.tex
	"$(R_HOME)/bin/R" CMD pdflatex pRoloc-imputation.tex

ml:
	cp pRoloc-ml.Rnw.knitr pRoloc-ml.Rnw
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit2pdf('pRoloc-ml.Rnw');"
	bibtex pRoloc-ml
	"$(R_HOME)/bin/R" CMD pdflatex pRoloc-ml.tex
	"$(R_HOME)/bin/R" CMD pdflatex pRoloc-ml.tex


.PHONY: clean allclean 

clean:	
	rm -f $(LATEXFILES)
	rm -f .Rhistory
	rm -f *~
	rm -f */*~
	cp pRoloc-tutorial.Rnw.stub pRoloc-tutorial.Rnw
	## cp pRoloc-imputation.Rnw.stub pRoloc-imputation.Rnw
	## cp pRoloc-ml.Rnw.stub pRoloc-ml.Rnw

allclean: clean
	rm -f  pRoloc-tutorial.pdf
	rm -f  pRoloc-imputation.pdf
	rm -f  pRoloc-ml.pdf
	rm -rf figure/
	rm -rf cache/
