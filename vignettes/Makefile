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
	make tutorial

setvars:
ifeq (${R_HOME},)
R_HOME= $(shell R RHOME)
endif

tutorial:
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit2pdf('pRoloc-tutorial.Rnw');"
	bibtex pRoloc-tutorial
	"$(R_HOME)/bin/R" CMD pdflatex pRoloc-tutorial.tex
	"$(R_HOME)/bin/R" CMD pdflatex pRoloc-tutorial.tex

clean:	
	rm -f $(LATEXFILES)
	rm -f .Rhistory
	rm -f *~
	rm -f */*~
	# rm -rf figure/
	# rm -rf cache/

pdfclean:
	rm -r pRoloc.pdf
	rm tutorial.pdf
