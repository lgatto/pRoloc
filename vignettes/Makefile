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
	pRoloc-tutorial.tex\
	HUPO_2011_poster.tex

all: 
	make poster
	make ml
	make tut

setvars:
ifeq (${R_HOME},)
R_HOME= $(shell R RHOME)
endif


ml: pRoloc-ml.Rnw
	"$(R_HOME)/bin/R" CMD Sweave --engine=knitr::knitr --pdf pRoloc-ml.Rnw	

tut: pRoloc-tutorial.Rnw
	"$(R_HOME)/bin/R" CMD Sweave --engine=knitr::knitr --pdf pRoloc-tutorial.Rnw

poster: 
	cp vigsrc/HUPO_2011_poster.pdf .
	cp vigsrc/HUPO_2014_poster.pdf .

.PHONY: clean allclean poster

clean:	
	rm -f $(LATEXFILES)
	rm -f .Rhistory
	rm -f pRoloc-ml.R
	rm -f pRoloc-tutorial.R
	rm -f *~
	rm -f */*~


allclean: clean
	rm -f pRoloc-tutorial.pdf
	rm -f pRoloc-ml.pdf
	rm -rf figure/
	rm -rf cache/
