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
	make tut
	make tl
	make ml


ml: pRoloc-ml.Rnw
	pdflatex pRoloc-ml.tex
	bibtex pRoloc-ml
	pdflatex pRoloc-ml.tex
	pdflatex pRoloc-ml.tex

tut: pRoloc-tutorial.Rnw
	pdflatex pRoloc-tutorial.tex
	bibtex pRoloc-tutorial
	pdflatex pRoloc-tutorial.tex
	pdflatex pRoloc-tutorial.tex

tl: pRoloc-transfer-learning.Rnw
	pdflatex pRoloc-transfer-learning.tex
	bibtex pRoloc-transfer-learning
	pdflatex pRoloc-transfer-learning.tex
	pdflatex pRoloc-transfer-learning.tex

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
