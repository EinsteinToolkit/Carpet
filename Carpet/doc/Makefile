# $Header: /home/eschnett/C/carpet/Carpet/Carpet/doc/Makefile,v 1.1 2002/11/29 10:43:43 schnetter Exp $

NAME=documentation

all: ${NAME}.dvi ${NAME}.ps ${NAME}.pdf

${NAME}.dvi: carpet.bib

%.dvi: %.tex
	latex $*
	bibtex $*
	latex $*
	latex $*

%.ps: %.dvi
	dvips $*

%.pdf: %.dvi
	pdflatex $*
	thumbpdf $*
	pdflatex $*

clean:
	$(RM) ${NAME}.aux ${NAME}.bbl ${NAME}.blg ${NAME}.dvi ${NAME}.log ${NAME}.out ${NAME}.pdf ${NAME}.ps ${NAME}.tpt

.PSEUDO: all clean
