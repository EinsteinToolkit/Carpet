NAME=scheduling

all: ${NAME}.dvi ${NAME}.ps ${NAME}.ps.gz ${NAME}.pdf

${NAME}.dvi: carpet.bib

%.dvi: %.tex
	latex $*
	bibtex $*
	latex $*
	latex $*

%.ps: %.dvi
	dvips $* -o ${NAME}.ps

%.pdf: %.dvi
	pdflatex $*
	thumbpdf $*
	pdflatex $*

%.gz: %
	gzip --best -c $* > $*.gz

clean:
	$(RM) ${NAME}.aux ${NAME}.bbl ${NAME}.blg ${NAME}.dvi ${NAME}.log ${NAME}.out ${NAME}.pdf ${NAME}.ps ${NAME}.ps.gz ${NAME}.tpt

.PSEUDO: all clean
