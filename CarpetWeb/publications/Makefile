SECTIONS = background refereed book report thesis webpage

ALLSECTIONS = $(SECTIONS) internal

all: carpet-publications.pdf publications.html

carpet-publications.pdf: carpet-publications.tex carpet-publications.bib $(SECTIONS:%=citations-%.bibitems)
	rubber -d carpet-publications

publications.html: $(SECTIONS:%=citations-%.html-updated)

.PRECIOUS: $(SECTIONS:%=citations-%.tex)
citations-%.tex: citations.template
	sed -e 's/STATUS/$*/g' citations.template > $@

.PRECIOUS: $(SECTIONS:%=citations-%.bib)
citations-%.bib: carpet-publications.bib
	bibtool -q -- 'select.by.string{status "$*"}' -- 'sort{on}' -- 'sort.reverse{on}' -- 'sort.format{%N(receiveddate)}' < carpet-publications.bib > $@

# this produces also %.aux and %.bbl
.PRECIOUS: $(SECTIONS:%=citations-%.pdf)
citations-%.pdf: citations-%.tex citations-%.bib
	rubber -d citations-$*

citations-%.bibitems: citations-%.pdf
	cat citations-$*.bbl | grep -v thebibliography | sed -e 's/bibitem/blurb/' > $@

# http://www.spinellis.gr/sw/textproc/bib2xhtml/
# Sorted alphabetically:
# ./bib2xhtml -r -s empty -u $*.aux publications.html
# ./bib2xhtml -s unsort -u $*.aux publications.html
# ./bib2xhtml -s unsortlist -u $*.aux publications.html
citations-%.html-updated: citations-%.pdf
	./bib2xhtml -s unsortlist -u citations-$*.aux publications.html
	: > $@

EMPTY :=
SPACE := $(EMPTY) $(EMPTY)
check:
	bibtool -q -- 'select.non{status "$(subst $(SPACE),\|,$(ALLSECTIONS))"}' < carpet-publications.bib

clean:
	rubber -d --clean carpet-publications
	$(RM) $(SECTIONS:%=citations-%.*)

.PHONY: all check clean
