FILE=results-ranger

all: $(FILE).pdf $(FILE).png

%.png: %.eps
	convert $*.eps $*.png

%.pdf: %.eps
	epstopdf $*.eps

%.eps: %.gnuplot
	gnuplot < $*.gnuplot

clean:
	$(RM) $(FILE).eps $(FILE).pdf $(FILE).png

.PSEUDO: all clean
