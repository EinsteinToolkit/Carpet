FILE=results-wavetoy-abe

all: $(FILE).pdf $(FILE).png

%.png: %.eps
	convert -rotate 90 $*.eps $*.png

%.pdf: %.eps
	epstopdf $*.eps

%.eps: %.gnuplot
	gnuplot < $*.gnuplot

clean:
	$(RM) $(FILE).eps $(FILE).pdf $(FILE).png

.PSEUDO: all clean
