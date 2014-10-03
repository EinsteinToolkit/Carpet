all: results-old.pdf results.pdf results-best.pdf results-scaled.pdf results-sicortex.pdf results-surveyor.pdf results-openmp-node.pdf results-openmp-procs.pdf results-whatif.pdf

results.eps results-best.eps: results-franklin.out results-hlrb2.out results-kraken.out results-queenbee.out results-ranger.out
results-scaled.eps: results-franklin.out results-hlrb2.out results-kraken.out results-queenbee.out results-ranger.out results-sicortex.out results-surveyor.out
results-sicortex.eps: results-sicortex.out
results-surveyor.eps: results-surveyor.out
results-openmp-node.eps: results-ranger-openmp-node.out
results-openmp-procs.eps: results-ranger-openmp-procs.out
results-whatif.eps: results-whatif-orig.out results-whatif-bw.out results-whatif-bwlat.out

results-ranger-openmp-procs.out results-ranger-openmp-node.out: results2openmp results-ranger.out
	./results2openmp

results-whatif-orig.out results-whatif-bw.out results-whatif-bwlat.out: results2whatif results-ranger.out
	./results2whatif

%.eps: %.gnuplot
	gnuplot < $*.gnuplot

%.pdf: %.eps
	epstopdf $*.eps

rsync:
	-rsync -Paz franklin.nersc.gov:/project/projectdirs/m152/schnette/franklin/Calpha/results-franklin.out .
	-rsync -Paz -e 'ssh numrel02.cct.lsu.edu ssh' lu78kog@hlrb2.lrz-muenchen.de:Calpha/results-hlrb2.out .
	-rsync -Paz kraken-pwd.nics.utk.edu:work/Calpha/results-kraken.out .
	-rsync -Paz queenbee.loni.org:work/Calpha/results-queenbee.out .
	-rsync -Paz ranger.tacc.utexas.edu:work/Calpha/results-ranger.out .

clean:
	rm -f results-old.pdf results.pdf results-best.pdf results-scaled.pdf results-sicortex.pdf results-surveyor.pdf results-openmp-node.pdf results-openmp-procs.pdf results-whatif.pdf
	rm -f results-old.eps results.eps results-best.eps results-scaled.eps results-sicortex.eps results-surveyor.eps results-openmp-node.eps results-openmp-procs.eps results-whatif.eps

.PHONY: all rsync clean

.PRECIOUS: results-old.eps results.eps results-best.eps results-scaled.eps results-sicortex.eps results-surveyor.eps results-openmp-node.eps results-openmp-procs.eps results-whatif.eps
