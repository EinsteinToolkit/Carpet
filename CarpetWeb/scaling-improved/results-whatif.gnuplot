set terminal postscript eps enhanced color colortext
set output "results-whatif.eps"

set size 0.5

set logscale x
set xtics (1, 4, 16, 64, 256, "1k" 1024, "4k" 4096, "16k" 16384)
set xrange [0.7:24576]
set yrange [0:]

set title "Cactus Benchmark (Kraken, NT=1)"
set xlabel "number of cores"
#set ylabel "time per RHS evaluation [{/Symbol m}s]"
set ylabel "wall time increase"

set key bottom

# 25^3 grid points per core
# 128 time steps
# 9 refinement levels
# a factor 2 because of the Berger-Oliger subcycling
# 4 Runge-Kutta substeps per time step
# 1e6 us per second
#p \
#"< awk '/^Bench_McLachlan_Carpet_9lev/ { if ($2==128 && $3==25 && ($6== 4 || $6==$5) && $7== 1 && $8==50) print; }' results-kraken.out"   u 5:($9/25**3/128/2/4*1e6) t "original"          w lp lt 1 lw 3, \
#"< awk '/^Bench_McLachlan_Carpet_9lev/ { if ($2==128 && $3==25 && ($6== 4 || $6==$5) && $7== 1 && $8==51) print; }' results-kraken.out"   u 5:($9/25**3/128/2/4*1e6) t "bandwidth"         w lp lt 2 lw 3, \
#"< awk '/^Bench_McLachlan_Carpet_9lev/ { if ($2==128 && $3==25 && ($6== 4 || $6==$5) && $7== 1 && $8==52) print; }' results-kraken.out"   u 5:($9/25**3/128/2/4*1e6) t "bandwidth+latency" w lp lt 3 lw 3

p \
"< paste results-whatif-orig.out results-whatif-orig.out"  u 5:($19/$9) t "original"             w lp lt 1 lw 3, \
"< paste results-whatif-orig.out results-whatif-bw.out"    u 5:($19/$9) t "increased BW"         w lp lt 2 lw 3, \
"< paste results-whatif-orig.out results-whatif-bwlat.out" u 5:($19/$9) t "increased BW and LAT" w lp lt 3 lw 3
