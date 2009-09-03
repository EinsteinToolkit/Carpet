set terminal postscript eps enhanced color colortext
set output "results-best.eps"

set size 0.5

set logscale x
set xtics (1, 4, 16, 64, 256, "1k" 1024, "4k" 4096, "16k" 16384)
set xrange [0.7:24576]
set yrange [0:]

set title "Cactus Benchmark"
set xlabel "number of cores"
set ylabel "time per RHS evaluation [{/Symbol m}s]"

set key bottom

# 25^3 grid points per core
# 128 time steps
# 9 refinement levels
# a factor 2 because of the Berger-Oliger subcycling
# 4 Runge-Kutta substeps per time step
# 1e6 us per second
p \
"< awk '/^Bench_McLachlan_Carpet_9lev/ { if ($2==128 && $3==25                       && $7== 1 && $8== 0) print; }' results-franklin.out" u 5:($9/25**3/128/2/4*1e6) t "Franklin"         w lp lt 1 lw 3, \
"< awk '/^Bench_McLachlan_Carpet_9lev/ { if ($2==128 && $3==25                       && $7== 1 && $8==10) print; }' results-hlrb2.out"    u 5:($9/25**3/128/2/4*1e6) t "HLRB II"          w lp lt 2 lw 3, \
"< awk '/^Bench_McLachlan_Carpet_9lev/ { if ($2==128 && $3==25 && ($6== 4 || $6==$5) && $7== 1 && $8==50) print; }' results-kraken.out"   u 5:($9/25**3/128/2/4*1e6) t "Kraken (NT=1)"    w lp lt 3 lw 3, \
"< awk '/^Bench_McLachlan_Carpet_9lev/ { if ($2==128 && $3==25 && ($6== 8 || $6==$5) && $7== 2 && $8== 0) print; }' results-queenbee.out" u 5:($9/25**3/128/2/4*1e6) t "Queen Bee (NT=2)" w lp lt 4 lw 3, \
"< awk '/^Bench_McLachlan_Carpet_9lev/ { if ($2==128 && $3==25 && ($6==16 || $6==$5) && $7== 1 && $8==71) print; }' results-ranger.out"   u 5:($9/25**3/128/2/4*1e6) t "Ranger (NT=1)"    w lp lt 5 lw 3
