set terminal postscript eps enhanced color colortext
set output "results-sicortex.eps"

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
"< awk '/^Bench_McLachlan_Carpet_9lev/ { if ($2==128 && $3==20 && ($6== 6 || $6==$5) && $7== 1 && $8== 9) print; }' results-sicortex.out" u 5:($9/20**3/128/2/4*1e6) t "SiCortex (NT=1)" w lp lt 1 lw 3
