set terminal postscript eps enhanced color colortext
set output "results-old.eps"

set size 0.5

set logscale x
set xtics (1, 4, 16, 64, 256, "1k" 1024, "4k" 4096, "16k" 16384)
set xrange [0.7:4096]
set yrange [0:30]

set title "Cactus Benchmark"
set xlabel "number of cores"
set ylabel "time per RHS evaluation [{/Symbol m}s]"

set parametric
set trange [0:100]

set key bottom

# 25^3 grid points per core
# 128 time steps
# 9 refinement levels
# a factor 2 because of the Berger-Oliger subcycling
# 4 Runge-Kutta substeps per time step
# 1e6 us per second
p \
"KL9_CCT_vacuum.asc" u 1:($4/25**3/128/2/4*1e6) t "Ranger, Sept. 2008" w lp lw 3, \
16, t t "node boundary"
