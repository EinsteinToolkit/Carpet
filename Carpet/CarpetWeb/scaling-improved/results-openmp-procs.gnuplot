set terminal postscript eps enhanced color colortext
set output "results-openmp-procs.eps"

set size 0.5

set logscale x
set xtics 2
set xrange [0.7:24]
set yrange [0:]

set title "Cactus Benchmark (using 16 cores)"
set xlabel "number of allocated nodes"
set ylabel "time per RHS evaluation [{/Symbol m}s]"

set key bottom

# 25^3 grid points per core
# 128 time steps
# 9 refinement levels
# a factor 2 because of the Berger-Oliger subcycling
# 4 Runge-Kutta substeps per time step
# 1e6 us per second
p \
"results-ranger-openmp-procs.out" u ($5/$6):($9/25**3/128/2/4*1e6) t "Ranger (varying NT)" w lp lw 3
