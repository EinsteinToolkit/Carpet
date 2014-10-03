set terminal postscript eps enhanced color colortext
set output "openmp.eps"

set size 0.5

set logscale x
set logscale y

set title "Cactus Benchmark"
set xlabel "number of MPI processes"
set ylabel "number of OpenMP threads"
set xrange [0.2:24]
set yrange [0.2:24]
set xtics 2
set ytics 2
set cbrange [0:]
set cblabel "time per RHS evaluation [{/Symbol m}s]"

set pm3d map
set pm3d
#set pm3d corners2color c1

# 25^3 grid points per core
# 128 time steps
# 9 refinement levels
# a factor 2 because of the Berger-Oliger subcycling
# 4 Runge-Kutta substeps per time step
# 1e6 us per second
sp "< awk '/^Bench_McLachlan_Carpet_9lev/ { if ($2==128 && $3==25 && ($6==$5 || $6==16) && $8==71) print; }' results-ranger.out | sort -n -k 7 | awk 'BEGIN { print 0, 0, 0, 0, 0.5, 0, 0.5, 0, 100; print 0, 0, 0, 0, 0.5, 0, 1, 0, 100; print 0, 0, 0, 0, 0.5, 0, 2, 0, 100; print 0, 0, 0, 0, 0.5, 0, 4, 0, 100; print 0, 0, 0, 0, 0.5, 0, 8, 0, 100; print 0, 0, 0, 0, 0.5, 0, 16, 0, 100; } { if ($7!=old7) { print \"\"; print 0, 0, 0, 0, 0.5, 0, old7, 0, 100; } old7=$7; print; }';" u ($5/$7):($7):($9/25**3/128/2/4*1e6) notitle

# TODO: show also ppn-used vs. threads, for mpiprocs=1
