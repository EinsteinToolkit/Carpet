set terminal postscript eps enhanced color colortext
set output "results-openmp-node.eps"

set size 0.5

set logscale x
set xtics 2
set xrange [0.7:24]
set yrange [0:]

set title "Cactus Benchmark (using 1 node)"
set xlabel "number of MPI processes"
set ylabel "time per RHS evaluation [{/Symbol m}s]"

set key bottom

# 25^3 grid points per core
# 128 time steps
# 9 refinement levels
# a factor 2 because of the Berger-Oliger subcycling
# 4 Runge-Kutta substeps per time step
# 1e6 us per second
p \
"results-ranger-openmp-node.out" u 5:($9/25**3/128/2/4*1e6) t "Ranger (varying NT)" w lp lw 3
#p \
#"results-ranger-openmp-node.out" u 5:($7== 1 ? $9/25**3/128/2/4*1e6 : 0/0) t "Ranger (NT=1)"  w lp lw 3, \
#"results-ranger-openmp-node.out" u 5:($7== 2 ? $9/25**3/128/2/4*1e6 : 0/0) t "Ranger (NT=2)"  w lp lw 3, \
#"results-ranger-openmp-node.out" u 5:($7== 4 ? $9/25**3/128/2/4*1e6 : 0/0) t "Ranger (NT=4)"  w lp lw 3, \
#"results-ranger-openmp-node.out" u 5:($7== 8 ? $9/25**3/128/2/4*1e6 : 0/0) t "Ranger (NT=8)"  w lp lw 3, \
#"results-ranger-openmp-node.out" u 5:($7==16 ? $9/25**3/128/2/4*1e6 : 0/0) t "Ranger (NT=16)" w lp lw 3
