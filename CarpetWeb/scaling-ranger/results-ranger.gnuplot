set terminal postscript enhanced colour eps
set output "results-ranger.eps"
set size 0.65

set title "Weak Scaling on Ranger"
set key left

set xlabel "cores"
set log x
set xrange [16:4096]
set xtics 16, 4, 4096

set ylabel "time per grid point [ms]"
set yrange [0:0.05]

second = 1000   # 1 second has 1000 milliseconds
numrhs = 4      # 4 rhs evaluations per time step
time_1lev(steps,size,walltime) = walltime / (numrhs *    steps  * size**3) * second
time_9lev(steps,size,walltime) = walltime / (numrhs * (2*steps) * size**3) * second

# plot \
# "< ./loopprocs ./findrow benchmark=Bench_Ccatie_PUGH           machine=ranger tpp=1 results-ranger.out" u 5:(time_1lev($2,$3,$9)) t "Ccatie, PUGH, NT=1"               w lp lt 1 lw 3, \
# "< ./loopprocs ./findrow benchmark=Bench_McLachlan_PUGH        machine=ranger tpp=1 results-ranger.out" u 5:(time_1lev($2,$3,$9)) t "McLachlan, PUGH, NT=1"            w lp lt 2 lw 3, \
# "< ./loopprocs ./findrow benchmark=Bench_McLachlan_PUGH        machine=ranger tpp=4 results-ranger.out" u 5:(time_1lev($2,$3,$9)) t "McLachlan, PUGH, NT=4"            w lp lt 3 lw 3, \
# "< ./loopprocs ./findrow benchmark=Bench_Ccatie_Carpet_1lev    machine=ranger tpp=1 results-ranger.out" u 5:(time_1lev($2,$3,$9)) t "Ccatie, Carpet unigrid, NT=1"     w lp lt 4 lw 3, \
# "< ./loopprocs ./findrow benchmark=Bench_McLachlan_Carpet_1lev machine=ranger tpp=1 results-ranger.out" u 5:(time_1lev($2,$3,$9)) t "McLachlan, Carpet unigrid, NT=1"  w lp lt 5 lw 3, \
# "< ./loopprocs ./findrow benchmark=Bench_McLachlan_Carpet_1lev machine=ranger tpp=4 results-ranger.out" u 5:(time_1lev($2,$3,$9)) t "McLachlan, Carpet unigrid, NT=4"  w lp lt 6 lw 3, \
# "< ./loopprocs ./findrow benchmark=Bench_Ccatie_Carpet_9lev    machine=ranger tpp=1 results-ranger.out" u 5:(time_9lev($2,$3,$9)) t "Ccatie, Carpet 9 levels, NT=1"    w lp lt 7 lw 3, \
# "< ./loopprocs ./findrow benchmark=Bench_McLachlan_Carpet_9lev machine=ranger tpp=1 results-ranger.out" u 5:(time_9lev($2,$3,$9)) t "McLachlan, Carpet 9 levels, NT=1" w lp lt 8 lw 3, \
# "< ./loopprocs ./findrow benchmark=Bench_McLachlan_Carpet_9lev machine=ranger tpp=4 results-ranger.out" u 5:(time_9lev($2,$3,$9)) t "McLachlan, Carpet 9 levels, NT=4" w lp lt 9 lw 3

plot \
"< ./loopprocs ./findrow benchmark=Bench_McLachlan_PUGH        machine=ranger tpp=4 results-ranger.out" u 5:(time_1lev($2,$3,$9)) t "McLachlan, PUGH, NT=4"            w lp lt 1 lw 3, \
"< ./loopprocs ./findrow benchmark=Bench_McLachlan_Carpet_1lev machine=ranger tpp=4 results-ranger.out" u 5:(time_1lev($2,$3,$9)) t "McLachlan, Carpet unigrid, NT=4"  w lp lt 2 lw 3, \
"< ./loopprocs ./findrow benchmark=Bench_McLachlan_Carpet_9lev machine=ranger tpp=4 results-ranger.out" u 5:(time_9lev($2,$3,$9)) t "McLachlan, Carpet 9 levels, NT=4" w lp lt 3 lw 3
