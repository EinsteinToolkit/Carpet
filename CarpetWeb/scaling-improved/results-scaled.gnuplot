set terminal postscript eps enhanced color colortext
set output "results-scaled.eps"

set size 0.5

set logscale x
set xtics (1, 4, 16, 64, 256, "1k" 1024, "4k" 4096, "16k" 16384)
set xrange [0.7:24576]
set yrange [0:]

set title "Cactus Benchmark"
set xlabel "number of cores"
set ylabel "1000 cycles per RHS evaluation"

set key bottom

# 25^3 grid points per core
# 128 time steps
# 9 refinement levels
# a factor 2 because of the Berger-Oliger subcycling
# 4 Runge-Kutta substeps per time step
# 1e6 us per second

#             MHz flop/cycle
# Franklin : 2300 4
# HLRB II  : 1600 4
# Kraken   : 2300 4
# Queen Bee: 2330 4
# Ranger   : 2300 4
# SiCortex :  700 2
# Surveyor :  850 4

p \
"< awk '/^Bench_McLachlan_Carpet_9lev/ { if ($2==128 && $3==25                       && $7== 1 && $8== 0) print; }' results-franklin.out" u 5:($9/25**3/128/2/4*1e3*(2300*4)) t "Franklin"         w lp lt 1 lw 3, \
"< awk '/^Bench_McLachlan_Carpet_9lev/ { if ($2==128 && $3==25                       && $7== 1 && $8==10) print; }' results-hlrb2.out"    u 5:($9/25**3/128/2/4*1e3*(1600*4)) t "HLRB II"          w lp lt 2 lw 3, \
"< awk '/^Bench_McLachlan_Carpet_9lev/ { if ($2==128 && $3==25 && ($6== 4 || $6==$5) && $7== 1 && $8==50) print; }' results-kraken.out"   u 5:($9/25**3/128/2/4*1e3*(2300*4)) t "Kraken (NT=1)"    w lp lt 3 lw 3, \
"< awk '/^Bench_McLachlan_Carpet_9lev/ { if ($2==128 && $3==25 && ($6== 8 || $6==$5) && $7== 2 && $8== 0) print; }' results-queenbee.out" u 5:($9/25**3/128/2/4*1e3*(2330*4)) t "Queen Bee (NT=2)" w lp lt 4 lw 3, \
"< awk '/^Bench_McLachlan_Carpet_9lev/ { if ($2==128 && $3==25 && ($6==16 || $6==$5) && $7== 1 && $8==71) print; }' results-ranger.out"   u 5:($9/25**3/128/2/4*1e3*(2300*4)) t "Ranger (NT=1)"    w lp lt 5 lw 3, \
"< awk '/^Bench_McLachlan_Carpet_9lev/ { if ($2==128 && $3==20 && ($6== 6 || $6==$5) && $7== 1 && $8== 9) print; }' results-sicortex.out" u 5:($9/20**3/128/2/4*1e3*( 700*2)) t "SiCortex (NT=1)"  w lp lt 6 lw 3, \
"< awk '/^Bench_McLachlan_Carpet_9lev/ { if ($2==128 && $3==13 && ($6== 4 || $6==$5) && $7== 1 && $8==12) print; }' results-surveyor.out" u 5:($9/13**3/128/2/4*1e3*( 850*4)) t "Surveyor (NT=1)"  w lp lt 7 lw 3
