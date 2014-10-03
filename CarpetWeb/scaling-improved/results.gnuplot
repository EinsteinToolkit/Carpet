set terminal postscript eps enhanced color colortext
set output "results.eps"

#set size 0.5

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
"< awk '/^Bench_McLachlan_Carpet_9lev/ { if ($2==128 && $3==25                       && $7== 1 && $8==10) print; }' results-hlrb2.out"    u 5:($9/25**3/128/2/4*1e6) t "HLRB II"          w lp lt 1 lw 3, \
"< awk '/^Bench_McLachlan_Carpet_9lev/ { if ($2==128 && $3==25 && ($6== 4 || $6==$5) && $7== 1 && $8==50) print; }' results-kraken.out"   u 5:($9/25**3/128/2/4*1e6) t "Kraken (NT=1)"    w lp lt 2 lw 3, \
"< awk '/^Bench_McLachlan_Carpet_9lev/ { if ($2==128 && $3==25 && ($6== 8 || $6==$5) && $7== 1 && $8== 0) print; }' results-queenbee.out" u 5:($9/25**3/128/2/4*1e6) t "Queen Bee (NT=1)" w lp lt 3 lw 3, \
"< awk '/^Bench_McLachlan_Carpet_9lev/ { if ($2==128 && $3==25 && ($6== 8 || $6==$5) && $7== 2 && $8== 0) print; }' results-queenbee.out" u 5:($9/25**3/128/2/4*1e6) t "Queen Bee (NT=2)" w lp lt 3 lw 3, \
"< awk '/^Bench_McLachlan_Carpet_9lev/ { if ($2==128 && $3==25 && ($6== 8 || $6==$5) && $7== 4 && $8== 0) print; }' results-queenbee.out" u 5:($9/25**3/128/2/4*1e6) t "Queen Bee (NT=4)" w lp lt 3 lw 3, \
"< awk '/^Bench_McLachlan_Carpet_9lev/ { if ($2==128 && $3==25 && ($6== 8 || $6==$5) && $7== 8 && $8== 0) print; }' results-queenbee.out" u 5:($9/25**3/128/2/4*1e6) t "Queen Bee (NT=8)" w lp lt 3 lw 3, \
"< awk '/^Bench_McLachlan_Carpet_9lev/ { if ($2==128 && $3==25 && ($6==16 || $6==$5) && $7== 1 && $8==71) print; }' results-ranger.out"   u 5:($9/25**3/128/2/4*1e6) t "Ranger (NT=1)"    w lp lt 4 lw 3, \
"< awk '/^Bench_McLachlan_Carpet_9lev/ { if ($2==128 && $3==25 && ($6==16 || $6==$5) && $7== 2 && $8==71) print; }' results-ranger.out"   u 5:($9/25**3/128/2/4*1e6) t "Ranger (NT=2)"    w lp lt 4 lw 3, \
"< awk '/^Bench_McLachlan_Carpet_9lev/ { if ($2==128 && $3==25 && ($6==16 || $6==$5) && $7== 4 && $8==71) print; }' results-ranger.out"   u 5:($9/25**3/128/2/4*1e6) t "Ranger (NT=4)"    w lp lt 4 lw 3

# HLRB 2: have results up to 2048 cores; runs with 4080 cores failed,
# reason unknown.  Job with updated executable and fewer output files
# submitted; still fails for unknown reason.  Could be out of memory.

# Intrepid: no allocation.

# Kraken: have results up to 8192 cores; runs with up 18048 cores did
# fail, reason unknown.  (out of memory?  did work with NT=2 and
# NT=4.)

# Queen Bee: have results up to 1024 cores; runs with up to 2048 cores
# are submitted.

# Ranger: have results up to 12288 cores; runs with more cores depend
# on permissions from system administrators.
