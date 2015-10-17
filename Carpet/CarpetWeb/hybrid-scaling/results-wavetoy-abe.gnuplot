#set terminal postscript color enhanced eps
set terminal postscript color enhanced
set size 0.55

set output "results-wavetoy-abe.eps"

set title "100^3 WaveToy on Abe"

set key bottom
set xlabel "processors"
set ylabel "wall time / grid point [{/Symbol m}s]"

set log x
set xrange [1:8192]
set xtics 1, 4
set yrange [0:0.1]

gridpoints=100**3
time(x)=x*1.0e+6/gridpoints

p \
"< < results-wavetoy-abe.out grep 'Bench_WaveToy_PUGH_100l'        | grep ' 22 '" u 3:(time($5/1000.0)) t "PUGH, MPI only"              w lp lt 1 lw 3, \
"< < results-wavetoy-abe.out grep 'Bench_WaveToy_Carpet_1lev_100l' | grep ' 22 '" u 3:(time($5/1023.0)) t "Carpet, 1 level, MPI only"   w lp lt 2 lw 3, \
"< < results-wavetoy-abe.out grep 'Bench_WaveToy_Carpet_1lev_100l' | grep ' 24 '" u 3:(time($5/1023.0)) t "Carpet, 1 level, MPI+OpenMP" w lp lt 3 lw 3
