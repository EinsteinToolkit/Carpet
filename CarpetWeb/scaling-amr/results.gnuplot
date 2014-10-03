set size 0.5
set terminal postscript enhanced eps color

set log x
set xlabel "cores"
set ylabel "time per grid point [{/Symbol m}s]"

set key bottom right
set title "McLachlan/PUGH Scaling"
set output "results-pugh.eps"
p [][0:20] \
"< < results-franklin.out awk '/Bench_McLachlan_PUGH/ { if ($3==65 && $6== 2 && $7==1) print; }' | ./minimise.sh" u 5:(1.0e6*$9/(4*$2*$3**3)) t "Franklin"  w lp lw 3, \
"< < results-queenbee.out awk '/Bench_McLachlan_PUGH/ { if ($3==65 && $6== 8 && $7==8) print; }' | ./minimise.sh" u 5:(1.0e6*$9/(4*$2*$3**3)) t "Queen Bee" w lp lw 3, \
"< < results-ranger.out   awk '/Bench_McLachlan_PUGH/ { if ($3==65 && $6==16 && $7==4) print; }' | ./minimise.sh" u 5:(1.0e6*$9/(4*$2*$3**3)) t "Ranger"    w lp lw 3
set output
!epstopdf results-pugh.eps

set key bottom right
set title "McLachlan/Carpet unigrid Scaling"
set output "results-carpet-1lev.eps"
p [][0:20] \
"< < results-franklin.out awk '/Bench_McLachlan_Carpet_1lev/ { if ($3==65 && $6== 2 && $7==1) print; }' | ./minimise.sh" u 5:(1.0e6*$9/(4*$2*$3**3)) t "Franklin"  w lp lw 3, \
"< < results-queenbee.out awk '/Bench_McLachlan_Carpet_1lev/ { if ($3==65 && $6== 8 && $7==8) print; }' | ./minimise.sh" u 5:(1.0e6*$9/(4*$2*$3**3)) t "Queen Bee" w lp lw 3, \
"< < results-ranger.out   awk '/Bench_McLachlan_Carpet_1lev/ { if ($3==65 && $6==16 && $7==4) print; }' | ./minimise.sh" u 5:(1.0e6*$9/(4*$2*$3**3)) t "Ranger"    w lp lw 3
!epstopdf results-carpet-1lev.eps

set key top left
set title "McLachlan/Carpet AMR Scaling"
set output "results-carpet-9lev.eps"
p [][0:200] \
"< < results-franklin.out awk '/Bench_McLachlan_Carpet_9lev/ { if ($3==25 && $6== 2 && $7==1) print; }' | ./minimise.sh" u 5:(1.0e6*$9/(4*$2*$3**3)) t "Franklin"  w lp lw 3, \
"< < results-queenbee.out awk '/Bench_McLachlan_Carpet_9lev/ { if ($3==25 && $6== 8 && $7==8) print; }' | ./minimise.sh" u 5:(1.0e6*$9/(4*$2*$3**3)) t "Queen Bee" w lp lw 3, \
"< < results-ranger.out   awk '/Bench_McLachlan_Carpet_9lev/ { if ($3==25 && $6==16 && $7==4) print; }' | ./minimise.sh" u 5:(1.0e6*$9/(4*$2*$3**3)) t "Ranger"    w lp lw 3
!epstopdf results-carpet-9lev.eps
