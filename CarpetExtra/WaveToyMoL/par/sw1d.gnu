# $Header: /home/eschnett/C/carpet/Carpet/CarpetExtra/WaveToyMoL/par/sw1d.gnu,v 1.2 2003/06/30 17:28:51 schnetter Exp $

set grid

dt=0.05
f(x)=cos(2*pi*(x+t*dt))

t=200


p [0:1] "sw1d_1l_0020/phi.xl" i t u 10:13 w l, "sw1d_1l_0040/phi.xl" i t u 10:13 w l, "sw1d_1l_0080/phi.xl" i t u 10:13 w l, f(x)

p [0:1] "sw1d_1l_0020/phi.xl" i t u 10:($13-f($10)) w l, "sw1d_1l_0040/phi.xl" i t u 10:(4*($13-f($10))) w l, "sw1d_1l_0080/phi.xl" i t u 10:(16*($13-f($10))) w l



p [0:1] "sw1d_2l_0020/phi.xl" i 2*t u 10:13 w l, "sw1d_2l_0040/phi.xl" i 2*t u 10:13 w l, "sw1d_2l_0080/phi.xl" i 2*t u 10:13 w l, "sw1d_2l_0160/phi.xl" i 2*t u 10:13 w l, "sw1d_2l_0320/phi.xl" i 2*t u 10:13 w l, f(x)

p [0:1] "sw1d_2l_0020/phi.xl" i 2*t u 10:($13-f($10)) w l, "sw1d_2l_0040/phi.xl" i 2*t u 10:(4*($13-f($10))) w l, "sw1d_2l_0080/phi.xl" i 2*t u 10:(16*($13-f($10))) w l, "sw1d_2l_0160/phi.xl" i 2*t u 10:(64*($13-f($10))) w l, "sw1d_2l_0320/phi.xl" i 2*t u 10:(256*($13-f($10))) w l





set grid


dt=0.1
A=1.0/16.0

f(x)=A*t*dt*sin(2*pi*x)*cos(2*pi*dt*t)



t=400

p [0:1] "sw1d_1l_0020/phierror.xl" i t u 10:13 w lp, "sw1d_1l_0040/phierror.xl" i 2*t u 10:(4*$13) w lp, "sw1d_1l_0080/phierror.xl" i 4*t u 10:(16*$13) w lp



p [0:1] "sw1d_2l_0020/phi.xl" i 3*t ev :::2::2 u 10:13 w l, "sw1d_2l_0040/phi.xl" i 3*2*t ev :::2::2 u 10:13 w l, "sw1d_2l_0080/phi.xl" i 3*4*t ev :::2::2 u 10:13 w l, "sw1d_2l_0160/phi.xl" i 3*8*t ev :::2::2 u 10:13 w l, "sw1d_2l_0320/phi.xl" i 3*16*t ev :::2::2 u 10:13 w l, cos(2*pi*x)
p [0:1] "sw1d_2l_20/phierror.xl" i 3*t u 10:13 w l, "sw1d_2l_40/phierror.xl" i 3*2*t u 10:(4*$13) w l, "sw1d_2l_80/phierror.xl" i 3*4*t u 10:(16*$13) w l, "sw1d_2l_160/phierror.xl" i 3*8*t u 10:(64*$13) w l, "sw1d_2l_320/phierror.xl" i 3*16*t u 10:(256*$13) w l

p [0:1] "sw1d_2l_0020/phi.xl" i 2*t u 10:13 w l, "sw1d_2l_0040/phi.xl" i 2*t u 10:13 w l, "sw1d_2l_0080/phi.xl" i 2*t u 10:13 w l, "sw1d_2l_0160/phi.xl" i 2*t u 10:13 w l, "sw1d_2l_0320/phi.xl" i 2*t u 10:13 w l, cos(2*pi*x)
p [0:1] "sw1d_2l_0020/phierror.xl" i 2*t u 10:13 w l, "sw1d_2l_0040/phierror.xl" i 2*t u 10:(4*$13) w l, "sw1d_2l_0080/phierror.xl" i 2*t u 10:(16*$13) w l, "sw1d_2l_0160/phierror.xl" i 2*t u 10:(64*$13) w l, "sw1d_2l_0320/phierror.xl" i 2*t u 10:(256*$13) w l
p [0:1] "sw1d_2l_0020/phierror.xl" i 2*t u 10:($13-f($10)/20) w l, "sw1d_2l_0040/phierror.xl" i 2*t u 10:(4*($13-f($10)/40)) w l, "sw1d_2l_0080/phierror.xl" i 2*t u 10:(16*($13-f($10)/80)) w l, "sw1d_2l_0160/phierror.xl" i 2*t u 10:(64*($13-f($10)/160)) w l, "sw1d_2l_0320/phierror.xl" i 2*t u 10:(256*($13-f($10)/320)) w l



p [:] "< tail +6 sw1d_1l_0020/phierror_norm2.xg" t "unigrid error, dx=1/20" w l, "< tail +6 sw1d_1l_0040/phierror_norm2.xg" u 1:(4*$2) t "4 * unigrid error, dx=1/40" w l, "< tail +6 sw1d_2l_0020/phierror_norm2.xg" t "FMR error, dx=1/20" w l, "< tail +6 sw1d_2l_0040/phierror_norm2.xg" u 1:(4*$2) t "4 * FMR error, dx=1/40" w l, "< tail +6 sw1d_2l_0080/phierror_norm2.xg" u 1:(16*$2) t "16 * FMR error, dx=1/80" w l, "< tail +6 sw1d_2l_0160/phierror_norm2.xg" u 1:(64*$2) t "64 * FMR error, dx=1/160" w l, "< tail +6 sw1d_2l_0320/phierror_norm2.xg" u 1:(256*$2) t "256 * FMR error, dx=1/320" w l
, "< tail +6 sw1d_1l_0080/phierror_norm2.xg" u 1:(16*$2) t "16 * unigrid error, dx=1/80" w l

p [:5] "< tail +6 sw1d_1l_20/phierror_norm2.xg" t "unigrid error, dx=1/20" w l, "< tail +6 sw1d_1l_40/phierror_norm2.xg" u 1:(4*$2) t "4 * unigrid error, dx=1/40" w l, "< tail +6 sw1d_1l_80/phierror_norm2.xg" u 1:(16*$2) t "16 * unigrid error, dx=1/80" w l, "< tail +6 sw1d_2l_20/phierror_norm2.xg" ev 2 t "FMR error, dx=1/20" w l, "< tail +6 sw1d_2l_40/phierror_norm2.xg" ev 2 u 1:(4*$2) t "4 * FMR error, dx=1/40" w l, "< tail +6 sw1d_2l_80/phierror_norm2.xg" ev 2 u 1:(16*$2) t "16 * FMR error, dx=1/80" w l, "< tail +6 sw1d_2l_160/phierror_norm2.xg" ev 2 u 1:(64*$2) t "64 * FMR error, dx=1/160" w l, "< tail +6 sw1d_2l_320/phierror_norm2.xg" ev 2 u 1:(256*$2) t "256 * FMR error, dx=1/320" w l

p [:5] "< tail +6 sw1d_1l_20/psierror_norm2.xg" w l, "< tail +6 sw1d_1l_40/psierror_norm2.xg" u 1:(4*$2) w l, "< tail +6 sw1d_1l_80/psierror_norm2.xg" u 1:(16*$2) w l, "< tail +6 sw1d_2l_20/psierror_norm2.xg" ev 2 w l, "< tail +6 sw1d_2l_40/psierror_norm2.xg" ev 2 u 1:(4*$2) w l, "< tail +6 sw1d_2l_80/psierror_norm2.xg" ev 2 u 1:(16*$2) w l, "< tail +6 sw1d_2l_160/psierror_norm2.xg" ev 2 u 1:(64*$2) w l
