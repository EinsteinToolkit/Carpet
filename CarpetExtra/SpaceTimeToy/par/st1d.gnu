# $Header: /home/eschnett/C/carpet/Carpet/CarpetExtra/SpaceTimeToy/par/st1d.gnu,v 1.1 2003/06/18 18:24:29 schnetter Exp $

set grid

dt=0.05

p [0:1] "st1d_1l_0020/phi.xl" i t u 10:13 w l, "st1d_1l_0040/phi.xl" i t u 10:13 w l, "st1d_1l_0080/phi.xl" i t u 10:13 w l, cos (2*pi*(x+dt*t))
p [0:1] "st1d_2l_0020/phi.xl" i 2*t u 10:13 w l, "st1d_2l_0040/phi.xl" i 2*t u 10:13 w l, "st1d_2l_0080/phi.xl" i 2*t u 10:13 w l, "st1d_2l_0160/phi.xl" i 2*t u 10:13 w l, "st1d_2l_0320/phi.xl" i 2*t u 10:13 w l, cos (2*pi*(x+dt*t))

p [0:1] "st1d_1l_0020/phi.xl" i t u 10:($13-cos(2*pi*($10+dt*t))) w l, "st1d_1l_0040/phi.xl" i t u 10:(4*($13-cos(2*pi*($10+dt*t)))) w l, "st1d_1l_0080/phi.xl" i t u 10:(16*($13-cos(2*pi*($10+dt*t)))) w l
p [0:1] "st1d_2l_0020/phi.xl" i 2*t u 10:($13-cos(2*pi*($10+dt*t))) w l, "st1d_2l_0040/phi.xl" i 2*t u 10:(4*($13-cos(2*pi*($10+dt*t)))) w l, "st1d_2l_0080/phi.xl" i 2*t u 10:(16*($13-cos(2*pi*($10+dt*t)))) w l, "st1d_2l_0160/phi.xl" i 2*t u 10:(64*($13-cos(2*pi*($10+dt*t)))) w l, "st1d_2l_0320/phi.xl" i 2*t u 10:(256*($13-cos(2*pi*($10+dt*t)))) w l
