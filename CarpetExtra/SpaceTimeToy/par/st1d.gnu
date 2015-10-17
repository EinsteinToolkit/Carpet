# $Header:$

set grid

dt=0.05
f(x)=cos(2*pi*(x+t*dt))

t=200


p [0:1] "st1d_1l_0020/phi.xl" i t u 10:13 w l, "st1d_1l_0040/phi.xl" i t u 10:13 w l, "st1d_1l_0080/phi.xl" i t u 10:13 w l, f(x)

p [0:1] "st1d_1l_0020/phi.xl" i t u 10:($13-f($10)) w l, "st1d_1l_0040/phi.xl" i t u 10:(4*($13-f($10))) w l, "st1d_1l_0080/phi.xl" i t u 10:(16*($13-f($10))) w l



p [0:1] "st1d_2l_0020/phi.xl" i 2*t u 10:13 w l, "st1d_2l_0040/phi.xl" i 2*t u 10:13 w l, "st1d_2l_0080/phi.xl" i 2*t u 10:13 w l, "st1d_2l_0160/phi.xl" i 2*t u 10:13 w l, "st1d_2l_0320/phi.xl" i 2*t u 10:13 w l, f(x)

p [0:1] "st1d_2l_0020/phi.xl" i 2*t u 10:($13-f($10)) w l, "st1d_2l_0040/phi.xl" i 2*t u 10:(4*($13-f($10))) w l, "st1d_2l_0080/phi.xl" i 2*t u 10:(16*($13-f($10))) w l, "st1d_2l_0160/phi.xl" i 2*t u 10:(64*($13-f($10))) w l, "st1d_2l_0320/phi.xl" i 2*t u 10:(256*($13-f($10))) w l



set grid

dt=0.05

p [0:1] "st1d_1l_0020/phi.xl" i t u 10:13 w l, "st1d_1l_0040/phi.xl" i t u 10:13 w l, "st1d_1l_0080/phi.xl" i t u 10:13 w l, cos (2*pi*(x+dt*t))
p [0:1] "st1d_2l_0020/phi.xl" i 2*t u 10:13 w l, "st1d_2l_0040/phi.xl" i 2*t u 10:13 w l, "st1d_2l_0080/phi.xl" i 2*t u 10:13 w l, "st1d_2l_0160/phi.xl" i 2*t u 10:13 w l, "st1d_2l_0320/phi.xl" i 2*t u 10:13 w l, cos (2*pi*(x+dt*t))

p [0:1] "st1d_1l_0020/phi.xl" i t u 10:($13-cos(2*pi*($10+dt*t))) w l, "st1d_1l_0040/phi.xl" i t u 10:(4*($13-cos(2*pi*($10+dt*t)))) w l, "st1d_1l_0080/phi.xl" i t u 10:(16*($13-cos(2*pi*($10+dt*t)))) w l
p [0:1] "st1d_2l_0020/phi.xl" i 2*t u 10:($13-cos(2*pi*($10+dt*t))) w l, "st1d_2l_0040/phi.xl" i 2*t u 10:(4*($13-cos(2*pi*($10+dt*t)))) w l, "st1d_2l_0080/phi.xl" i 2*t u 10:(16*($13-cos(2*pi*($10+dt*t)))) w l, "st1d_2l_0160/phi.xl" i 2*t u 10:(64*($13-cos(2*pi*($10+dt*t)))) w l, "st1d_2l_0320/phi.xl" i 2*t u 10:(256*($13-cos(2*pi*($10+dt*t)))) w l
