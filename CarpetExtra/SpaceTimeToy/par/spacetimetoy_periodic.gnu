# $Header:$

set grid
set style data linespoints

kx=0.222222222222222
ky=0.151515151515152
kz=0.222222222222222
omega=sqrt(kx**2+ky**2+kz**2)
dt=0.3
phi(x)=cos(pi*(kx*x+omega*dt*t))
psi(x)=-pi*omega*sin(pi*(kx*x+omega*dt*t))

t=0



# unigrid function
p [-4.5:4.5] "spacetimetoy_periodic_coarse/phi.xl" i t u 10:13, "spacetimetoy_periodic/phi.xl" i 2*t u 10:13, "spacetimetoy_periodic_fine/phi.xl" i 4*t u 10:13, phi(x)

# unigrid error
p [-4.5:4.5] "spacetimetoy_periodic_coarse/phi.xl" i t u 10:($13-phi($10)), "spacetimetoy_periodic/phi.xl" i 2*t u 10:(4*($13-phi($10))), "spacetimetoy_periodic_fine/phi.xl" i 4*t u 10:(16*($13-phi($10)))



# refinement function
p [-4.5:4.5] "spacetimetoy_periodic_coarse_rl2/phi.xl" i 3*t u 10:13, "spacetimetoy_periodic_rl2/phi.xl" i 6*t u 10:13, "spacetimetoy_periodic_fine_rl2/phi.xl" i 12*t u 10:13, phi(x)

# refinement error
p [-4.5:4.5] "spacetimetoy_periodic_coarse_rl2/phi.xl" i 3*t u 10:($13-phi($10)), "spacetimetoy_periodic_rl2/phi.xl" i 6*t u 10:(4*($13-phi($10))), "spacetimetoy_periodic_fine_rl2/phi.xl" i 12*t u 10:(16*($13-phi($10)))
