# initialise
kx=0.222222222222222
ky=0.151515151515152
kz=0.222222222222222
omega=sqrt(kx**2 + ky**2 + kz**2)



# without refinement

# solutions
p "hydrotoy_periodic/u.zl" u ($8-1):($1==time&&$3==0?$9:0/0) w lp, "hydrotoy_periodic_coarse/u.zl" u (2*($8-1)):(2*$1==time&&$3==0?$9:0/0) w lp, cos((kz*(0.3*x-9) + omega*(0.15*time))*pi) w l

# errors
p "hydrotoy_periodic/u.zl" u ($8-1):($1==time&&$3==0?$9-cos((kz*(0.3*($8-1)-9) + omega*(0.15*time))*pi):0/0) w lp, "hydrotoy_periodic_coarse/u.zl" u (2*($8-1)):(2*$1==time&&$3==0?$9-cos((kz*(0.3*(2*($8-1))-9) + omega*(0.15*time))*pi):0/0) w lp



# with refinement

# solutions
p "hydrotoy_periodic_rl2/u.zl" u ($8-2):($1==time&&$3==0?$9:0/0) w lp, "hydrotoy_periodic_coarse_rl2/u.zl" u (2*($8-2)):(2*$1==time&&$3==0?$9:0/0) w lp, cos((kz*(0.15*x-9) + omega*(0.075*time))*pi) w l

# errors
p "hydrotoy_periodic_rl2/u.zl" u ($8-2):($1==time&&$3==0?$9-cos((kz*(0.15*($8-2)-9) + omega*(0.075*time))*pi):0/0) w lp, "hydrotoy_periodic_coarse_rl2/u.zl" u (2*($8-2)):(2*$1==time&&$3==0?$9-cos((kz*(0.15*(2*($8-2))-9) + omega*(0.075*time))*pi):0/0) w lp
