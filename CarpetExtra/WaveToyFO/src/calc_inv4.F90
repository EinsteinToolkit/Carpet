! $Header: /home/eschnett/C/carpet/Carpet/CarpetExtra/WaveToyFO/src/Attic/calc_inv4.F90,v 1.1 2004/05/07 20:55:47 schnetter Exp $

#include "cctk.h"

subroutine WaveToyFO_calc_inv4 (a, b)
  use matinv
  implicit none
  CCTK_REAL a(4,4), b(4,4)
  call calc_inv4 (a, b)
end subroutine WaveToyFO_calc_inv4
