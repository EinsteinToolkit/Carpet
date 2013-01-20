#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

#include "loopcontrol.h"



subroutine TestLoopControlFortran_all(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  
  CCTK_LOOP3_ALL_DECLARE(all3)
  integer   :: i,j,k
  
  call CCTK_INFO("TestLoopControlFortran_all")
  fsum_all = 0
  CCTK_LOOP3_ALL(all3, i,j,k)
     fsum_all = fsum_all + r(i,j,k)
  CCTK_ENDLOOP3_ALL(all3)
end subroutine TestLoopControlFortran_all



subroutine TestLoopControlFortran_int(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  
  CCTK_LOOP3_INT_DECLARE(int3)
  integer   :: i,j,k
  
  call CCTK_INFO("TestLoopControlFortran_int")
  fsum_int = 0
  CCTK_LOOP3_INT(int3, i,j,k)
     fsum_int = fsum_int + r(i,j,k)
  CCTK_ENDLOOP3_INT(int3)
end subroutine TestLoopControlFortran_int



subroutine TestLoopControlFortran_bnd(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  
  CCTK_LOOP3_BND_DECLARE(bnd3)
  integer   :: i,j,k
  integer   :: ni,nj,nk
  
  call CCTK_INFO("TestLoopControlFortran_bnd")
  fsum_bnd = 0
  CCTK_LOOP3_BND(bnd3, i,j,k, ni,nj,nk)
     fsum_bnd = fsum_bnd + r(i,j,k)
  CCTK_ENDLOOP3_BND(bnd3)
end subroutine TestLoopControlFortran_bnd



subroutine TestLoopControlFortran_intbnd(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  
  CCTK_LOOP3_INTBND_DECLARE(intbnd3)
  integer   :: i,j,k
  integer   :: ni,nj,nk
  
  call CCTK_INFO("TestLoopControlFortran_intbnd")
  fsum_intbnd = 0
  CCTK_LOOP3_INTBND(intbnd3, i,j,k, ni,nj,nk)
     fsum_intbnd = fsum_intbnd + r(i,j,k)
  CCTK_ENDLOOP3_INTBND(intbnd3)
end subroutine TestLoopControlFortran_intbnd



subroutine TestLoopControlFortran(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  
  call TestLoopControlFortran_all(CCTK_PASS_FTOF)
  call TestLoopControlFortran_int(CCTK_PASS_FTOF)
  call TestLoopControlFortran_bnd(CCTK_PASS_FTOF)
  call TestLoopControlFortran_intbnd(CCTK_PASS_FTOF)
end subroutine TestLoopControlFortran
