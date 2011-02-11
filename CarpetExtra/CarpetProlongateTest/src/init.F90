#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"



subroutine CarpetProlongateTest_Init (CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  CCTK_INT,  parameter :: izero = 0
  
  ! Add +1 to coordinates so that the domain is not symmetric about
  ! the origin (which may accidendally cancel out some errors)
  u = 0 + &
       1 * density(x+1,y+1,z+1, power_x,izero  ,izero  ) + &
       2 * density(x+1,y+1,z+1, izero  ,power_y,izero  ) + &
       3 * density(x+1,y+1,z+1, izero  ,izero  ,power_z) + &
       4 * density(x+1,y+1,z+1, power_x,power_y,izero  ) + &
       5 * density(x+1,y+1,z+1, power_x,izero  ,power_z) + &
       6 * density(x+1,y+1,z+1, izero  ,power_y,power_z) + &
       7 * density(x+1,y+1,z+1, power_x,power_y,power_z)
  u0 = u
  
  uscaled = u * product(cctk_delta_space(:))
  
contains
  
  ! Calculate an integrand that is polynomial in the coordinates
  ! (xx,yy,zz), and which is zero if integrated over [0,2]^3
  elemental function density (xx, yy, zz, px, py, pz) result (res)
    CCTK_REAL, intent(in) :: xx, yy, zz
    CCTK_INT,  intent(in) :: px, py, pz
    CCTK_REAL :: res
    
    CCTK_REAL, parameter :: xmin = 0
    CCTK_REAL, parameter :: ymin = 0
    CCTK_REAL, parameter :: zmin = 0
    CCTK_REAL, parameter :: xmax = 2
    CCTK_REAL, parameter :: ymax = 2
    CCTK_REAL, parameter :: zmax = 2
    
    CCTK_REAL :: integral, volume
    
    integral = &
         (xmax**(px+1) - xmin**(px+1)) / (px+1) * &
         (ymax**(py+1) - ymin**(py+1)) / (py+1) * &
         (zmax**(pz+1) - zmin**(pz+1)) / (pz+1)
    volume = (xmax - xmin) * (ymax - ymin) * (zmax - zmin)
    
    res = xx**px * yy**py * zz**pz - integral / volume
    
  end function density
  
end subroutine CarpetProlongateTest_Init



subroutine CarpetProlongateTest_Diff (CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  
  du = u - u0
end subroutine CarpetProlongateTest_Diff



subroutine CarpetProlongateTest_NormInit (CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  
  errornorm = 0
end subroutine CarpetProlongateTest_NormInit

subroutine CarpetProlongateTest_NormCalc (CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  
  errornorm = max (errornorm, maxval(abs(du)))
end subroutine CarpetProlongateTest_NormCalc

subroutine CarpetProlongateTest_NormReduce (CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  
  integer   :: max_handle
  integer   :: ierr
  CCTK_REAL :: tmp
  
  call CCTK_ReductionArrayHandle (max_handle, "maximum")
  if (max_handle < 0) then
     call CCTK_WARN (CCTK_WARN_ABORT, "Could not obtain norm handle")
  end if
  call CCTK_ReduceLocScalar &
       (ierr, cctkGH, -1, max_handle, errornorm, tmp, CCTK_VARIABLE_REAL)
  if (ierr/=0) then
     call CCTK_WARN (CCTK_WARN_ABORT, "Could not evaluate norm")
  end if
  errornorm = tmp
end subroutine CarpetProlongateTest_NormReduce
