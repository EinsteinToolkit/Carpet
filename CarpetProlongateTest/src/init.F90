#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"



module CarpetProlongateTest
  implicit none
  
contains
  
  ! Set the time scaling factor to a sum of powers, ensuring that the
  ! scaling factor is not zero for t=0. Also, for power_t=0, the
  ! scaling factor is 1 and has no effect.
  function density_time_scale (tt) result (tscale)
    DECLARE_CCTK_PARAMETERS
    CCTK_REAL, intent(in) :: tt
    CCTK_REAL :: tscale
    integer   :: n
    
    tscale = 1
    do n = 1, power_t
       tscale = tscale + tt ** n
    end do
  end function density_time_scale
  
  elemental function density_sum (xx, yy, zz) result (res)
    DECLARE_CCTK_PARAMETERS
    CCTK_INT, parameter :: izero = 0
    CCTK_REAL, intent(in) :: xx, yy, zz
    CCTK_REAL :: res
    
    res = 0 + &
         1 * density(xx,yy,zz, power_x,izero  ,izero  ) + &
         2 * density(xx,yy,zz, izero  ,power_y,izero  ) + &
         3 * density(xx,yy,zz, izero  ,izero  ,power_z) + &
         4 * density(xx,yy,zz, power_x,power_y,izero  ) + &
         5 * density(xx,yy,zz, power_x,izero  ,power_z) + &
         6 * density(xx,yy,zz, izero  ,power_y,power_z) + &
         7 * density(xx,yy,zz, power_x,power_y,power_z)
  end function density_sum
  
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
  
end module CarpetProlongateTest



subroutine CarpetProlongateTest_Init (CCTK_ARGUMENTS)
  use CarpetProlongateTest
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
!!$  print '("CarpetProlongateTest_Init it=",i0," t=",g0.6, " rl=",i0)', &
!!$       cctk_iteration, cctk_time, nint(log(dble(cctk_levfac(1)))/log(2.0d0))
  
  ! Add +1 to coordinates so that the domain is not symmetric about
  ! the origin (which may accidentally cancel out some errors)
  u = density_time_scale(cctk_time) * density_sum(x+1,y+1,z+1)
  
  uscaled = u * product(cctk_delta_space(:))
end subroutine CarpetProlongateTest_Init



subroutine CarpetProlongateTest_Diff (CCTK_ARGUMENTS)
  use CarpetProlongateTest
  implicit none
  DECLARE_CCTK_ARGUMENTS
  
  u0 = density_time_scale(cctk_time) * density_sum(x+1,y+1,z+1)
  du = u - u0
end subroutine CarpetProlongateTest_Diff



subroutine CarpetProlongateTest_InterpInit (CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS_CarpetProlongateTest_InterpInit
  DECLARE_CCTK_PARAMETERS
  
  CCTK_REAL :: dx, dy, dz
  integer   :: lbnd(3), lsh(3), gsh(3)
  integer   :: i,j,k
  integer   :: ierr
  
!!$  print '("CarpetProlongateTest_InterpInit it=",i0," t=",g0.6, " rl=",i0)', &
!!$       cctk_iteration, cctk_time, nint(log(dble(cctk_levfac(1)))/log(2.0d0))
  
  call CCTK_GrouplbndGN (ierr, cctkGH, 3, lbnd, "CarpetProlongateTest::interp_difference")
  if (ierr/=0) call CCTK_WARN (0, "internal error")
  call CCTK_GrouplshGN  (ierr, cctkGH, 3, lsh , "CarpetProlongateTest::interp_difference")
  if (ierr/=0) call CCTK_WARN (0, "internal error")
  call CCTK_GroupgshGN  (ierr, cctkGH, 3, gsh , "CarpetProlongateTest::interp_difference")
  if (ierr/=0) call CCTK_WARN (0, "internal error")
  
  dx = (interp_xmax - interp_xmin) / (gsh(1)-1)
  dy = (interp_ymax - interp_ymin) / (gsh(2)-1)
  dz = (interp_zmax - interp_zmin) / (gsh(3)-1)
  
  do k=1,lsh(3)
     do j=1,lsh(2)
        do i=1,lsh(1)
           interp_x(i,j,k) = interp_xmin + (lbnd(1)+i-1) * dx
           interp_y(i,j,k) = interp_ymin + (lbnd(2)+j-1) * dy
           interp_z(i,j,k) = interp_zmin + (lbnd(3)+k-1) * dz
        end do
     end do
  end do
end subroutine CarpetProlongateTest_InterpInit

subroutine CarpetProlongateTest_Interp (CCTK_ARGUMENTS)
  use cctk
  implicit none
  DECLARE_CCTK_ARGUMENTS_CarpetProlongateTest_Interp
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  
  integer :: lsh(3)
  
  integer   :: len_interpolator
  integer   :: len_interpolator_options
  character :: fort_interpolator*100
  character :: fort_interpolator_options*1000
  integer   :: coord_handle
  integer   :: interp_handle
  integer   :: options_table
  
  integer      :: coord_type
  CCTK_POINTER :: coords(3)
  integer      :: ninputs
  CCTK_INT     :: inputs(1)
  integer      :: input ! CCTK_VarIndex needs an integer target

  integer      :: noutputs
  CCTK_INT     :: output_types(1)
  CCTK_POINTER :: outputs(1)
  integer      :: npoints
  
  integer :: ierr
  
!!$  print '("CarpetProlongateTest_Interp it=",i0," t=",g0.6, " rl=",i0)', &
!!$       cctk_iteration, cctk_time, nint(log(dble(cctk_levfac(1)))/log(2.0d0))
  
  call CCTK_GrouplshGN (ierr, cctkGH, 3, lsh, "CarpetProlongateTest::interp_difference")
  if (ierr/=0) call CCTK_WARN (0, "internal error")
  
  ! Get coordinate system
  call CCTK_CoordSystemHandle (coord_handle, "cart3d")
  if (coord_handle<0) then
     call CCTK_WARN (0, "Could not obtain coordinate system handle")
  end if
  
  ! Get interpolator
  call CCTK_FortranString &
       (len_interpolator, interpolator, fort_interpolator)
  call CCTK_InterpHandle (interp_handle, fort_interpolator)
  if (interp_handle<0) then
     call CCTK_WARN (0, "Could not obtain interpolator handle")
  end if
  
  ! Get interpolator options
  call CCTK_FortranString &
       (len_interpolator_options, interpolator_options, fort_interpolator_options)
  call Util_TableCreateFromString (options_table, fort_interpolator_options)
  if (options_table<0) then
     call CCTK_WARN (0, "Could not create interpolator option table")
  end if
  
  npoints = product(lsh)
  
  coord_type = CCTK_VARIABLE_REAL
  coords = (/ CCTK_PointerTo(interp_x), CCTK_PointerTo(interp_y), CCTK_PointerTo(interp_z) /)
  
  ninputs = 1
  call CCTK_VarIndex(input, "CarpetProlongateTest::u")
  inputs(1) = input
  if (any(inputs<0)) then
     call CCTK_WARN (0, "Could not obtain interpolation variable index")
  end if
  
  noutputs     = 1
  output_types = (/ CCTK_VARIABLE_REAL /)
  outputs      = (/ CCTK_PointerTo(interp_u) /)
  
  call CCTK_InterpGridArrays &
       (ierr, cctkGH, 3, &
       interp_handle, options_table, coord_handle, &
       npoints, coord_type, coords, &
       ninputs, inputs, &
       noutputs, output_types, outputs)
  
  if (ierr/=0) then
     call CCTK_WARN (0, "Interpolator error")
  end if
end subroutine CarpetProlongateTest_Interp

subroutine CarpetProlongateTest_InterpDiff (CCTK_ARGUMENTS)
  use CarpetProlongateTest
  implicit none
  DECLARE_CCTK_ARGUMENTS_CarpetProlongateTest_InterpDiff
  
!!$  print '("CarpetProlongateTest_InterpDiff it=",i0," t=",g0.6, " rl=",i0)', &
!!$       cctk_iteration, cctk_time, nint(log(dble(cctk_levfac(1)))/log(2.0d0))
  
  interp_u0 = &
       density_time_scale(cctk_time) * &
       density_sum(interp_x+1,interp_y+1,interp_z+1)
  interp_du = interp_u - interp_u0
end subroutine CarpetProlongateTest_InterpDiff



subroutine CarpetProlongateTest_NormInit (CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS_CarpetProlongateTest_NormInit
  
  errornorm = 0
  interp_errornorm = 0
end subroutine CarpetProlongateTest_NormInit

subroutine CarpetProlongateTest_NormCalc (CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  
  errornorm = max (errornorm, maxval(abs(du)))
  interp_errornorm = max (interp_errornorm, maxval(abs(interp_du)))
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
  
  call CCTK_ReduceLocScalar &
       (ierr, cctkGH, -1, max_handle, interp_errornorm, tmp, CCTK_VARIABLE_REAL)
  if (ierr/=0) then
     call CCTK_WARN (CCTK_WARN_ABORT, "Could not evaluate norm")
  end if
  interp_errornorm = tmp
end subroutine CarpetProlongateTest_NormReduce
