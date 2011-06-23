#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine AMR_set_level_mask (CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  
!!$  character*1000 :: msg
!!$  
!!$  write (msg, '("Origin: ",3g25.15)') CCTK_ORIGIN_SPACE(:)
!!$  call CCTK_INFO (msg)
!!$  write (msg, '("Delta:  ",3g25.15)') CCTK_DELTA_SPACE(:)
!!$  call CCTK_INFO (msg)
  
  if (any(level_mask<0 .or. level_mask>100)) then
     call CCTK_WARN(CCTK_WARN_ABORT, "level mask contains strange values")
  end if
  
  if (CCTK_EQUALS(radius_type, "diamond")) then
     ! L_1 norm
     level_mask = max_radius / max(abs(x)+abs(y)+abs(z), min_radius)
  else if (CCTK_EQUALS(radius_type, "sphere")) then
     ! L_2 norm
     level_mask = max_radius / max(r, min_radius)
  else if (CCTK_EQUALS(radius_type, "box")) then
     ! L_inf norm
     level_mask = max_radius / max(max(abs(x),abs(y),abs(z)), min_radius)
  else
     call CCTK_WARN(CCTK_WARN_ABORT, "internal error")
  end if
end subroutine AMR_set_level_mask
