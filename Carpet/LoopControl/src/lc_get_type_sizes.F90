#include "cctk.h"
#include "cctk_Functions.h"

subroutine lc_get_fortran_type_sizes (lc_control_size, lc_statmap_size)
  
  use loopcontrol_types
  
  implicit none
  
  DECLARE_CCTK_FUNCTIONS
  
  integer, intent(out) :: lc_control_size, lc_statmap_size
  
  type(lc_control_t), dimension(2) :: lc_lc
  type(lc_statmap_t), dimension(2) :: lc_lm
  
  lc_control_size = CCTK_PointerTo(lc_lc(2)) - CCTK_PointerTo(lc_lc(1))
  lc_statmap_size = CCTK_PointerTo(lc_lm(2)) - CCTK_PointerTo(lc_lm(1))
  
end subroutine lc_get_fortran_type_sizes
