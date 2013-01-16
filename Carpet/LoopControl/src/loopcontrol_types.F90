#include "cctk.h"

#include "loopcontrol.h"



module loopcontrol_types
  
  implicit none
  
  ! Note: These types must correspond to the corresponding C types
  ! declared in loopcontrol.h
  
  type, bind(C) :: lc_vec_t
     CCTK_POINTER :: v(LC_DIM)
  end type lc_vec_t
  
  type, bind(C) :: lc_ivec_t
     integer :: v(LC_DIM)
  end type lc_ivec_t
  
  type, bind(C) :: lc_space_t
     type(lc_vec_t)  :: min, max, step, pos
     type(lc_ivec_t) :: count, idx
  end type lc_space_t
  
  type, bind(C) :: lc_control_t
     type(lc_vec_t)   :: ash
     type(lc_space_t) :: loop
     type(lc_space_t) :: thread
     CCTK_POINTER     :: thread_idx_ptr
     integer          :: thread_done
     type(lc_space_t) :: coarse
     type(lc_space_t) :: fine
     type(lc_space_t) :: fine_thread
     CCTK_POINTER     :: fine_thread_info_ptr
     CCTK_POINTER     :: selftest_array
  end type lc_control_t
  
end module loopcontrol_types
