#include "cctk.h"

#include "loopcontrol.h"



module loopcontrol
  use loopcontrol_types
  implicit none
  
  interface
     
     subroutine lc_stats_init(stats, name)
       use loopcontrol_types
       implicit none
       CCTK_POINTER :: stats
       character(*) :: name
     end subroutine lc_stats_init
     
     subroutine lc_control_init( &
          control, stats, &
          imin, jmin, kmin, &
          imax, jmax, kmax, &
          iash, jash, kash, &
          di, dj, dk)
       use loopcontrol_types
       implicit none
       type(lc_control_t) :: control
       CCTK_POINTER       :: stats
       CCTK_POINTER       :: imin, jmin, kmin
       CCTK_POINTER       :: imax, jmax, kmax
       CCTK_POINTER       :: iash, jash, kash
       CCTK_POINTER       :: di, dj, dk
     end subroutine lc_control_init

     subroutine lc_control_finish(control, stats)
       use loopcontrol_types
       implicit none
       type(lc_control_t) :: control
       CCTK_POINTER       :: stats
     end subroutine lc_control_finish
     
     subroutine lc_thread_init(control)
       use loopcontrol_types
       implicit none
       type(lc_control_t) :: control
     end subroutine lc_thread_init
     
     logical function lc_thread_done(control)
       use loopcontrol_types
       implicit none
       type(lc_control_t) :: control
     end function lc_thread_done
     
     subroutine lc_thread_step(control)
       use loopcontrol_types
       implicit none
       type(lc_control_t) :: control
     end subroutine lc_thread_step
     
     subroutine lc_get_fortran_type_sizes(type_sizes)
       use loopcontrol_types
       implicit none
       CCTK_POINTER :: type_sizes(4)
     end subroutine lc_get_fortran_type_sizes
     
  end interface
  
end module loopcontrol
