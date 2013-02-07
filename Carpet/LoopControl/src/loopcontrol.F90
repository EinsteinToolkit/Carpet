#include "cctk.h"

#include "loopcontrol.h"



module loopcontrol
  use loopcontrol_types
  implicit none
  
  interface
     
     subroutine lc_stats_init(stats, line, file, name)
       use loopcontrol_types
       implicit none
       CCTK_POINTER :: stats
       integer      :: line
       character(*) :: file
       character(*) :: name
     end subroutine lc_stats_init
     
     subroutine lc_control_init( &
          control, stats, &
          imin, jmin, kmin, &
          imax, jmax, kmax, &
          iash, jash, kash, &
          istr)
       use loopcontrol_types
       implicit none
       type(lc_control_t) :: control
       CCTK_POINTER       :: stats
       integer            :: imin, jmin, kmin
       integer            :: imax, jmax, kmax
       integer            :: iash, jash, kash
       integer            :: istr
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
     
     integer function lc_thread_done(control)
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
