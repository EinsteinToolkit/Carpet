#include "cctk.h"

#include "loopcontrol.h"



module loopcontrol
  use loopcontrol_types
  implicit none
  
  interface
     
     subroutine LC_descr_init(descr, line, file, name)
       use loopcontrol_types
       implicit none
       CCTK_POINTER :: descr
       integer      :: line
       character(*) :: file
       character(*) :: name
     end subroutine LC_descr_init
     
     subroutine LC_control_init( &
          control, descr, &
          imin, jmin, kmin, &
          imax, jmax, kmax, &
          iash, jash, kash, &
          ialn, ioff, &
          istr)
       use loopcontrol_types
       implicit none
       type(lc_control_t) :: control
       CCTK_POINTER       :: descr
       integer            :: imin, jmin, kmin
       integer            :: imax, jmax, kmax
       integer            :: iash, jash, kash
       integer            :: ialn, ioff
       integer            :: istr
     end subroutine LC_control_init

     subroutine LC_control_finish(control, descr)
       use loopcontrol_types
       implicit none
       type(lc_control_t) :: control
       CCTK_POINTER       :: descr
     end subroutine LC_control_finish
     
     subroutine LC_thread_init(control)
       use loopcontrol_types
       implicit none
       type(lc_control_t) :: control
     end subroutine LC_thread_init
     
     integer function LC_thread_done(control)
       use loopcontrol_types
       implicit none
       type(lc_control_t) :: control
     end function LC_thread_done
     
     subroutine LC_thread_step(control)
       use loopcontrol_types
       implicit none
       type(lc_control_t) :: control
     end subroutine LC_thread_step
     
     subroutine LC_get_fortran_type_sizes(type_sizes)
       use loopcontrol_types
       implicit none
       CCTK_POINTER :: type_sizes(4)
     end subroutine LC_get_fortran_type_sizes
     
  end interface
  
end module loopcontrol
