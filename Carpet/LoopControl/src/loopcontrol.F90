module loopcontrol
  
  use loopcontrol_types
  
  interface
     
     subroutine lc_statmap_init (lc_lm, name)
       use loopcontrol_types
       implicit none
       type (lc_statmap_t) :: lc_lm
       character(*)        :: name
     end subroutine lc_statmap_init
     
     subroutine lc_control_init (lc_lc, lc_lm, &
          imin,jmin,kmin, imax,jmax,kmax, ilsh,jlsh,klsh)
       use loopcontrol_types
       implicit none
       type (lc_control_t) :: lc_lc
       type (lc_statmap_t) :: lc_lm
       integer, intent(in) :: imin, jmin, kmin
       integer, intent(in) :: imax, jmax, kmax
       integer, intent(in) :: ilsh, jlsh, klsh
     end subroutine lc_control_init
     
     subroutine lc_control_finish (lc_lc)
       use loopcontrol_types
       implicit none
       type (lc_control_t) :: lc_lc
     end subroutine lc_control_finish
     
  end interface
  
end module loopcontrol
