module loopcontrol
  
  use loopcontrol_types
  
  interface
     
     subroutine lc_statmap_init (initialised, lc_lm, name)
       use loopcontrol_types
       implicit none
       integer, intent(out) :: initialised
       type (lc_statmap_t) :: lc_lm
       character(*)        :: name
     end subroutine lc_statmap_init
     
     subroutine lc_control_init (lc_lc, lc_lm, &
          imin,jmin,kmin, imax,jmax,kmax, ilsh,jlsh,klsh,di)
       use loopcontrol_types
       implicit none
       type (lc_control_t) :: lc_lc
       type (lc_statmap_t) :: lc_lm
       integer, intent(in) :: imin, jmin, kmin
       integer, intent(in) :: imax, jmax, kmax
       integer, intent(in) :: ilsh, jlsh, klsh
       integer, intent(in) :: di
     end subroutine lc_control_init
     
     subroutine lc_control_finish (lc_lc)
       use loopcontrol_types
       implicit none
       type (lc_control_t) :: lc_lc
     end subroutine lc_control_finish

     subroutine lc_get_fortran_type_sizes (sum_of_type_sizes)
       use loopcontrol_types
       implicit none
       integer, intent(out) :: sum_of_type_sizes
     end subroutine lc_get_fortran_type_sizes

  end interface
  
end module loopcontrol
