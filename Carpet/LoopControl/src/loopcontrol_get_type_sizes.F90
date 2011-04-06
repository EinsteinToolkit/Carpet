     subroutine lc_get_fortran_type_sizes (sum_of_type_sizes)
       
       use loopcontrol_types
       
       implicit none
       
       integer, intent(out) :: sum_of_type_sizes

       integer :: lc_control_size, lc_statmap_size

       type(lc_control_t) :: lc_lc
       type(lc_statmap_t) :: lc_lm

       INQUIRE(IOLENGTH=lc_control_size) lc_lc
       INQUIRE(IOLENGTH=lc_statmap_size) lc_lm

       sum_of_type_sizes = lc_control_size + lc_statmap_size

     end subroutine lc_get_fortran_type_sizes

