      subroutine CCTK_WARN(a,b)

      implicit none

      integer a
      character*100 b
      
      write(*,*) b
      stop

      end subroutine CCTK_WARN
