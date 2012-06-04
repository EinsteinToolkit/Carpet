      subroutine CCTK_WARN(a,b)

      implicit none

      integer a
      character*(*) b
      
      write(*,*) b
      stop

      end subroutine CCTK_WARN
