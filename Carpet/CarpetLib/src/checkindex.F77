c     -*-Fortran-*-
c     $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/checkindex.F77,v 1.1 2003/11/05 16:18:39 schnetter Exp $
      
#include "cctk.h"



      subroutine checkindex (i,j,k, di,dj,dk, imax,jmax,kmax, where)
      implicit none
      integer i,j,k
      integer di,dj,dk
      integer imax,jmax,kmax
      character*(*) where
      character*1000 msg
      
      if (      i.lt.1 .or. i+di-1.gt.imax
     $     .or. j.lt.1 .or. j+dj-1.gt.jmax
     $     .or. k.lt.1 .or. k+dk-1.gt.kmax) then
         write (msg, '(a, " array index out of bounds: shape is (",i4,",",i4,",",i4,"), index is (",i4,",",i4,",",i4,"), extent is (",i4,",",i4,",",i4,")")')
     $        where, imax,jmax,kmax, i,j,k, di,dj,dk
         call CCTK_WARN (0, msg(1:len_trim(msg)))
      end if
      end
