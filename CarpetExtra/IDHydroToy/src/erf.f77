      FUNCTION erf(x)
      implicit none
      REAL*8 erf,x
CU    USES gammp
      REAL*8 gammp
      if(x.lt.0)then
        erf=-gammp(.5d0,x**2)
      else
        erf=gammp(.5d0,x**2)
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software t4.
