      FUNCTION gammp(a,x)
      implicit none
      REAL*8 a,gammp,x
CU    USES gcf,gser
      REAL*8 gammcf,gamser,gln
      if(x.lt.0.or.a.le.0)pause 'bad arguments in gammp'
      if(x.lt.a+1)then
        call gser(gamser,a,x,gln)
        gammp=gamser
      else
        call gcf(gammcf,a,x,gln)
        gammp=1-gammcf
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software t4.
