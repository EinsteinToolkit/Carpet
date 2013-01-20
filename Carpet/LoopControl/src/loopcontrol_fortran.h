/* -*-f90-*- */

#ifndef LOOPCONTROL_FORTRAN_H
#define LOOPCONTROL_FORTRAN_H

#include "cctk.h"




#define LC_COARSE_DECLARE(name, D)                              \
   CCTK_POINTER :: name/**/_cmin/**/D, name/**/_cmax/**/D,      \
                   name/**/_cstep/**/D, name/**/_cpos/**/D
#define LC_COARSE_OMP_PRIVATE(name, D)          \
   name/**/_cmin/**/D, name/**/_cmax/**/D,      \
   name/**/_cstep/**/D, name/**/_cpos/**/D
#define LC_COARSE_SETUP(name, D)                                        \
   name/**/_control%coarse%min%v(D) = name/**/_control%thread%pos%v(D)  && \
   name/**/_control%coarse%max%v(D) =                                   \
     min(name/**/_control%thread%max%v(D),                              \
         name/**/_control%coarse%min%v(D) +                             \
         name/**/_control%thread%step%v(D))                             && \
   name/**/_cmin/**/D = name/**/_control%coarse%min%v(D)                && \
   name/**/_cmax/**/D = name/**/_control%coarse%max%v(D)                && \
   name/**/_cstep/**/D = name/**/_control%coarse%step%v(D)
#define LC_COARSE_LOOP(name, D)                                         \
   do name/**/_cpos/**/D = name/**/_cmin/**/D, name/**/_cmax/**/D,      \
                           name/**/_cstep/**/D

#define LC_FINE_DECLARE(name, I, D)                             \
   CCTK_POINTER :: name/**/_fmin/**/D, name/**/_fmax/**/D,      \
                   name/**/_fstep/**/D, I
#define LC_FINE_OMP_PRIVATE(name, I, D)         \
   name/**/_fmin/**/D, name/**/_fmax/**/D,      \
   name/**/_fstep/**/D, I
#define LC_FINE_SETUP(name, D)                                  \
   name/**/_control%fine%min%v(D) = name/**/_cpos/**/D          && \
   name/**/_control%fine%max%v(D) =                             \
     min(name/**/_control%coarse%max%v(D),                      \
         name/**/_control%fine%min%v(D) +                       \
         name/**/_control%coarse%step%v(D))                     && \
   name/**/_fmin/**/D = name/**/_control%fine%min%v(D)          && \
   name/**/_fmax/**/D = name/**/_control%fine%max%v(D)          && \
   name/**/_fstep/**/D = name/**/_control%fine%step%v(D)
#define LC_FINE_LOOP(name, I, D)                                        \
   do I = name/**/_fmin/**/D, name/**/_fmax/**/D, name/**/_fstep/**/D



#define LC_DECLARE3(name, i,j,k)                                        \
   CCTK_POINTER :: name/**/_ash1, name/**/_ash2, name/**/_ash3          && \
   CCTK_POINTER :: name/**/_align1, name/**/_align2, name/**/_align3    && \
   CCTK_POINTER, save :: name/**/_stats = 0                             && \
   type(lc_control_t) :: name/**/_control                               && \
   LC_COARSE_DECLARE(name, 1)                                           && \
   LC_COARSE_DECLARE(name, 2)                                           && \
   LC_COARSE_DECLARE(name, 3)                                           && \
   LC_FINE_DECLARE(name, i, 1)                                          && \
   LC_FINE_DECLARE(name, j, 2)                                          && \
   LC_FINE_DECLARE(name, k, 3)

#define LC_OMP_PRIVATE(name, i,j,k)                     \
   name/**/_ash1, name/**/_ash2, name/**/_ash3,         \
   name/**/_align1, name/**/_align2, name/**/_align3,   \
   name/**/_control,                                    \
   LC_COARSE_OMP_PRIVATE(name, 1),                      \
   LC_COARSE_OMP_PRIVATE(name, 2),                      \
   LC_COARSE_OMP_PRIVATE(name, 3),                      \
   LC_FINE_OMP_PRIVATE(name, i, 1),                     \
   LC_FINE_OMP_PRIVATE(name, j, 2),                     \
   LC_FINE_OMP_PRIVATE(name, k, 3)



#define LC_LOOP3STR(name, i,j,k, imin_,jmin_,kmin_, imax_,jmax_,kmax_,  \
                    iash_,jash_,kash_, imin,imax, di_)                  \
   name/**/_ash1 = (iash_)                                              && \
   name/**/_ash2 = (jash_)                                              && \
   name/**/_ash3 = (kash_)                                              && \
   name/**/_align1 = (di_)                                              && \
   name/**/_align2 = 1                                                  && \
   name/**/_align3 = 1                                                  && \
                                                                        && \
   call lc_stats_init(name/**/stats, #name)                             && \
   call lc_control_init(name/**/control, name/**/stats,                 \
                        (imin_), (jmin_), (kmin_),                      \
                        (imax_), (jmax_), (kmax_),                      \
                        name/**/ash1, name/**/ash2, name/**/ash3,       \
                        name/**/align1, name/**/align2, name/**/align3) && \
                                                                        && \
   /* Multithreading */                                                 && \
   call lc_thread_init(name/**/control)                                 && \
   do while (.not. lc_thread_done(name/**/control))                     && \
                                                                        && \
      /* Coarse loops */                                                && \
      LC_COARSE_SETUP(3)                                                && \
      LC_COARSE_SETUP(2)                                                && \
      LC_COARSE_SETUP(1)                                                && \
      LC_COARSE_LOOP(3)                                                 && \
      LC_COARSE_LOOP(2)                                                 && \
      LC_COARSE_LOOP(1)                                                 && \
                                                                        && \
         /* Fine loops */                                               && \
         LC_FINE_SETUP(3)                                               && \
         LC_FINE_SETUP(2)                                               && \
         LC_FINE_SETUP(1)                                               && \
         LC_FINE_LOOP(3)                                                && \
         LC_FINE_LOOP(2)                                                && \
         LC_FINE_LOOP(1)

#define LC_ENDLOOP3STR(name)                                    && \
         end do                                                 && \
         end do                                                 && \
         end do                                                 && \
      end do                                                    && \
      end do                                                    && \
      end do                                                    && \
      call lc_thread_step(name/**/control)                      && \
   end do                                                       && \
   call lc_control_finish(name/**/control, name/**/stats)



#define LC_LOOP3(name, i,j,k, imin,jmin,kmin, imax,jmax,kmax, iash,jash,kash) \
  LC_LOOP3STR(name, i,j,k, imin,jmin,kmin, imax,jmax,kmax, iash,jash,kash, 1)
#define LC_ENDLOOP3(name)                       \
  LC_ENDLOOP3STR(name)



/* Replace CCTK_LOOP macros */
#if (!defined CCTK_LOOP3STR_NORMAL_DECLARE ||           \
     !defined CCTK_ENDLOOP3STR_NORMAL_OMP_PRIVATE ||    \
     !defined CCTK_LOOP3STR_NORMAL ||                   \
     !defined CCTK_ENDLOOP3STR_NORMAL)
#  error "internal error"
#endif
#undef CCTK_LOOP3STR_NORMAL_DECLARE
#undef CCTK_LOOP3STR_NORMAL_OMP_PRIVATE
#undef CCTK_LOOP3STR_NORMAL
#undef CCTK_ENDLOOP3STR_NORMAL
#define CCTK_LOOP3STR_NORMAL_DECLARE     LC_LOOP3_DECLARE
#define CCTK_LOOP3STR_NORMAL_OMP_PRIVATE LC_LOOP3_OMP_PRIVATE
#define CCTK_LOOP3STR_NORMAL             LC_LOOP3
#define CCTK_ENDLOOP3STR_NORMAL          LC_ENDLOOP3



#endif  /* #ifndef LOOPCONTROL_FORTRAN_H */
