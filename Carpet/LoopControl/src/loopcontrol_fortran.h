/* -*-f90-*- */

#ifndef LOOPCONTROL_FORTRAN_H
#define LOOPCONTROL_FORTRAN_H



#define LC_DECLARE3(name, i,j,k)                        &&\
type (lc_statmap_t), save :: name/**/_lm                &&\
logical, save :: name/**/_initialised = .false.         &&\
type (lc_control_t) :: name/**/_lc                      &&\
integer :: name/**/_ii, name/**/_jj, name/**/_kk        &&\
integer :: name/**/_imax, name/**/_jmax, name/**/_kmax  &&\
integer :: i, j, k



#define LC_PRIVATE3(name)                       \
name/**/_lc,                                    \
name/**/_imax, name/**/_jmax, name/**/_kmax



#define LC_LOOP3(name, i,j,k, imin,jmin,kmin, imax,jmax,kmax, ilsh,jlsh,klsh) &&\
if (.not. name/**/_initialised) then                                    &&\
!$omp single                                                            &&\
   call lc_statmap_init (name/**/_lm, "name")                           &&\
!$omp end single                                                        &&\
!$omp single                                                            &&\
   /* Set this flag only after initialising */                          &&\
   name/**/_initialised = .true.                                        &&\
!$omp end single                                                        &&\
end if                                                                  &&\
call lc_control_init (name/**/_lc, name/**/_lm, imin,jmin,kmin, imax,jmax,kmax, ilsh,jlsh,klsh) &&\
                                                                        &&\
/* Coarse loop */                                                       &&\
do name/**/_kk = name/**/_lc%kkmin + 1, name/**/_lc%kkmax, name/**/_lc%kkstep &&\
   name/**/_kmax = min (name/**/_kk - 1 + name/**/_lc%kkstep, name/**/_lc%kkmax) &&\
   do name/**/_jj = name/**/_lc%jjmin + 1, name/**/_lc%jjmax, name/**/_lc%jjstep &&\
      name/**/_jmax = min (name/**/_jj - 1 + name/**/_lc%jjstep, name/**/_lc%jjmax) &&\
      do name/**/_ii = name/**/_lc%iimin + 1, name/**/_lc%iimax, name/**/_lc%iistep &&\
         name/**/_imax = min (name/**/_ii - 1 + name/**/_lc%iistep, name/**/_lc%iimax) &&\
                                                                        &&\
         /* Fine loop */                                                &&\
         do k = name/**/_kk, name/**/_kmax                              &&\
            do j = name/**/_jj, name/**/_jmax                           &&\
               do i = name/**/_ii, name/**/_imax
                  


#define LC_ENDLOOP3(name)                       &&\
               end do                           &&\
            end do                              &&\
         end do                                 &&\
                                                &&\
      end do                                    &&\
   end do                                       &&\
end do                                          &&\
                                                &&\
call lc_control_finish (name/**/_lc)



#endif  /* #ifndef LOOPCONTROL_FORTRAN_H */
