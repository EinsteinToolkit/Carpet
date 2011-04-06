/* -*-f90-*- */

#ifndef LOOPCONTROL_FORTRAN_H
#define LOOPCONTROL_FORTRAN_H



#define LC_DECLARE3(name, i,j,k)                        &&\
type (lc_statmap_t), save :: name/**/_lm                &&\
integer, save :: name/**/_initialised = 0               &&\
type (lc_control_t) :: name/**/_lc                      &&\
integer :: name/**/_ii, name/**/_jj, name/**/_kk        &&\
integer :: name/**/_imin, name/**/_jmin, name/**/_kmin  &&\
integer :: name/**/_imax, name/**/_jmax, name/**/_kmax  &&\
integer :: name/**/_ipos, name/**/_ioffset, name/**/_di &&\
integer :: i, j, k



#define LC_PRIVATE3(name) \
name/**/_lc, \
name/**/_imin, name/**/_jmin, name/**/_kmin, \
name/**/_imax, name/**/_jmax, name/**/_kmax, \
name/**/_ipos, name/**/_ioffset, name/**/_di 



#define LC_LOOP3(name, i,j,k, imin,jmin,kmin, imax,jmax,kmax, ilsh,jlsh,klsh) \
  LC_LOOP3VEC(name, i,j,k, imin,jmin,kmin, imax,jmax,kmax, ilsh,jlsh,klsh, 1)
#define LC_ENDLOOP3(name)                       \
  LC_ENDLOOP3VEC(name)

#define LC_LOOP3VEC(name, i,j,k, imin_,jmin_,kmin_, imax_,jmax_,kmax_, ilsh_,jlsh_,klsh_, di_) &&\
    if (name/**/_initialised .eq. 0) then                                                      &&\
      call lc_statmap_init (name/**/_initialised, name/**/_lm, "name")                         &&\
    end if                                                                                     &&\
    name/**/_di = (di_)                                                                        &&\
    call lc_control_init (name/**/_lc, name/**/_lm,                                              \
                     (imin_), (jmin_), (kmin_),                                                  \
                     (imax_), (jmax_), (kmax_),                                                  \
                     (ilsh_), (jlsh_), (klsh_),                                                  \
                     name/**/_di)                                                              &&\
                                                                                               &&\
    /* Coarse loop */                                                                          &&\
name/**/_lc_loop3vec:                                                                            \
     do name/**/_kk = name/**/_lc%kkmin + 1, name/**/_lc%kkmax, name/**/_lc%kkstep             &&\
      name/**/_kmin = name/**/_kk - 1 + name/**/_lc%kkkkmin                                    &&\
      name/**/_kmax = min (name/**/_kk - 1 + name/**/_lc%kkkkmax, name/**/_lc%kkmax)           &&\
                                                                                               &&\
      do name/**/_jj = name/**/_lc%jjmin + 1, name/**/_lc%jjmax, name/**/_lc%jjstep            &&\
        name/**/_jmin = name/**/_jj - 1 + name/**/_lc%jjjjmin                                  &&\
        name/**/_jmax = min (name/**/_jj - 1 + name/**/_lc%jjjjmax, name/**/_lc%jjmax)         &&\
                                                                                               &&\
        do name/**/_ii = name/**/_lc%iimin + 1, name/**/_lc%iimax, name/**/_lc%iistep          &&\
          name/**/_imin = name/**/_ii - 1 + name/**/_lc%iiiimin                                &&\
          name/**/_imax = min (name/**/_ii - 1 + name/**/_lc%iiiimax, name/**/_lc%iimax)       &&\
                                                                                               &&\
          /* Fine loop */                                                                      &&\
          do k = name/**/_kmin + 1, name/**/_kmax                                              &&\
            do j = name/**/_jmin + 1, name/**/_jmax                                            &&\
              LC_PRELOOP_STATEMENTS                                                            &&\
                name/**/_ipos =                                                                  \
                  name/**/_imin + name/**/_lc%ilsh * (j - 1 + name/**/_lc%jlsh * (k-1))        &&\
                /* round down to the next multiple of lc_di, assuming that lc_di is a power */ &&\
                /* of 2, and integers are encoded as two's-complement */                       &&\
                name/**/_ioffset = iand(name/**/_ipos, - name/**/_di) - name/**/_ipos          &&\
                do i = name/**/_imin + name/**/_ioffset + 1, name/**/_imax, name/**/_di



#define LC_ENDLOOP3VEC(name)                    &&\
               end do                           &&\
            end do                              &&\
            LC_POSTLOOP_STATEMENTS              &&\
         end do                                 &&\
                                                &&\
      end do                                    &&\
   end do                                       &&\
end do name/**/_lc_loop3vec                     &&\
                                                &&\
call lc_control_finish (name/**/_lc)

/* Pre- and post loop statements are inserted around the innermost
   loop, which is executed serially. By default these are empty. */
#define LC_PRELOOP_STATEMENTS
#define LC_POSTLOOP_STATEMENTS

#endif  /* #ifndef LOOPCONTROL_FORTRAN_H */
