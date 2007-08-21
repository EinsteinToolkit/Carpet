#ifndef LC_LOOPCONTROL_H
#define LC_LOOPCONTROL_H

/* This file uses the namespace LC_* for macros and lc_* for C
   identifiers.  */

#include <cctk.h>

#ifdef CCODE

#include <cctk_Arguments.h>



#if 0
/* The most simple implementation */

#define LC_LOOP3(name, i,j,k, imin,jmin,kmin, imax,jmax,kmax, ilsh,jlsh,klsh) \
  do {                                                                  \
    int const lc_imin = (imin);                                         \
    int const lc_jmin = (jmin);                                         \
    int const lc_kmin = (kmin);                                         \
    int const lc_imax = (imax);                                         \
    int const lc_jmax = (jmax);                                         \
    int const lc_kmax = (kmax);                                         \
    for (int k = lc_kmin; k < lc_kmax; ++k)                             \
      for (int j = lc_jmin; j < lc_jmax; ++j)                           \
        for (int i = lc_imin; i < lc_imax; ++i)

#define LC_ENDLOOP3                             \
  } while (0)

#endif



/* A topology */
struct lc_topology_t {
  int ni, nj, nk;
};

/* A tiling specification */
struct lc_tiling_t {
  int npoints;
};



/* Statistics for one control parameter set (thread topology and
   tiling specification) of one user parameter set of one loop */
struct lc_stattime_t {
  struct lc_stattime_t * next;
  
  /* Keys */
  
  int topology;
  int inthreads, jnthreads, knthreads;
  int tiling[3];
  int inpoints, jnpoints, knpoints;
  
  /* Data */
  
  /* Statistics */
  double time_count;            /* number of calls and threads */
  double time_setup_sum, time_setup_sum2; /* time spent setting up loops */
  double time_calc_sum, time_calc_sum2;   /* time spent iterating */
};

/* Statistics for one user parameter set (number of threads and number
   of iterations) of one loop */
struct lc_statset_t {
  struct lc_statset_t * next;
  
  /* Keys */
  
  int num_threads;
  int npoints[3];
  
  /* Data */
  
  /* Thread topologies */
  struct lc_topology_t * restrict topologies;
  int ntopologies;
  
  /* Tiling specifications */
  struct lc_tiling_t * restrict tilings[3];
  int ntilings[3];
  
  struct lc_stattime_t * stattime_list;
};

/* Statistics for one loop (one source code location) */
struct lc_statmap_t {
  struct lc_statmap_t * next;   /* for linked list */
  
  /* Name */
  char const * restrict name;
  
  struct lc_statset_t * statset_list;
};

/* Linked list of all loop statistics structures */
extern struct lc_statmap_t * lc_statmap_list;



struct lc_control_t {
  struct lc_statmap_t * restrict statmap;
  struct lc_statset_t * restrict statset;
  struct lc_stattime_t * restrict stattime;
  
  /* Copy of arguments (useful for debugging) */
  int imin, jmin, kmin;
  int imax, jmax, kmax;
  int ilsh, jlsh, klsh;
  
  /* Control settings for thread parallelism (useful for debugging) */
  int iiimin, jjjmin, kkkmin;
  int iiimax, jjjmax, kkkmax;
  int iiistep, jjjstep, kkkstep;
  
  /* Control settings for current thread (useful for debugging) */
  int thread_num;
  int iii, jjj, kkk;
  
  /* Control settings for tiling loop */
  int iimin, jjmin, kkmin;
  int iimax, jjmax, kkmax;
  int iistep, jjstep, kkstep;
  
  /* Timing statistics */
  double time_setup_begin, time_calc_begin;
};



static inline
int
lc_min (int const i, int const j)
{
  return i < j ? i : j;
}

static inline
int
lc_max (int const i, int const j)
{
  return i > j ? i : j;
}



void
lc_statmap_init (struct lc_statmap_t * restrict ls,
                 char const * restrict name);

void
lc_control_init (struct lc_control_t * restrict lc,
                 struct lc_statmap_t * restrict lm,
                 int imin, int jmin, int kmin,
                 int imax, int jmax, int kmax,
                 int ilsh, int jlsh, int klsh);

void
lc_control_finish (struct lc_control_t * restrict lc);



#define LC_LOOP3(name, i,j,k, imin,jmin,kmin, imax,jmax,kmax, ilsh,jlsh,klsh) \
  do {                                                                  \
    static struct lc_statmap_t lc_lm;                                   \
    static int lc_initialised = 0;                                      \
    if (! lc_initialised) {                                             \
      _Pragma ("omp single") {                                          \
        lc_statmap_init (& lc_lm, #name);                               \
      }                                                                 \
      _Pragma ("omp single") {                                          \
        /* Set this flag only after initialising */                     \
        lc_initialised = 1;                                             \
      }                                                                 \
    }                                                                   \
    struct lc_control_t lc_lc;                                          \
    lc_control_init (& lc_lc, & lc_lm,                                  \
                     imin,jmin,kmin, imax,jmax,kmax, ilsh,jlsh,klsh);   \
                                                                        \
    /* Coarse loop */                                                   \
    for (int lc_kk = lc_lc.kkmin;                                       \
         lc_kk < lc_lc.kkmax;                                           \
         lc_kk += lc_lc.kkstep)                                         \
    {                                                                   \
      int const lc_kmax = lc_min (lc_kk + lc_lc.kkstep, lc_lc.kkmax);   \
                                                                        \
      for (int lc_jj = lc_lc.jjmin;                                     \
           lc_jj < lc_lc.jjmax;                                         \
           lc_jj += lc_lc.jjstep)                                       \
      {                                                                 \
        int const lc_jmax = lc_min (lc_jj + lc_lc.jjstep, lc_lc.jjmax); \
                                                                        \
        for (int lc_ii = lc_lc.iimin;                                   \
             lc_ii < lc_lc.iimax;                                       \
             lc_ii += lc_lc.iistep)                                     \
        {                                                               \
          int const lc_imax = lc_min (lc_ii + lc_lc.iistep, lc_lc.iimax); \
                                                                        \
          /* Fine loop */                                               \
          for (int k = lc_kk; k < lc_kmax; ++k)                         \
            for (int j = lc_jj; j < lc_jmax; ++j)                       \
              for (int i = lc_ii; i < lc_imax; ++i)

#define LC_ENDLOOP3(name)                       \
  }                                             \
  }                                             \
  }                                             \
  lc_control_finish (& lc_lc);                  \
  } while (0)



void
lc_printstats (CCTK_ARGUMENTS);

#endif



#ifdef FCODE

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

#endif

#endif  /* ifndef LC_LOOPCONTROL_H */
