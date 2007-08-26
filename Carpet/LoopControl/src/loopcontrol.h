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
typedef struct lc_topology_t {
  int ni, nj, nk;
} lc_topology_t;

/* A tiling specification */
typedef struct lc_tiling_t {
  int npoints;
} lc_tiling_t;



/* For simulated annealing */
typedef struct lc_auto_state_t lc_auto_state_t;



/* Statistics for one control parameter set (thread topology and
   tiling specification) of one user parameter set of one loop */
typedef struct lc_stattime_t {
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
} lc_stattime_t;

/* Statistics for one user parameter set (number of threads and number
   of iterations) of one loop */
typedef struct lc_statset_t {
  struct lc_statset_t * next;
  
  /* Keys */
  
  int num_threads;
  int npoints[3];
  
  /* Data */
  
  /* Thread topologies */
  lc_topology_t * restrict topologies;
  int ntopologies;
  
  /* Tiling specifications */
  lc_tiling_t * restrict tilings[3];
  int ntilings[3];
  
  /* Simulated annealing state */
  lc_auto_state_t * auto_state;
  
  lc_stattime_t * stattime_list;
} lc_statset_t;

/* Statistics for one loop (one source code location) */
typedef struct lc_statmap_t {
  struct lc_statmap_t * next;   /* for linked list */
  
  /* Name */
  char const * restrict name;
  
  lc_statset_t * statset_list;
} lc_statmap_t;

/* Linked list of all loop statistics structures */
extern lc_statmap_t * lc_statmap_list;



typedef struct lc_control_t {
  lc_statmap_t * restrict statmap;
  lc_statset_t * restrict statset;
  lc_stattime_t * restrict stattime;
  
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
} lc_control_t;



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
lc_statmap_init (lc_statmap_t * restrict ls,
                 char const * restrict name);

void
lc_control_init (lc_control_t * restrict lc,
                 lc_statmap_t * restrict lm,
                 int imin, int jmin, int kmin,
                 int imax, int jmax, int kmax,
                 int ilsh, int jlsh, int klsh);

void
lc_control_finish (lc_control_t * restrict lc);



#define LC_LOOP3(name, i,j,k, imin,jmin,kmin, imax,jmax,kmax, ilsh,jlsh,klsh) \
  do {                                                                  \
    static lc_statmap_t lc_lm;                                          \
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
    lc_control_t lc_lc;                                                 \
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
#  include "loopcontrol_fortran.h"
#endif

#endif  /* ifndef LC_LOOPCONTROL_H */
