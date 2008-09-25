#ifndef LC_LOOPCONTROL_H
#define LC_LOOPCONTROL_H

/* This file uses the namespace LC_* for macros and lc_* for C
   identifiers.  */

#include <cctk.h>

#ifdef CCODE

#include <cctk_Arguments.h>

#ifdef __cplusplus
extern "C" {
#endif



#ifdef __cplusplus
#  ifdef CCTK_CXX_RESTRICT
#    define restrict CCTK_CXX_RESTRICT
#  endif
#endif



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

#define LC_ENDLOOP3(name)                       \
  } while (0)

#endif



#if 0
/* Vector size */
#define LC_VECTORSIZE 2         /* Correct for double precision on
                                   Intel */
#endif



/* A topology */
typedef struct lc_topology_t {
  int nthreads[3];
} lc_topology_t;

/* A tiling specification */
typedef struct lc_tiling_t {
  int npoints;
} lc_tiling_t;



typedef struct lc_state_t {
  int topology;
  int tiling[3];
} lc_state_t;



/* For simulated annealing */
typedef struct lc_auto_state_t lc_auto_state_t;

/* For hill climbing */
typedef struct lc_hill_state_t lc_hill_state_t;



/* Statistics for one control parameter set (thread topology and
   tiling specification) of one user parameter set of one loop */
typedef struct lc_stattime_t {
  struct lc_stattime_t * next;
  
  /* Keys */
  
  lc_state_t state;
  int inthreads, jnthreads, knthreads;
  int inpoints, jnpoints, knpoints;
  
  /* Data */
  
  /* Statistics */
  double time_count;            /* number of calls and threads */
  double time_setup_sum, time_setup_sum2; /* time spent setting up loops */
  double time_calc_sum, time_calc_sum2;   /* time spent iterating */
  
  double last_updated;          /* wall time tag */
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
  int * restrict topology_ntilings[3]; /* [dim][topology] */
  
  /* Simulated annealing state */
  lc_auto_state_t * auto_state;
  
  /* Hill climbing state */
  lc_hill_state_t * hill_state;
  
  lc_stattime_t * stattime_list;
  
  /* Statistics */
  double time_count;            /* number of calls and threads */
  double time_setup_sum, time_setup_sum2; /* time spent setting up loops */
  double time_calc_sum, time_calc_sum2;   /* time spent iterating */
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



static inline
int
lc_state_valid (lc_statset_t const * restrict const ls,
                lc_state_t const * restrict const state)
{
  if (state->topology >= 0 && state->topology < ls->ntopologies) {
    int const * restrict const ntilings =
      ls->topology_ntilings[state->topology];
    return (state->tiling[0] >= 0 && state->tiling[0] < ntilings[0] &&
            state->tiling[1] >= 0 && state->tiling[1] < ntilings[1] &&
            state->tiling[2] >= 0 && state->tiling[2] < ntilings[2]);
  }
  return 0;
}

static inline
int
lc_state_equal (lc_state_t const * restrict const state1,
                lc_state_t const * restrict const state2)
{
  return (state1->topology == state2->topology &&
          state1->tiling[0] == state2->tiling[0] &&
          state1->tiling[1] == state2->tiling[1] &&
          state1->tiling[2] == state2->tiling[2]);
}



void
lc_stattime_init (lc_stattime_t * restrict const lt,
                  lc_statset_t * restrict const ls,
                  lc_state_t const * restrict const state);

lc_stattime_t *
lc_stattime_find (lc_statset_t const * restrict const ls,
                  lc_state_t const * restrict const state);

lc_stattime_t *
lc_stattime_find_create (lc_statset_t * restrict const ls,
                         lc_state_t const * restrict const state);



/* TODO: introduce type for num_threads and npoints[3] */
void
lc_statset_init (lc_statset_t * restrict const ls,
                 lc_statmap_t * restrict const lm,
                 int const num_threads,
                 int const npoints[3]);

lc_statset_t *
lc_statset_find (lc_statmap_t const * restrict const lm,
                 int const num_threads,
                 int const npoints[3]);

lc_statset_t *
lc_statset_find_create (lc_statmap_t * restrict const lm,
                        int const num_threads,
                        int const npoints[3]);



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
lc_statmap_init (int * restrict initialised,
                 lc_statmap_t * restrict ls,
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
    static int lc_initialised = 0;                                      \
    static lc_statmap_t lc_lm;                                          \
    if (! lc_initialised) {                                             \
      lc_statmap_init (& lc_initialised, & lc_lm, #name);               \
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
          for (int k = lc_kk; k < lc_kmax; ++k) {                       \
            for (int j = lc_jj; j < lc_jmax; ++j) {                     \
              int const lc_imin = lc_ii;                                \
              LC_PRELOOP_STATEMENTS                                     \
              {                                                         \
                for (int i = lc_imin; i < lc_imax; ++i) {

#define LC_ENDLOOP3(name)                       \
  }                                             \
                }                               \
              LC_POSTLOOP_STATEMENTS            \
              }                                 \
            }                                   \
          }                                     \
        }                                       \
      }                                         \
    lc_control_finish (& lc_lc);                \
    } while (0)

/* Pre- and post loop statements are inserted around the innermost
   loop, which is executed serially.  By default these are empty.  */
#define LC_PRELOOP_STATEMENTS
#define LC_POSTLOOP_STATEMENTS



void
lc_printstats (CCTK_ARGUMENTS);

#ifdef __cplusplus
}
#endif

#endif



#ifdef FCODE
#  include "loopcontrol_fortran.h"
#endif

#endif  /* ifndef LC_LOOPCONTROL_H */
