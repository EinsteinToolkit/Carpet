#ifndef LC_LOOPCONTROL_H
#define LC_LOOPCONTROL_H

/* This file uses the namespace LC_* for macros and lc_* for C
   identifiers. */

#include <cctk.h>



#ifdef CCODE

#ifdef __cplusplus
extern "C" {
#endif

#ifdef __cplusplus
#  ifdef CCTK_CXX_RESTRICT
#    define restrict CCTK_CXX_RESTRICT
#  endif
#endif



/* A topology */
typedef struct lc_topology_t {
  int nthreads[2][3];           /* [0:outer|1:inner][ijk] */
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
  int inithreads, jnithreads, knithreads;
  int inpoints, jnpoints, knpoints;
  
  /* Data */
  
  /* Statistics */
  /* number of calls and threads */
  double time_count, time_count_init;
  /* time spent setting up loops */
  double time_setup_sum, time_setup_sum2;
  /* time spent iterating */
  double time_calc_sum, time_calc_sum2; 
  double time_calc_init;        /* time for first calculation */
  
  /* wall time tag */
  double last_updated;
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
  /* number of calls and threads */
  double time_count, time_count_init;
  /* time spent setting up loops */
  double time_setup_sum, time_setup_sum2;
  /* time spent iterating */
  double time_calc_sum, time_calc_sum2; 
  double time_calc_init;        /* time for first calculation */
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
  CCTK_ATTRIBUTE_PURE;
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
  CCTK_ATTRIBUTE_PURE;
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
                  lc_state_t const * restrict const state)
  CCTK_ATTRIBUTE_PURE;

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
                 int const npoints[3])
  CCTK_ATTRIBUTE_PURE;

lc_statset_t *
lc_statset_find_create (lc_statmap_t * restrict const lm,
                        int const num_threads,
                        int const npoints[3]);



typedef struct lc_control_t {
  lc_statmap_t * restrict statmap;
  lc_statset_t * restrict statset;
  lc_stattime_t * restrict stattime;
  
  /* Copy of arguments (useful for debugging) */
  /* Full domain */
  int imin, jmin, kmin;
  int imax, jmax, kmax;
  int ilsh, jlsh, klsh;
  int di;
  
  /* Control settings for thread parallelism (useful for debugging) */
  /* Outer thread decomposition of full domain */
  int iiimin, jjjmin, kkkmin;
  int iiimax, jjjmax, kkkmax;
  int iiistep, jjjstep, kkkstep;
  
  /* Control settings for current thread (useful for debugging) */
  int thread_num;
  /* Location of this thread in full domain */
  int iii, jjj, kkk;
  /* Index (not location!) of this thread in loop tile */
  int iiii, jjjj, kkkk;
  
  /* Control settings for tiling loop */
  /* Loop tiling decomposition in this thread's domain */
  int iimin, jjmin, kkmin;
  int iimax, jjmax, kkmax;
  int iistep, jjstep, kkstep;
  
  /* Control settings for inner thread parallelism */
  /* Inner thread decomposition, as offsets (!) to loop tiling */
  int iiiimin, jjjjmin, kkkkmin;
  int iiiimax, jjjjmax, kkkkmax;
  int iiiistep, jjjjstep, kkkkstep;
  
  /* Timing statistics */
  double time_setup_begin, time_calc_begin;
  
  /* Self check */
  char * restrict selftest_count;
} lc_control_t;



static inline
int
lc_min (int const i, int const j)
  CCTK_ATTRIBUTE_CONST;
static inline
int
lc_min (int const i, int const j)
{
  return i < j ? i : j;
}

static inline
int
lc_max (int const i, int const j)
  CCTK_ATTRIBUTE_CONST;
static inline
int
lc_max (int const i, int const j)
{
  return i > j ? i : j;
}

/* Align by shifting to the right if necessary */
static inline
int
lc_align (int const i, int const di)
  CCTK_ATTRIBUTE_CONST;
static inline
int
lc_align (int const i, int const di)
{
  return (i + di - 1) / di * di;
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
                 int ilsh, int jlsh, int klsh,
                 int di);

void
lc_control_selftest (lc_control_t * restrict lc,
                     int imin, int imax, int j, int k);

void
lc_control_finish (lc_control_t * restrict lc);



#define LC_LOOP3(name, i,j,k, imin,jmin,kmin, imax,jmax,kmax, ilsh,jlsh,klsh) \
  LC_LOOP3VEC(name, i,j,k, imin,jmin,kmin, imax,jmax,kmax, ilsh,jlsh,klsh, 1)
#define LC_ENDLOOP3(name)                       \
  LC_ENDLOOP3VEC(name)

#define LC_LOOP3VEC(name, i,j,k, imin_,jmin_,kmin_, imax_,jmax_,kmax_, ilsh_,jlsh_,klsh_, di_) \
  do {                                                                  \
    typedef int lc_loop3vec_##name;                                     \
    static int lc_initialised = 0;                                      \
    static lc_statmap_t lc_lm;                                          \
    if (! lc_initialised) {                                             \
      lc_statmap_init (& lc_initialised, & lc_lm, #name);               \
    }                                                                   \
    int const lc_di = (di_);                                            \
    lc_control_t lc_lc;                                                 \
    lc_control_init (& lc_lc, & lc_lm,                                  \
                     (imin_), (jmin_), (kmin_),                         \
                     (imax_), (jmax_), (kmax_),                         \
                     (ilsh_), (jlsh_), (klsh_),                         \
                     lc_di);                                            \
    int const lc_do_selftest = lc_lc.selftest_count != 0;               \
                                                                        \
    /* Coarse loop */                                                   \
    for (int lc_kk = lc_lc.kkmin;                                       \
         lc_kk < lc_lc.kkmax;                                           \
         lc_kk += lc_lc.kkstep)                                         \
    {                                                                   \
      int const lc_kmin = lc_kk + lc_lc.kkkkmin;                        \
      int const lc_kmax =                                               \
        lc_min (lc_kk + lc_min (lc_lc.kkkkmax, lc_lc.kkstep),           \
                lc_lc.kkmax);                                           \
                                                                        \
      for (int lc_jj = lc_lc.jjmin;                                     \
           lc_jj < lc_lc.jjmax;                                         \
           lc_jj += lc_lc.jjstep)                                       \
      {                                                                 \
        int const lc_jmin = lc_jj + lc_lc.jjjjmin;                      \
        int const lc_jmax =                                             \
          lc_min (lc_jj + lc_min (lc_lc.jjjjmax, lc_lc.jjstep),         \
                  lc_lc.jjmax);                                         \
                                                                        \
        for (int lc_ii = lc_lc.iimin;                                   \
             lc_ii < lc_lc.iimax;                                       \
             lc_ii += lc_lc.iistep)                                     \
        {                                                               \
          int const lc_imin = lc_ii + lc_lc.iiiimin;                    \
          int const lc_imax =                                           \
            lc_min (lc_ii + lc_min (lc_lc.iiiimax, lc_lc.iistep),       \
                    lc_lc.iimax);                                       \
                                                                        \
          /* Fine loop */                                               \
          for (int k = lc_kmin; k < lc_kmax; ++k) {                     \
            for (int j = lc_jmin; j < lc_jmax; ++j) {                   \
              LC_PRELOOP_STATEMENTS                                     \
              {                                                         \
                if (CCTK_BUILTIN_EXPECT(lc_do_selftest, 0)) {           \
                  lc_control_selftest (& lc_lc, lc_imin, lc_imax, j, k); \
                }                                                       \
                int const lc_ipos =                                     \
                  lc_imin + lc_lc.ilsh * (j + lc_lc.jlsh * k);          \
                int const lc_ioffset = (lc_ipos & - lc_di) - lc_ipos;   \
                for (int i = lc_imin + lc_ioffset; i < lc_imax; i += lc_di) { 

#define LC_ENDLOOP3VEC(name)                                    \
                }                                               \
              }                                                 \
              LC_POSTLOOP_STATEMENTS                            \
            }                                                   \
          }                                                     \
        }                                                       \
      }                                                         \
    }                                                           \
    lc_control_finish (& lc_lc);                                \
    typedef lc_loop3vec_##name lc_ensure_proper_nesting;        \
  } while (0)

/* Pre- and post loop statements are inserted around the innermost
   loop, which is executed serially. By default these are empty. */
#define LC_PRELOOP_STATEMENTS
#define LC_POSTLOOP_STATEMENTS



/* Replace CCTK_LOOP macros */
#undef CCTK_LOOP3
#undef CCTK_ENDLOOP3
#define CCTK_LOOP3    LC_LOOP3
#define CCTK_ENDLOOP3 LC_ENDLOOP3



#ifdef __cplusplus
}
#endif

#endif



#ifdef FCODE
#  include "loopcontrol_fortran.h"
#endif

#endif  /* ifndef LC_LOOPCONTROL_H */
