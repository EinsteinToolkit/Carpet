#include <assert.h>
#include <errno.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sys/time.h>

#ifdef _OPENMP
#  include <omp.h>
#endif

#include <cctk.h>
#include <cctk_Parameters.h>

/* #ifdef HAVE_TGMATH_H */
/* #  include <tgmath.h> */
/* #endif */

#include "loopcontrol.h"

#include "lc_auto.h"
#include "lc_hill.h"



#ifndef _OPENMP
/* Replacements for some OpenMP routines if OpenMP is not available */

static inline
int
omp_get_thread_num (void)
{
  return 0;
}

static inline
int
omp_get_num_threads (void)
{
  return 1;
}

static inline
double
omp_get_wtime (void)
{
  struct timeval tv;
  gettimeofday (& tv, NULL);
  return tv.tv_sec + 1.0e-6 * tv.tv_usec;
}

#endif



/* Linked list of all loop statistics structures */
lc_statmap_t * lc_statmap_list = NULL;



/* Find all possible thread topologies */
/* This finds all possible thread topologies which can be expressed as
   NIxNJxNK.  More complex topologies, e.g. based on a recursive
   subdiviston, are not considered (and cannot be expressed in the
   data structures used in LoopControl).  I think more complex
   topologies are not necessary, since the number of treads is usually
   quite small and contains many small factors in its prime
   decomposition.  */
static
void
find_thread_topologies (lc_topology_t * restrict const topologies,
                        const int maxntopologies,
                        int * restrict const ntopologies,
                        int const nthreads)
{
  * ntopologies = 0;
  for (int nk=1; nk<=nthreads; ++nk) {
    if (nthreads % nk == 0) {
      for (int nj=1; nj<=nthreads/nk; ++nj) {
        if (nthreads % (nj*nk) == 0) {
          int const ni = nthreads/(nj*nk);
          if (nthreads == ni*nj*nk) {
            assert (* ntopologies < maxntopologies);
            topologies[* ntopologies].nthreads[0] = ni;
            topologies[* ntopologies].nthreads[1] = nj;
            topologies[* ntopologies].nthreads[2] = nk;
            ++ * ntopologies;
          }
        }
      }
    }
  }
}



/* Find "good" tiling specifications */
/* This calculates a subset of all possible thread specifications.
   One aim is to reduce the search space by disregarding some
   specifications.  The other aim is to distribute the specifications
   "equally", so that one does not have to spend much effort
   investigating tiling specifications with very similar properties.
   For example, if there are 200 grid points, then half of the
   possible tiling specifications consists of splitting the domain
   into two subdomains with [100+N, 100-N] points.  This is avoided by
   covering all possible tiling specifications in exponentially
   growing step sizes.  */
#if 0
static
int tiling_compare (const void * const a, const void * const b)
{
  lc_tiling_t const * const aa = a;
  lc_tiling_t const * const bb = b;
  return aa->npoints - bb->npoints;
}

static
void
find_tiling_specifications (lc_tiling_t * restrict const tilings,
                            const int maxntilings,
                            int * restrict const ntilings,
                            int const npoints)
{
  /* In order to reduce the number of possible tilings, require that
     the step sizes differ by more than 10%.  */
  double const distance_factor = 1.1;
  /* Determine the "good" step sizes in two passes: first small step
     sizes from 1 up to snpoints, then large step sizes from npoints
     down to snpoints+1.  */
  int const snpoints = floor (sqrt (npoints));
  /* For N grid points and a minimum spacing factor F, there are at
     most log(N) / log(F) possible tilings.  There will be fewer,
     since the actual spacings will be rounded up to integers.  */
  
  * ntilings = 0;
  
  /* Small step sizes */
  int minnpoints = 0;
  for (int n=1; n<=snpoints; ++n) {
    if ((double) n > minnpoints * distance_factor) {
      assert (* ntilings < maxntilings);
      tilings[* ntilings].npoints = n;
      minnpoints = n;
      ++ * ntilings;
    }
  }
  
  /* Large step sizes */
  int maxnpoints = 1000000000;
  for (int n=npoints; n>snpoints; --n) {
    if (n * distance_factor < (double) maxnpoints) {
      assert (* ntilings < maxntilings);
      tilings[* ntilings].npoints = n;
      maxnpoints = n;
      ++ * ntilings;
    }
  }
  
  /* Sort */
  qsort (tilings, * ntilings, sizeof * tilings, tiling_compare);
}
#endif

static
void
find_tiling_specifications (lc_tiling_t * restrict const tilings,
                            const int maxntilings,
                            int * restrict const ntilings,
                            int const npoints)
{
  /* In order to reduce the number of possible tilings, require that
     the step sizes differ by more than 10%.  */
  double const distance_factor = 1.1;
  /* For N grid points and a minimum spacing factor F, there are at
     most log(N) / log(F) possible tilings.  There will be fewer,
     since the actual spacings will be rounded up to integers.  */
  
  * ntilings = 0;
  
  int minnpoints = 0;
  for (int n=1; n<npoints; ++n) {
    if ((double) n > minnpoints * distance_factor) {
      assert (* ntilings < maxntilings);
      tilings[* ntilings].npoints = n;
      minnpoints = n;
      ++ * ntilings;
    }
  }
  
  assert (* ntilings < maxntilings);
  /* step size should be at least 1, even if there are only 0
     points */
  tilings[* ntilings].npoints = lc_max (npoints, 1);
  ++ * ntilings;
}



/* Initialise control parameter set statistics */
void
lc_stattime_init (lc_stattime_t * restrict const lt,
                  lc_statset_t * restrict const ls,
                  lc_state_t const * restrict const state)
{
  DECLARE_CCTK_PARAMETERS;
  
  /* Check arguments */
  assert (lt);
  assert (ls);
  assert (state);
  
  /*** Topology ***************************************************************/
  
  lt->state.topology = state->topology;
  
  if (state->topology == -1) {
    
    /* User-specified topology */
    lt->inthreads = -1;
    lt->jnthreads = -1;
    lt->knthreads = -1;
    
  } else {
    
    assert (state->topology >= 0 && state->topology < ls->ntopologies);
    
    lt->inthreads = ls->topologies[lt->state.topology].nthreads[0];
    lt->jnthreads = ls->topologies[lt->state.topology].nthreads[1];
    lt->knthreads = ls->topologies[lt->state.topology].nthreads[2];
    
  }
  
  if (debug) {
    printf ("Thread topology #%d [%d,%d,%d]\n",
            lt->state.topology, lt->inthreads, lt->jnthreads, lt->knthreads);
  }
  
  /* Assert thread topology consistency */
  if (lt->state.topology != -1) {
    assert (lt->inthreads >= 1);
    assert (lt->jnthreads >= 1);
    assert (lt->knthreads >= 1);
    assert (lt->inthreads * lt->jnthreads * lt->knthreads == ls->num_threads);
  }
  
  /*** Tilings ****************************************************************/
  
  for (int d=0; d<3; ++d) {
    lt->state.tiling[d] = state->tiling[d];
    if (state->tiling[d] != -1) {
      assert (state->tiling[d] >= 0 &&
              state->tiling[d] < ls->topology_ntilings[d][lt->state.topology]);
    }
  }
  
  if (state->tiling[0] != -1) {
    lt->inpoints = ls->tilings[0][lt->state.tiling[0]].npoints;
  }
  if (state->tiling[1] != -1) {
    lt->jnpoints = ls->tilings[1][lt->state.tiling[1]].npoints;
  }
  if (state->tiling[2] != -1) {
    lt->knpoints = ls->tilings[2][lt->state.tiling[2]].npoints;
  }
  
  if (debug) {
    printf ("Tiling stride [%d,%d,%d]\n",
            lt->inpoints, lt->jnpoints, lt->knpoints);
  }
  
  /* Assert tiling specification consistency */
  if (state->tiling[0] != -1) {
    assert (lt->inpoints > 0);
  }
  if (state->tiling[1] != -1) {
    assert (lt->jnpoints > 0);
  }
  if (state->tiling[2] != -1) {
    assert (lt->knpoints > 0);
  }
  
  
  
  /* Initialise statistics */
  lt->time_count      = 0.0;
  lt->time_setup_sum  = 0.0;
  lt->time_setup_sum2 = 0.0;
  lt->time_calc_sum   = 0.0;
  lt->time_calc_sum2  = 0.0;
  
  lt->last_updated = 0.0;       /* never updated */
  
  
  
  /* Append to loop statistics list */
  lt->next = ls->stattime_list;
  ls->stattime_list = lt;
}

lc_stattime_t *
lc_stattime_find (lc_statset_t const * restrict const ls,
                  lc_state_t const * restrict const state)
{
  assert (ls);
  
  lc_stattime_t * lt;
  
  for (lt = ls->stattime_list; lt; lt = lt->next) {
    if (lc_state_equal (& lt->state, state)) {
      break;
    }
  }
  
  return lt;
}

lc_stattime_t *
lc_stattime_find_create (lc_statset_t * restrict const ls,
                         lc_state_t const * restrict const state)
{
  assert (ls);
  
  lc_stattime_t * lt;
  
  for (lt = ls->stattime_list; lt; lt = lt->next) {
    if (lc_state_equal (& lt->state, state)) {
      break;
    }
  }
  
  if (! lt) {
    lt = malloc (sizeof * lt);
    lc_stattime_init (lt, ls, state);
  }
  
  assert (lt);
  return lt;
}



/* Initialise user parameter set statistics */
void
lc_statset_init (lc_statset_t * restrict const ls,
                 lc_statmap_t * restrict const lm,
                 int const num_threads,
                 int const npoints[3])
{
  DECLARE_CCTK_PARAMETERS;
  
  /* Check arguments */
  assert (ls);
  assert (lm);
  assert (num_threads >= 1);
  for (int d=0; d<3; ++d) {
    assert (npoints[d] >= 0);
  }
  
  /*** Threads ****************************************************************/
  
  ls->num_threads = num_threads;
  
  /* For up to 1024 threads, there are at most 270 possible
     topologies */
  int const maxntopologies = 1000;
  if (debug) {
    printf ("Running on %d threads\n", ls->num_threads);
  }
  ls->topologies = malloc (maxntopologies * sizeof * ls->topologies);
  find_thread_topologies
    (ls->topologies, maxntopologies, & ls->ntopologies, ls->num_threads);
#if 0
  ls->topologies =
    realloc (ls->topologies, ls->ntopologies * sizeof * ls->topologies);
#endif
  if (debug) {
    printf ("Found %d possible thread topologies\n", ls->ntopologies);
    for (int n = 0; n < ls->ntopologies; ++n) {
      printf ("   %2d: %2d %2d %2d\n",
              n,
              ls->topologies[n].nthreads[0],
              ls->topologies[n].nthreads[1],
              ls->topologies[n].nthreads[2]);
    }
  }
  assert (ls->ntopologies > 0);
  
  /*** Tilings ****************************************************************/
  
  for (int d=0; d<3; ++d) {
    ls->npoints[d] = npoints[d];
  }
  
  /* For up to 1000000 grid points, there are at most 126 possible
     tilings (assuming a minimum distance of 10%) */
  int const maxntilings = 1000;
  for (int d=0; d<3; ++d) {
    if (debug) {
      printf ("Dimension %d: %d points\n", d, ls->npoints[d]);
    }
    ls->tilings[d] = malloc (maxntilings * sizeof * ls->tilings[d]);
    find_tiling_specifications
      (ls->tilings[d], maxntilings, & ls->ntilings[d], ls->npoints[d]);
    ls->topology_ntilings[d] =
      malloc (ls->ntopologies * sizeof * ls->topology_ntilings[d]);
    for (int n = 0; n < ls->ntopologies; ++n) {
      int tiling;
      for (tiling = 1; tiling < ls->ntilings[d]; ++tiling) {
        if (ls->tilings[d][tiling].npoints * ls->topologies[n].nthreads[d] >
            ls->npoints[d])
        {
          break;
        }
      }
      if (tiling == 0) {
        /* Always allow at least one tiling */
        tiling = 1;
      }
      ls->topology_ntilings[d][n] = tiling;
    }
    if (debug) {
      printf ("   Found %d possible tilings\n", ls->ntilings[d]);
      printf ("     ");
      for (int n = 0; n < ls->ntilings[d]; ++n) {
        printf (" %d", ls->tilings[d][n].npoints);
      }
      printf ("\n");
    }
  }
  
  
  
  /* Simulated annealing state */
  ls->auto_state = NULL;
  
  /* Hill climbing state */
  ls->hill_state = NULL;
  
  
  
  /* Initialise list */
  ls->stattime_list = NULL;
  
  /* Initialise statistics */
  ls->time_count      = 0.0;
  ls->time_setup_sum  = 0.0;
  ls->time_setup_sum2 = 0.0;
  ls->time_calc_sum   = 0.0;
  ls->time_calc_sum2  = 0.0;
  
  /* Append to loop statistics list */
  ls->next = lm->statset_list;
  lm->statset_list = ls;
}

lc_statset_t *
lc_statset_find (lc_statmap_t const * restrict const lm,
                 int const num_threads,
                 int const npoints[3])
{
  assert (lm);
  
  lc_statset_t * ls;
  
  for (ls = lm->statset_list; ls; ls = ls->next) {
    if (ls->num_threads == num_threads &&
        ls->npoints[0] == npoints[0] &&
        ls->npoints[1] == npoints[1] &&
        ls->npoints[2] == npoints[2])
    {
      break;
    }
  }
  
  return ls;
}

lc_statset_t *
lc_statset_find_create (lc_statmap_t * restrict const lm,
                        int const num_threads,
                        int const npoints[3])
{
  assert (lm);
  
  lc_statset_t * ls;
  
  for (ls = lm->statset_list; ls; ls = ls->next) {
    if (ls->num_threads == num_threads &&
        ls->npoints[0] == npoints[0] &&
        ls->npoints[1] == npoints[1] &&
        ls->npoints[2] == npoints[2])
    {
      break;
    }
  }
  
  if (! ls) {
    ls = malloc (sizeof * ls);
    lc_statset_init (ls, lm, num_threads, npoints);
  }
  
  assert (ls);
  return ls;
}
                         


/* Initialise loop statistics */
void
lc_statmap_init (int * restrict const initialised,
                 lc_statmap_t * restrict const lm,
                 char const * restrict const name)
{
  /* Check arguments */
  assert (initialised);
  assert (lm);
  
#pragma omp single
  {
    
    /* Set name */
    lm->name = strdup (name);
    
    /* Initialise list */
    lm->statset_list = NULL;
    
    /* Append to loop statistics list */
    lm->next = lc_statmap_list;
    lc_statmap_list = lm;
    
  }
  
#pragma omp single
  {
    /* Set this flag only after initialising */
    * initialised = 1;
  }
}



void
lc_control_init (lc_control_t * restrict const lc,
                 lc_statmap_t * restrict const lm,
                 int const imin, int const jmin, int const kmin,
                 int const imax, int const jmax, int const kmax,
                 int const ilsh, int const jlsh, int const klsh)
{
  DECLARE_CCTK_PARAMETERS;
  
  /* Check arguments */
  assert (lc);
  
  /* Timer */
  lc->time_setup_begin = omp_get_wtime();
  
  /* Check arguments */
  assert (imin >= 0 && imax <= ilsh && ilsh >= 0);
  assert (jmin >= 0 && jmax <= jlsh && jlsh >= 0);
  assert (kmin >= 0 && kmax <= klsh && klsh >= 0);
  
  /* Copy arguments */
  lc->imin = imin;
  lc->jmin = jmin;
  lc->kmin = kmin;
  lc->imax = imax;
  lc->jmax = jmax;
  lc->kmax = kmax;
  lc->ilsh = ilsh;
  lc->jlsh = jlsh;
  lc->klsh = klsh;
  
  
  
  lc_statset_t * restrict ls;
#pragma omp single copyprivate (ls)
  {
    /* Get number of threads */
    int const num_threads = omp_get_num_threads();
    
    /* Calculate number of points */
    int npoints[3];
    npoints[0] = lc_max (imax - imin, 0);
    npoints[1] = lc_max (jmax - jmin, 0);
    npoints[2] = lc_max (kmax - kmin, 0);
    
    ls = lc_statset_find_create (lm, num_threads, npoints);
  }
  
  
  
  lc_stattime_t * restrict lt;
#pragma omp single copyprivate (lt)
  {
    
    lc_state_t state;
    
    /* Select topology */
    
    if (lc_inthreads != -1 || lc_jnthreads != -1 || lc_knthreads != -1)
    {
      /* User-specified thread topology */
      
      if (lc_inthreads == -1 || lc_jnthreads == -1 || lc_knthreads == -1) {
        CCTK_VWarn (CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Illegal thread topology [%d,%d,%d] specified",
                    (int)lc_inthreads, (int)lc_jnthreads, (int)lc_knthreads);
      }
      if (lc_inthreads * lc_jnthreads * lc_knthreads != ls->num_threads) {
        CCTK_VWarn (CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Specified thread topology [%d,%d,%d] is not compatible with the number of threads %d",
                    (int)lc_inthreads, (int)lc_jnthreads, (int)lc_knthreads,
                    ls->num_threads);
      }
      
      state.topology = -1;
      
    } else {
      
      /* Split in the k direction */
      
      for (state.topology = ls->ntopologies - 1;
           state.topology >= 0;
           -- state.topology)
      {
        int have_tilings = 1;
        for (int d=0; d<3; ++d) {
          have_tilings = have_tilings &&
            ls->topology_ntilings[d][state.topology] > 0;
        }
        if (have_tilings) break;
      }
      if (state.topology < 0) {
        assert (0);
        CCTK_WARN (CCTK_WARN_ABORT, "grid too small");
      }
    }
    
    /* Select tiling */
    
    if (lc_inpoints != -1) {
      /* User-specified tiling */
      state.tiling[0] = -1;
    } else {
      /* as many points as possible */
      assert (state.topology >= 0);
      state.tiling[0] = ls->topology_ntilings[0][state.topology] - 1;
    }
    
    if (lc_jnpoints != -1) {
      /* User-specified tiling */
      state.tiling[1] = -1;
    } else {
      if (cycle_j_tilings) {
        /* cycle through all tilings */
        static int count = 0;
        assert (state.topology >= 0);
        state.tiling[1] = (count ++) % ls->topology_ntilings[1][state.topology];
      } else if (legacy_init) {
        /* as many points as possible */
        assert (state.topology >= 0);
        state.tiling[1] = ls->topology_ntilings[1][state.topology] - 1;
      } else {
        /* as few points as possible */
        state.tiling[1] = 0;
      }
    }
    
    if (lc_knpoints != -1) {
      /* User-specified tiling */
      state.tiling[2] = -1;
    } else {
      /* as many points as possible */
      assert (state.topology >= 0);
      state.tiling[2] = ls->topology_ntilings[2][state.topology] - 1;
    }
    
    /* Use simulated annealing to find the best loop configuration */
    if (use_simulated_annealing) {
      lc_auto_init (ls, & state);
    }
    /* Use hill climbing to find the best loop configuration */
    if (use_random_restart_hill_climbing) {
      lc_hill_init (ls, & state);
    }
    
    /* Find or create database entry */
    
    lt = lc_stattime_find_create (ls, & state);
    
    /* Topology */
    
    if (state.topology == -1) {
      /* User-specified topology */
      lt->inthreads = lc_inthreads;
      lt->jnthreads = lc_jnthreads;
      lt->knthreads = lc_knthreads;
    }
    
    /* Assert thread topology consistency */
    assert (lt->inthreads >= 1);
    assert (lt->jnthreads >= 1);
    assert (lt->knthreads >= 1);
    assert (lt->inthreads * lt->jnthreads * lt->knthreads == ls->num_threads);
    
    /* Tilings */
    
    if (state.tiling[0] == -1) {
      /* User-specified tiling */
      lt->inpoints = lc_inpoints;
    }
    if (state.tiling[1] == -1) {
      /* User-specified tiling */
      lt->jnpoints = lc_jnpoints;
    }
    if (state.tiling[2] == -1) {
      /* User-specified tiling */
      lt->knpoints = lc_knpoints;
    }
    
    /* Assert tiling specification consistency */
    assert (lt->inpoints > 0);
    assert (lt->jnpoints > 0);
    assert (lt->knpoints > 0);
    
  } /* omp single */
  
  
  
  lc->statmap  = lm;
  lc->statset  = ls;
  lc->stattime = lt;
  
  

  /*** Threads ****************************************************************/
  
  
  /* Thread loop settings */
  lc->iiimin = imin;
  lc->jjjmin = jmin;
  lc->kkkmin = kmin;
  lc->iiimax = imax;
  lc->jjjmax = jmax;
  lc->kkkmax = kmax;
  lc->iiistep = (lc->iiimax - lc->iiimin + lt->inthreads-1) / lt->inthreads;
  lc->jjjstep = (lc->jjjmax - lc->jjjmin + lt->jnthreads-1) / lt->jnthreads;
  lc->kkkstep = (lc->kkkmax - lc->kkkmin + lt->knthreads-1) / lt->knthreads;
  
#if 0
  /* Correct threading for vectorisation (cache line size) */
  lc->iiistep =
    (lc->iiistep + LC_VECTORSIZE - 1) / LC_VECTORSIZE * LC_VECTORSIZE;
#endif
  
  /* Find location of current thread */
  lc->thread_num  = omp_get_thread_num();
  int c = lc->thread_num;
  int const ci = c % lt->inthreads; c /= lt->inthreads;
  int const cj = c % lt->jnthreads; c /= lt->jnthreads;
  int const ck = c % lt->knthreads; c /= lt->knthreads;
  assert (c == 0);
  lc->iii = lc->iiimin + ci * lc->iiistep;
  lc->jjj = lc->jjjmin + cj * lc->jjjstep;
  lc->kkk = lc->kkkmin + ck * lc->kkkstep;
  
  
  
  /*** Tilings ****************************************************************/
  
  /* Tiling loop settings */
  lc->iimin = lc->iii;
  lc->jjmin = lc->jjj;
  lc->kkmin = lc->kkk;
  lc->iimax = lc_min (lc->iii + lc->iiistep, lc->iiimax);
  lc->jjmax = lc_min (lc->jjj + lc->jjjstep, lc->jjjmax);
  lc->kkmax = lc_min (lc->kkk + lc->kkkstep, lc->kkkmax);
  lc->iistep = lt->inpoints;
  lc->jjstep = lt->jnpoints;
  lc->kkstep = lt->knpoints;
  
#if 0
  /* Correct tiling for vectorisation (cache line size) */
  lc->iistep = (lc->iistep + LC_VECTORSIZE - 1) / LC_VECTORSIZE * LC_VECTORSIZE;
#endif
  
  
  
  /****************************************************************************/
  
  /* Timer */
  lc->time_calc_begin = omp_get_wtime();
}



void
lc_control_finish (lc_control_t * restrict const lc)
{
  lc_stattime_t * restrict const lt = lc->stattime;
  lc_statset_t * restrict const ls = lc->statset;
  
  int ignore_iteration;
#pragma omp single copyprivate (ignore_iteration)
  {
    DECLARE_CCTK_PARAMETERS;
    ignore_iteration = ignore_initial_overhead && lt->time_count == 0.0;
  }
  
  /* Timer */
  double const time_calc_end   = omp_get_wtime();
  double const time_calc_begin = lc->time_calc_begin;
  
  double const time_setup_end   = time_calc_begin;
  double const time_setup_begin = lc->time_setup_begin;
  
  double const time_setup_sum  =
    ignore_iteration ? 0.0 : time_setup_end - time_setup_begin;
  double const time_setup_sum2 = pow (time_setup_sum, 2);
  
  double const time_calc_sum  = time_calc_end - time_calc_begin;
  double const time_calc_sum2 = pow (time_calc_sum, 2);

  /* Update statistics */
#pragma omp critical
  {
    lt->time_count += 1.0;
    
    lt->time_setup_sum  += time_setup_sum;
    lt->time_setup_sum2 += time_setup_sum2;
    
    lt->time_calc_sum  += time_calc_sum;
    lt->time_calc_sum2 += time_calc_sum2;
    
    ls->time_count += 1.0;
    
    ls->time_setup_sum  += time_setup_sum;
    ls->time_setup_sum2 += time_setup_sum2;
    
    ls->time_calc_sum  += time_calc_sum;
    ls->time_calc_sum2 += time_calc_sum2;
  }
  
#pragma omp master
  {
    lt->last_updated = time_calc_end;
  }
  
#pragma omp barrier
  
  {
    DECLARE_CCTK_PARAMETERS;
    if (use_simulated_annealing) {
#pragma omp single
      {
        lc_auto_finish (ls, lt);
      }
    }
    if (use_random_restart_hill_climbing) {
#pragma omp single
      {
        lc_hill_finish (ls, lt);
      }
    }
  }
}



static
double
avg (double const c, double const s)
{
  if (c == 0.0) return 0.0;
  return s / c;
}

static
double
stddev (double const c, double const s, double const s2)
{
  if (c == 0.0) return 0.0;
  return sqrt (s2 / c - pow (s / c, 2));
}



void
lc_printstats (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;
  
  int nmaps = 0;
  for (lc_statmap_t * lm = lc_statmap_list; lm; lm = lm->next) {
    printf ("statmap #%d \"%s\":\n",
            nmaps,
            lm->name);
    int nsets = 0;
    for (lc_statset_t * ls = lm->statset_list; ls; ls = ls->next) {
      printf ("   statset #%d nthreads=%d npoints=[%d,%d,%d]\n",
              nsets,
              ls->num_threads, ls->npoints[0], ls->npoints[1], ls->npoints[2]);
      double sum_count = 0.0;
      double sum_setup = 0.0;
      double sum_calc  = 0.0;
      double min_calc = DBL_MAX;
      int imin_calc = -1;
      int ntimes = 0;
      for (lc_stattime_t * lt = ls->stattime_list; lt; lt = lt->next) {
        printf ("      stattime #%d topology=%d [%d,%d,%d] tiling=[%d,%d,%d]\n",
                ntimes,
                lt->state.topology, lt->inthreads, lt->jnthreads, lt->knthreads,
                lt->inpoints, lt->jnpoints, lt->knpoints);
        double const count = lt->time_count;
        double const setup = lt->time_setup_sum / count;
        double const calc  = lt->time_calc_sum / count;
        printf ("         count: %g   setup: %g   calc: %g\n",
                count, setup, calc);
        sum_count += lt->time_count;
        sum_setup += lt->time_setup_sum;
        sum_calc  += lt->time_calc_sum;
        if (calc < min_calc) {
          min_calc = calc;
          imin_calc = ntimes;
        }
        ++ ntimes;
      }
      double const avg_calc = sum_calc / sum_count;
      printf ("      total count: %g   total setup: %g   total calc: %g\n",
              sum_count, sum_setup, sum_calc);
      printf ("      avg calc: %g   min calc: %g (#%d)\n",
              avg_calc, min_calc, imin_calc);
      ++ nsets;
    }
    ++ nmaps;
  }
}



CCTK_FCALL
void
CCTK_FNAME (lc_statmap_init) (int * restrict const initialised,
                              lc_statmap_t * restrict const lm,
                              ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (name);
  lc_statmap_init (initialised, lm, name);
  free (name);
}

CCTK_FCALL
void
CCTK_FNAME (lc_control_init) (lc_control_t * restrict const lc,
                              lc_statmap_t * restrict const lm,
                              int const * restrict const imin,
                              int const * restrict const jmin,
                              int const * restrict const kmin,
                              int const * restrict const imax,
                              int const * restrict const jmax,
                              int const * restrict const kmax,
                              int const * restrict const ilsh,
                              int const * restrict const jlsh,
                              int const * restrict const klsh)
{
  lc_control_init (lc, lm,
                   * imin - 1, * jmin - 1, * kmin - 1,
                   * imax, * jmax, * kmax,
                   * ilsh, * jlsh, * klsh);
}

CCTK_FCALL
void
CCTK_FNAME (lc_control_finish) (lc_control_t * restrict const lc)
{
  lc_control_finish (lc);
}
