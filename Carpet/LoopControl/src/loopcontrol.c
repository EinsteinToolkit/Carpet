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

#include <cctk_Parameters.h>

#include "loopcontrol.h"



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
struct lc_statmap_t * lc_statmap_list = NULL;



/* Find all possible thread topologies */
static
void
find_thread_topologies (struct lc_topology_t * restrict const topologies,
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
            topologies[* ntopologies].ni = ni;
            topologies[* ntopologies].nj = nj;
            topologies[* ntopologies].nk = nk;
            ++ * ntopologies;
          }
        }
      }
    }
  }
}



/* Find "good" tiling specifications */
static
int tiling_compare (const void * const a, const void * const b)
{
  struct lc_tiling_t const * const aa = a;
  struct lc_tiling_t const * const bb = b;
  return aa->npoints - bb->npoints;
}

static
void
find_tiling_specifications (struct lc_tiling_t * restrict const tilings,
                            const int maxntilings,
                            int * restrict const ntilings,
                            int const npoints)
{
  /* In order to reduce the number of possible tilings, require that
     the step sizes differ by more than 10% */
  double const distance_factor = 1.1;
  /* Determine the "good" step sizes in two passes: first small step
     sizes from 1 up to snpoints, then large step sizes from npoints
     down to snpoints+1 */
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



/* Initialise control parameter set statistics */
static
void
lc_stattime_init (struct lc_stattime_t * restrict const lt,
                  struct lc_statset_t * restrict const ls,
                  int const topology,
                  int const tiling[3])
{
  DECLARE_CCTK_PARAMETERS;
  
  /* Check arguments */
  assert (lt);
  assert (ls);
  
  /*** Topology ***************************************************************/
  
  lt->topology = topology;
  
  if (topology == -1) {
    
    /* User-specified topology */
    lt->inthreads = -1;
    lt->jnthreads = -1;
    lt->knthreads = -1;
    
  } else {
    
    assert (topology >= 0 && topology < ls->ntopologies);
    
    lt->inthreads = ls->topologies[lt->topology].ni;
    lt->jnthreads = ls->topologies[lt->topology].nj;
    lt->knthreads = ls->topologies[lt->topology].nk;
    
  }
  
  if (debug) {
    printf ("Thread topology #%d [%d,%d,%d]\n",
            lt->topology, lt->inthreads, lt->jnthreads, lt->knthreads);
  }
  
  /* Assert thread topology consistency */
  if (lt->topology != -1) {
    assert (lt->inthreads >= 1);
    assert (lt->jnthreads >= 1);
    assert (lt->knthreads >= 1);
    assert (lt->inthreads * lt->jnthreads * lt->knthreads == ls->num_threads);
  }
  
  /*** Tilings ****************************************************************/
  
  for (int d=0; d<3; ++d) {
    lt->tiling[d] = tiling[d];
    if (tiling[d] != -1) {
      assert (tiling[d] >= 0 && tiling[d] < ls->ntilings[d]);
    }
  }
  
  if (tiling[0] != -1) {
    lt->inpoints = ls->tilings[0][lt->tiling[0]].npoints;
  }
  if (tiling[1] != -1) {
    lt->jnpoints = ls->tilings[1][lt->tiling[1]].npoints;
  }
  if (tiling[2] != -1) {
    lt->knpoints = ls->tilings[2][lt->tiling[2]].npoints;
  }
  
  if (debug) {
    printf ("Tiling stride [%d,%d,%d]\n",
            lt->inpoints, lt->jnpoints, lt->knpoints);
  }
  
  /* Assert tiling specification consistency */
  if (tiling[0] != -1) {
    assert (lt->inpoints > 0);
  }
  if (tiling[1] != -1) {
    assert (lt->jnpoints > 0);
  }
  if (tiling[2] != -1) {
    assert (lt->knpoints > 0);
  }
  
  
  
  /* Initialise statistics */
  lt->time_count      = 0.0;
  lt->time_setup_sum  = 0.0;
  lt->time_setup_sum2 = 0.0;
  lt->time_calc_sum   = 0.0;
  lt->time_calc_sum2  = 0.0;
  
  /* Append to loop statistics list */
/*   _Pragma ("omp critical") { */
    lt->next = ls->stattime_list;
    ls->stattime_list = lt;
/*   } */
}

static
struct lc_stattime_t *
lc_stattime_find (struct lc_statset_t * restrict const ls,
                  int const topology,
                  int const tiling[3])
{
  assert (ls);
  
  struct lc_stattime_t * lt;
  
  for (lt = ls->stattime_list; lt; lt = lt->next) {
    if (lt->topology == topology &&
        lt->tiling[0] == tiling[0] &&
        lt->tiling[1] == tiling[1] &&
        lt->tiling[2] == tiling[2])
    {
      break;
    }
  }
  
  if (! lt) {
    lt = malloc (sizeof * lt);
    lc_stattime_init (lt, ls, topology, tiling);
  }
  
  return lt;
}



/* Initialise user parameter set statistics */
static
void
lc_statset_init (struct lc_statset_t * restrict const ls,
                 struct lc_statmap_t * restrict const lm,
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
              ls->topologies[n].ni, ls->topologies[n].nj, ls->topologies[n].nk);
    }
  }
  
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
#if 0
    ls->tilings[d] =
      realloc (ls->tilings[d], ls->ntilings[d] * sizeof * ls->tilings[d]);
#endif
    if (debug) {
      printf ("   Found %d possible tilings\n", ls->ntilings[d]);
      printf ("     ");
      for (int n = 0; n < ls->ntilings[d]; ++n) {
        printf (" %d", ls->tilings[d][n].npoints);
      }
      printf ("\n");
    }
  }
  
  
  
  /* Initialise list */
  ls->stattime_list = NULL;
  
  /* Append to loop statistics list */
/*   _Pragma ("omp critical") { */
    ls->next = lm->statset_list;
    lm->statset_list = ls;
/*   } */
}

static
struct lc_statset_t *
lc_statset_find (struct lc_statmap_t * restrict const lm,
                 int const num_threads,
                 int const npoints[3])
{
  assert (lm);
  
  struct lc_statset_t * ls;
  
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
  
  return ls;
}
                         


/* Initialise loop statistics */
void
lc_statmap_init (struct lc_statmap_t * restrict const lm,
                 char const * restrict const name)
{
  /* Check arguments */
  assert (lm);
  
  /* Set name */
  lm->name = strdup (name);
  
  /* Initialise list */
  lm->statset_list = NULL;
  
  /* Append to loop statistics list */
  lm->next = lc_statmap_list;
  lc_statmap_list = lm;
}



void
lc_control_init (struct lc_control_t * restrict const lc,
                 struct lc_statmap_t * restrict const lm,
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
  
  
  
  struct lc_statset_t * restrict ls;
  _Pragma ("omp single copyprivate (ls)") {
    /* Get number of threads */
    int const num_threads = omp_get_num_threads();
    
    /* Calculate number of points */
    int npoints[3];
    npoints[0] = lc_max (imax - imin, 0);
    npoints[1] = lc_max (jmax - jmin, 0);
    npoints[2] = lc_max (kmax - kmin, 0);
    
    ls = lc_statset_find (lm, num_threads, npoints);
  }
  
  
  
  struct lc_stattime_t * restrict lt;
  _Pragma ("omp single copyprivate (lt)") {
    
    /* Select topology */
    
    int topology;
    
    if (lc_inthreads != -1 || lc_jnthreads != -1 || lc_knthreads != -1)
    {
      /* User-specified thread topology */
      
      if (lc_inthreads == -1 || lc_jnthreads == -1 || lc_knthreads == -1) {
        CCTK_VWarn (CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Illegal thread topology [%d,%d,%d] specified\n",
                    (int)lc_inthreads, (int)lc_jnthreads, (int)lc_knthreads);
      }
      if (lc_inthreads * lc_jnthreads * lc_knthreads != ls->num_threads) {
        CCTK_VWarn (CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Specified thread topology [%d,%d,%d] is not compatible with the number of threads %d\n",
                    (int)lc_inthreads, (int)lc_jnthreads, (int)lc_knthreads,
                    ls->num_threads);
      }
      
      topology = -1;
      
    } else {
      
      /* Split in the k direction */
      topology = ls->ntopologies - 1;
      
    }
    
    /* Select tiling */
    
    int tiling[3];
    
    if (lc_inpoints != -1) {
      /* User-specified tiling */
      tiling[0] = -1;
    } else {
      tiling[0] = ls->ntilings[0] - 1; /* as many points as possible */
    }
    
    if (lc_jnpoints != -1) {
      /* User-specified tiling */
      tiling[1] = -1;
    } else {
      if (cycle_j_tilings) {
        /* cycle through all tilings */
        static int count = 0;
        tiling[1] = (count ++) % ls->ntilings[1];
      } else {
        /* as few points as possible */
        tiling[1] = 0;
      }
    }
    
    if (lc_knpoints != -1) {
      /* User-specified tiling */
      tiling[2] = -1;
    } else {
      tiling[2] = ls->ntilings[2] - 1; /* as many points as possible */
    }
    
    /* Find or create database entry */
    
    lt = lc_stattime_find (ls, topology, tiling);
    
    /* Topology */
    
    if (topology == -1) {
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
    
    if (tiling[0] == -1) {
      /* User-specified tiling */
      lt->inpoints = lc_inpoints;
    }
    if (tiling[1] == -1) {
      /* User-specified tiling */
      lt->jnpoints = lc_jnpoints;
    }
    if (tiling[2] == -1) {
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
  
  
  
  /****************************************************************************/
  
  /* Timer */
  lc->time_calc_begin = omp_get_wtime();
}



void
lc_control_finish (struct lc_control_t * restrict const lc)
{
  /* Timer */
  double const time_calc_end   = omp_get_wtime();
  double const time_calc_begin = lc->time_calc_begin;
  
  double const time_setup_end   = time_calc_begin;
  double const time_setup_begin = lc->time_setup_begin;
  
  double const time_setup_sum  = time_setup_end - time_setup_begin;
  double const time_setup_sum2 = pow (time_setup_sum, 2);
  
  double const time_calc_sum  = time_calc_end - time_calc_begin;
  double const time_calc_sum2 = pow (time_calc_sum, 2);
  
  /* Update statistics */
  struct lc_stattime_t * restrict const lt = lc->stattime;
  _Pragma ("omp critical") {
    lt->time_count += 1.0;
    
    lt->time_setup_sum  += time_setup_sum;
    lt->time_setup_sum2 += time_setup_sum2;
    
    lt->time_calc_sum  += time_calc_sum;
    lt->time_calc_sum2 += time_calc_sum2;
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
  if (! verbose) return;
  
  int nmaps = 0;
  for (struct lc_statmap_t * lm = lc_statmap_list; lm; lm = lm->next) {
    printf ("statmap #%d \"%s\":\n",
            nmaps,
            lm->name);
    int nsets = 0;
    for (struct lc_statset_t * ls = lm->statset_list; ls; ls = ls->next) {
      printf ("   statset #%d nthreads=%d npoints=[%d,%d,%d]\n",
              nsets,
              ls->num_threads, ls->npoints[0], ls->npoints[1], ls->npoints[2]);
      double sum_setup = 0.0;
      double sum_calc  = 0.0;
      double min_calc = DBL_MAX;
      int imin_calc = -1;
      int ntimes = 0;
      for (struct lc_stattime_t * lt = ls->stattime_list; lt; lt = lt->next) {
        printf ("      stattime #%d topology=%d [%d,%d,%d] tiling=[%d,%d,%d]\n",
                ntimes,
                lt->topology, lt->inthreads, lt->jnthreads, lt->knthreads,
                lt->inpoints, lt->jnpoints, lt->knpoints);
        printf ("         count: %g   setup: %g   calc: %g\n",
                lt->time_count, lt->time_setup_sum, lt->time_calc_sum);
        sum_setup += lt->time_setup_sum;
        sum_calc  += lt->time_calc_sum;
        if (lt->time_calc_sum < min_calc) {
          min_calc = lt->time_calc_sum;
          imin_calc = ntimes;
        }
        ++ ntimes;
      }
      printf ("      total setup: %g   total calc: %g   min calc: %g (#%d)\n",
              sum_setup, sum_calc, min_calc, imin_calc);
      ++ nsets;
    }
    ++ nmaps;
  }
}



CCTK_FCALL
void
CCTK_FNAME (lc_statmap_init) (struct lc_statmap_t * restrict const lm,
                              ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (name);
  lc_statmap_init (lm, name);
  free (name);
}

CCTK_FCALL
void
CCTK_FNAME (lc_control_init) (struct lc_control_t * restrict const lc,
                              struct lc_statmap_t * restrict const lm,
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
CCTK_FNAME (lc_control_finish) (struct lc_control_t * restrict const lc)
{
  lc_control_finish (lc);
}
