#include <assert.h>
#include <math.h>

#include <gsl/gsl_rng.h>

#include <cctk.h>
#include <cctk_Parameters.h>

/* #ifdef HAVE_TGMATH_H */
/* #  include <tgmath.h> */
/* #endif */

#include "lc_siman.h"

#include "lc_auto.h"



static lc_statset_t const * restrict statset;



static inline
int
step (int const oldval, int const maxval,
      gsl_rng const * restrict const rng, double const step_size)
{
  int const offset = (int) ((2 * gsl_rng_uniform_pos (rng) - 1) * step_size);
  return (oldval + offset + maxval) % maxval;
}



static
void
take_step (const gsl_rng *rng, void *xp_, double step_size)
{
  DECLARE_CCTK_PARAMETERS;
  lc_state_t * restrict const xp = xp_;
  if (gsl_rng_uniform (rng) < siman_probability_change_topology) {
    xp->topology = gsl_rng_uniform_int (rng, statset->ntopologies);
  }
  for (int d=0; d<3; ++d) {
    xp->tiling[d] = step (xp->tiling[d], statset->ntilings[d], rng, step_size);
  }
}

static
double
distance (void *xp_, void *yp_)
{
  lc_state_t const * restrict const xp = xp_;
  lc_state_t const * restrict const yp = yp_;
  double dist = 10 * (xp->topology != yp->topology);
  for (int d=0; d<3; ++d) {
    dist += fabs (xp->tiling[d] - yp->tiling[d]);
  }
  return dist / 2;              /* 2 = sqrt(4) */
}

static
void
print_position (void *xp_)
{
  lc_state_t const * restrict const xp = xp_;
  printf ("   %2d/{%3d,%3d,%3d}",
          xp->topology,
          xp->tiling[0], xp->tiling[1], xp->tiling[2]);
}



void
lc_auto_init (lc_statset_t * restrict const ls,
              lc_state_t * restrict const state)
{
  DECLARE_CCTK_PARAMETERS;
  assert (! cycle_j_tilings);   /* don't mix strategies */
  
  if (! ls->auto_state) {
    /* First call for this user parameter set of this loop */
    
    ls->auto_state = malloc (sizeof * ls->auto_state);
    
    /* Initialise simulated annealing parameters */
    ls->auto_state->siman_params.iters_fixed_T = siman_iters_fixed_T;
    ls->auto_state->siman_params.step_size     = siman_step_size;
    ls->auto_state->siman_params.k             = siman_k;
    ls->auto_state->siman_params.t_initial     = siman_T_initial;
    ls->auto_state->siman_params.mu_t          = siman_mu_T;    
    ls->auto_state->siman_params.t_min         = siman_T_min;
    
    /* Create random number generator */
    gsl_rng_env_setup();
    ls->auto_state->rng = gsl_rng_alloc (gsl_rng_default);
    
    /* Set initial state */
    ls->auto_state->state = * state;
    
    /* Initialise simulated annealing state */
    statset = ls;
    ls->auto_state->siman_state =
      lc_siman_solve (NULL,
                      ls->auto_state->rng,
                      NULL, 0.0,
                      take_step, distance, verbose ? print_position : NULL,
                      sizeof (lc_state_t),
                      ls->auto_state->siman_params);
    
  } else {
    /* Not the first call */
    
    if (ls->auto_state->siman_state) {
      /* The solver is still active: ask for the next state */
      statset = ls;
      ls->auto_state->siman_state =
        lc_siman_solve (ls->auto_state->siman_state,
                        ls->auto_state->rng,
                        & ls->auto_state->state, ls->auto_state->time,
                        take_step, distance, verbose ? print_position : NULL,
                        sizeof (lc_state_t),
                        ls->auto_state->siman_params);
    }
    
    /* Set thread topology and tiling specification */
    * state = ls->auto_state->state;
    
  } /* if not the first call */
}



void
lc_auto_finish (lc_statset_t * restrict const ls,
                lc_stattime_t const * restrict const lt)
{
  ls->auto_state->time =
    lt->time_calc_sum /
    (lt->time_count * ls->npoints[0] * ls->npoints[1] * ls->npoints[2]);
}
