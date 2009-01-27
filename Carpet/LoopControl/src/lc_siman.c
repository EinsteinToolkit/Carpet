/* Simulated annealing */

/* Adapted from GSL, the GNU Scientific Library, version 1.9 */

/* Copyright (C) 1996, 1997, 1998, 1999, 2000 Mark Galassi
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gsl/gsl_machine.h>
#include <gsl/gsl_rng.h>

#include <cctk.h>

/* #ifdef HAVE_TGMATH_H */
/* #  include <tgmath.h> */
/* #endif */

#include "lc_siman.h"

static inline 
double safe_exp (double x) /* avoid underflow errors for large uphill steps */
{ 
  return (x < GSL_LOG_DBL_MIN) ? 0.0 : exp(x);
}

/* this structure contains internal state information for
   lc_siman_solve */

typedef enum { state_initial, state_first, state_looping } lc_siman_location_t;

/* no typedef here; forward declcared in lc_siman.h */
struct lc_siman_state_t {
  lc_siman_location_t state;
  const gsl_rng *r;
  lc_siman_step_t take_step;
  lc_siman_metric_t distance;
  lc_siman_print_t print_position;
  size_t element_size;
  lc_siman_params_t params;
  void *x, *new_x, *best_x;
  double E, new_E, best_E;
  int i, done;
  double T;
  int n_evals, n_iter, n_accepts, n_rejects, n_eless;
};

/* implementation of a basic simulated annealing algorithm */

lc_siman_state_t *
lc_siman_solve (lc_siman_state_t *restrict state,
                const gsl_rng *restrict r,
                void *restrict x0_p, double E0,
                lc_siman_step_t take_step,
                lc_siman_metric_t distance,
                lc_siman_print_t print_position,
                size_t element_size,
                lc_siman_params_t params)
{
  if (!state) {
    state = malloc(sizeof *state);
    state->state = state_initial;
  }
  switch (state->state) {
  case state_initial: goto label_initial;
  case state_first  : goto label_first;
  case state_looping: goto label_looping;
  }
  abort();
  
 label_initial:;
  state->n_evals = 1;
  state->n_iter = 0;
  
  state->state = state_first;
  return state;
  
 label_first:;
  state->E = E0;
  
  state->x = malloc (element_size);
  memcpy (state->x, x0_p, element_size);
  state->new_x = malloc (element_size);
  state->best_x = malloc (element_size);
  memcpy (state->best_x, x0_p, element_size);
  
  state->best_E = state->E;
  
  state->T = params.t_initial;
  state->done = 0;
  
  if (print_position) {
    printf ("#-iter  #-evals   temperature     position   energy\n");
  }
  
  /* while (!done) */
 begin_while_notdone:;
  if (state->done) goto end_while_notdone;
  
  state->n_accepts = 0;
  state->n_rejects = 0;
  state->n_eless = 0;
  
  /* for (i = 0; i < params.iters_fixed_T; ++i) */
  state->i = 0;
 begin_for_i:;
  if (state->i >= params.iters_fixed_T) goto end_for_i;
  
  memcpy (state->new_x, state->x, element_size);
  
  take_step (r, state->new_x, params.step_size);
  memcpy (x0_p, state->new_x, element_size);
  
  state->state = state_looping;
  return state;
  
 label_looping:;
  state->new_E = E0;
  
  if (state->new_E <= state->best_E) {
    memcpy (state->best_x, state->new_x, element_size);
    state->best_E = state->new_E;
  }
  
  ++state->n_evals;             /* keep track of evaluations */
  /* now take the crucial step: see if the new point is accepted or
     not, as determined by the boltzman probability */
  if (state->new_E < state->E) {
    /* yay! take a step */
    memcpy (state->x, state->new_x, element_size);
    state->E = state->new_E;
    ++state->n_eless;
  } else if (gsl_rng_uniform(r) <
             safe_exp (-(state->new_E - state->E)/(params.k * state->T)))
  {
    /* yay! take a step */
    memcpy(state->x, state->new_x, element_size);
    state->E = state->new_E;
    ++state->n_accepts;
  } else {
    ++state->n_rejects;
  }
  
  ++state->i;
  goto begin_for_i;
 end_for_i:;
  
  if (print_position) {
    /* see if we need to print stuff as we go */
    /*       printf("%5d %12g %5d %3d %3d %3d", n_iter, T, n_evals, */
    /*           100*n_eless/n_steps, 100*n_accepts/n_steps, */
    /*           100*n_rejects/n_steps); */
    printf ("%5d   %7d  %12g", state->n_iter, state->n_evals, state->T);
    print_position (state->x);
    printf ("  %12g\n", state->E);
  }
  
  /* apply the cooling schedule to the temperature */
  /* FIXME: I should also introduce a cooling schedule for the iters */
  state->T /= params.mu_t;
  ++state->n_iter;
  if (state->T < params.t_min) {
    state->done = 1;
  }
  
  goto begin_while_notdone;
 end_while_notdone:;
  
  /* at the end, copy the result onto the initial point, so we pass it
     back to the caller */
  memcpy (x0_p, state->best_x, element_size);
  
  free (state->x);
  free (state->new_x);
  free (state->best_x);
  
  free (state);
  return NULL;
}
