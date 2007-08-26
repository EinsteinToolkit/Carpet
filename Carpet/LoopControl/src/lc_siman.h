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

#ifndef SIMAN_H
#define SIMAN_H

#include <stdlib.h>
#include <gsl/gsl_rng.h>

/* types for the function pointers passed to lc_siman_solve */

typedef void (*lc_siman_step_t) (const gsl_rng *r, void *xp, double step_size);
typedef double (*lc_siman_metric_t) (void *xp, void *yp);
typedef void (*lc_siman_print_t) (void *xp);

/* this structure contains all the information needed to structure the
   search, beyond the energy function, the step function and the
   initial guess. */

typedef struct {
  int iters_fixed_T;    /* how many iterations at each temperature? */
  double step_size;     /* max step size in the random walk */
  /* the following parameters are for the Boltzmann distribution */
  double k, t_initial, mu_t, t_min;
} lc_siman_params_t;

/* prototype for the workhorse function */

typedef struct lc_siman_state_t lc_siman_state_t;

lc_siman_state_t *
lc_siman_solve (lc_siman_state_t *restrict state,
                const gsl_rng *restrict r,
                void *restrict x0_p, double E,
                lc_siman_step_t take_step,
                lc_siman_metric_t distance,
                lc_siman_print_t print_position,
                size_t element_size,
                lc_siman_params_t params);

#endif  /* #ifndef SIMAN_H */
