#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include <cctk.h>
#include <cctk_Parameters.h>

#include "loopcontrol.h"

#include "lc_hill.h"



static inline
int
imin (int const a, int const b)
{
  return a < b ? a : b;
}

static inline
int
imax (int const a, int const b)
{
  return a > b ? a : b;
}

static
double
drand (void)
{
  return rand() / (RAND_MAX + 1.0);
}

static
int
irand (int const imaxval)
{
  return rand() / (RAND_MAX + 1.0) * imaxval;
}



static
double
time_for_stattime (lc_stattime_t const * restrict const lt)
{
  assert (lt);
  return lt->time_calc_sum / lt->time_count;
}

static
double
time_for_state (lc_statset_t const * restrict const ls,
                lc_state_t const * restrict const state)
{
  lc_stattime_t const * restrict const lt = lc_stattime_find (ls, state);
  assert (lt);
  return time_for_stattime (lt);
}



void
lc_hill_init (lc_statset_t * restrict const ls,
              lc_state_t * const new_state)
{
  DECLARE_CCTK_PARAMETERS;
  
  lc_hill_state_t * restrict lh = ls->hill_state;
  
  /* Initialise state */
  if (! lh) {
    if (verbose) {
      CCTK_INFO ("Hill climbing: Initialising");
    }
    ls->hill_state = malloc (sizeof * ls->hill_state);
    lh = ls->hill_state;
    lh->iteration = 0;
    lh->have_best = 0;
    lh->excursion_start = 0;
    lh->have_previous = 0;
    lh->state = * new_state;
    return;
  }
  
  /* If the overhead has become too large, do nothing.  */
  if (ls->time_setup_sum > maximum_setup_overhead * ls->time_calc_sum) {
    /* Stay at the old state.  */
    * new_state = lh->state;
    return;
  }
  
  ++ lh->iteration;
  
  if (verbose) {
    CCTK_VInfo (CCTK_THORNSTRING,
                "Hill climbing: iter %d, state %2d/{%2d,%2d,%2d}, time %g",
                lh->iteration,
                lh->state.topology,
                lh->state.tiling[0],
                lh->state.tiling[1],
                lh->state.tiling[2],
                lh->time);
  }
  
  /* Test whether we have a new best time */
  if (! lh->have_best || lh->time < lh->best_time) {
    /* Remember this state */
    if (verbose) {
      CCTK_INFO ("Hill climbing: This is a new best time");
    }
    lh->best = lh->state;
    lh->best_time = lh->time;
    lh->have_best = 1;
    lh->excursion_start = lh->iteration;
  } else if (lh->have_best && lc_state_equal (& lh->state, & lh->best)) {
    /* Update time for best state */
    if (verbose) {
      CCTK_INFO ("Hill climbing: Updating best time");
    }
    lh->best_time = lh->time;
  }
  
  /* Compare the time for the current state with the time for the
     previous state.  If the previous state was better, backtrack.  */
  if (lh->have_previous && lh->previous_time < lh->time) {
    if (verbose) {
      CCTK_INFO ("Hill climbing: Backtracking");
    }
    lh->state = lh->previous;
    lh->time = lh->previous_time;
    lh->have_previous = 0;
  }
    
  /* Give up if the current time is too bad */
  if (lh->have_best) {
    int const immediate_overhead = 
      lh->time > lh->best_time * (1.0 + immediate_overhead_threshold);
    int const delayed_overhead =
      lh->iteration > lh->excursion_start + overhead_threshold_delay &&
      lh->time > lh->best_time * (1.0 + delayed_overhead_threshold);
    if (immediate_overhead || delayed_overhead) {
      if (verbose) {
        CCTK_INFO ("Hill climbing: Reverting to best known state");
      }
      lh->excursion_start = lh->iteration;
      lh->state = lh->best;
      lh->time = lh->best_time;
    }
  }
  
  /* Age the state */
  lh->previous = lh->state;
  lh->previous_time = lh->time;
  lh->have_previous = 1;
  
 search:;
  
  /* Look which neighbours exist.  */
  typedef enum { nb_boundary, nb_missing, nb_exists } neighbour_t;
  neighbour_t  neighbours[3][2];
  lc_state_t   nb_state[3][2];
  double       nb_time[3][2];
  lc_state_t * nb_nonexist_state[6];
  int          num_nonexist_states = 0;
  lc_state_t * nb_minimum_time = NULL;
  double       minimum_time;
  for (int d=0; d<3; ++d) {
    for (int f=0; f<2; ++f) {
      nb_state[d][f] = lh->state;
      nb_state[d][f].tiling[d] += f ? + 1: -1;
      int const ntilings = ls->topology_ntilings[d][nb_state[d][f].topology];
      if (nb_state[d][f].tiling[d] < 0 ||
          nb_state[d][f].tiling[d] >= ntilings)
      {
        neighbours[d][f] = nb_boundary;
      } else {
        lc_stattime_t const * restrict const nb_lt =
          lc_stattime_find (ls, & nb_state[d][f]);
        if (! nb_lt) {
          neighbours[d][f] = nb_missing;
          nb_nonexist_state[num_nonexist_states++] = & nb_state[d][f];
        } else {
          neighbours[d][f] = nb_exists;
          nb_time[d][f] = time_for_stattime (nb_lt);
          if (! nb_minimum_time || nb_time[d][f] < minimum_time) {
            nb_minimum_time = & nb_state[d][f];
            minimum_time = nb_time[d][f];
          }
        }
      }
    }
  }
  
  /* If not all neighbours exist, then choose a random neighbour and
     move there.  */
  if (num_nonexist_states > 0) {
    if (verbose) {
      CCTK_INFO ("Hill climbing: Examining a new state");
    }
    int const choice = irand (num_nonexist_states);
    lh->state = * nb_nonexist_state[choice];
    * new_state = lh->state;
    return;
  }
  
  /* All neighbours exist.  Look whether we are in a local
     minimum.  */
  assert (nb_minimum_time);
  if (minimum_time >= lh->time) {
    /* We are in a local minimum.  */
    if (verbose) {
      CCTK_INFO ("Hill climbing: Local minimum reached");
    }
    
    /* Every so often take a small jump.  */
    if (drand() < probability_small_jump) {
      /* Be curious, go somewhere nearby.  */
      if (verbose) {
        CCTK_INFO ("Hill climbing: Making a small jump");
      }
      for (int ntries = 0; ntries < max_jump_attempts; ++ ntries) {
        lc_state_t try_state = lh->state;
        if (drand() < 0.25) {
          /* Change the topology, but not the tiling.  */
          try_state.topology = irand (ls->ntopologies);
          for (int d=0; d<3; ++d) {
            if (try_state.tiling[d] >=
                ls->topology_ntilings[d][try_state.topology])
            {
              /* The tiling doesn't fit for this new topology; don't
                 choose this topology.  */
              goto next_try;
            }
          }
        } else {
          /* Change the tiling a bit, but keep the topology */
          for (int d=0; d<3; ++d) {
            int const i0 =
              imax (try_state.tiling[d] - small_jump_distance, 0);
            int const i1 =
              imin (try_state.tiling[d] + small_jump_distance + 1,
                    ls->topology_ntilings[d][try_state.topology]);
            try_state.tiling[d] = i0 + irand (i1 - i0);
          }
        }
        if (! lc_stattime_find (ls, & try_state)) {
          lh->state = try_state;
          * new_state = lh->state;
          return;
        }
      next_try:;
      }
      /* Don't jump after all.  */
    }
    
    /* Every so often take a random jump.  */
    if (drand() < probability_random_jump) {
      /* Be adventurous, go somewhere unknown.  */
      if (verbose) {
        CCTK_INFO ("Hill climbing: Jumping randomly");
      }
      for (int ntries = 0; ntries < max_jump_attempts; ++ ntries) {
        lc_state_t try_state;
        try_state.topology = irand (ls->ntopologies);
        for (int d=0; d<3; ++d) {
          try_state.tiling[d] =
            irand (ls->topology_ntilings[d][try_state.topology]);
        }
        if (! lc_stattime_find (ls, & try_state)) {
          /* The new state is hitherto unknown, use it.  */
          lh->state = try_state;
          lh->excursion_start = lh->iteration;
          lh->have_previous = 0; /* disable backtracking */
          * new_state = lh->state;
          return;
        }
      }
      /* Don't jump after all.  */
    }
    
    /* If the current state is not the best state, give up and go
       back.  */
    if (! lc_state_equal (& lh->state, & lh->best)) {
      /* Revert to the best known state.  */
      if (verbose) {
        CCTK_INFO ("Hill climbing: Reverting to best known state");
      }
      lh->state = lh->best;
      lh->excursion_start = lh->iteration;
      lh->have_previous = 0;
      * new_state = lh->best;
      return;
    }
    
    /* Be content, do nothing.  */
    if (verbose) {
      CCTK_INFO ("Hill climbing: Resting");
    }
    * new_state = lh->state;
    return;
  }
  
  /* One of the neighbours is better.  Move to this neighbour, and
     continue the search there.  */
  if (verbose) {
    CCTK_INFO ("Hill climbing: Found a better neighbour, going there");
  }
  lh->state = * nb_minimum_time;
  lh->time = minimum_time;
  goto search;
}



void
lc_hill_finish (lc_statset_t * restrict const ls,
                lc_stattime_t const * restrict const lt)
{
  lc_hill_state_t * restrict const lh = ls->hill_state;
  
  lh->time = time_for_stattime (lt);
}
