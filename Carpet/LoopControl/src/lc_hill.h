#ifndef LC_HILL_H
#define LC_HILL_H

#include <cctk.h>

#include "loopcontrol.h"



/* no typedef here; forward declcared in loopcontrol.h */
struct lc_hill_state_t {
  int iteration;
  
  lc_state_t best;
  double best_time;
  int have_best;
  
  int excursion_start;          /* -1 if not on an excursion */
  
  lc_state_t previous;
  double previous_time;
  int have_previous;
  
  lc_state_t state;
  double time;
};



void
lc_hill_init (lc_statset_t * restrict const ls,
              lc_state_t * restrict const state);

void
lc_hill_finish (lc_statset_t * restrict const ls,
                lc_stattime_t const * restrict const lt);

#endif  /* #ifndef LC_HILL_H */
