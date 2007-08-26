#ifndef LC_AUTO_H
#define LC_AUTO_H

#include <cctk.h>

#include "lc_siman.h"
#include "loopcontrol.h"



typedef struct lc_auto_position_t {
  int topology;
  int tiling[3];
} lc_auto_position_t;

/* no typedef here; forward declcared in loopcontrol.h */
struct lc_auto_state_t {
  lc_siman_params_t siman_params;
  gsl_rng * rng;
  lc_siman_state_t * siman_state;
  
  lc_auto_position_t position;
  double time;
};



void
lc_auto_init (lc_statset_t * restrict const ls,
              int * restrict const topology,
              int tiling[3]);

void
lc_auto_finish (lc_statset_t * restrict const ls,
                lc_stattime_t const * restrict const lt);

#endif  /* #ifndef LC_AUTO_H */
