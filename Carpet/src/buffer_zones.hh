#ifndef CARPET_BUFFER_ZONES_HH
#define CARPET_BUFFER_ZONES_HH

#include <cassert>

namespace Carpet {

// Buffer zones on one face of one refinement level, excluding ghosts.
//
// A fine level takes time_substeps steps for each step of its parent.
// Each of those steps evaluates the RHS integrator_substeps times, and
// each evaluation reaches ghosts points in from the boundary.
// use_tapered_grids sets time_substeps to the number of fine steps in the
// parent interval, so the buffer can absorb every one of those steps. Its
// width is
//
//   time_substeps * (integrator_substeps * ghosts + additional) - ghosts
//
// One RK4 step with 3 ghosts is 9 zones. Two RK4 steps, the usual 2:1
// subcycling ratio, is 21 zones.
inline int buffer_zone_width(int const time_substeps,
                             int const integrator_substeps, int const ghosts,
                             int const additional) {
  assert(time_substeps >= 1);
  assert(integrator_substeps >= 1);
  assert(ghosts >= 0);
  int const width =
      time_substeps * (integrator_substeps * ghosts + additional) - ghosts;
  assert(width >= 0);
  return width;
}


// True when this fine step must not prolongate its boundary.
//
// iteration is cctk_iteration after AdvanceTime. do_every and
// parent_do_every are the iteration strides of this level and its parent.
// Carpet advances the parent first. The early fine steps of a parent
// interval land strictly between parent times, and prolongation there is
// time interpolation. The last fine step of the interval lands on the
// parent time, so prolongation copies.
inline bool skip_subcycling_prolongation(int const iteration,
                                         int const do_every,
                                         int const parent_do_every) {
  assert(iteration >= 1);
  assert(do_every >= 1);
  assert(parent_do_every >= do_every);
  assert(parent_do_every % do_every == 0);
  int const phase = (iteration - 1) % parent_do_every;
  int const coincident_phase = parent_do_every - do_every;
  return phase != coincident_phase;
}

// Whether this fine step skips boundary prolongation.
//
// A step between parent times is skipped when time_interpolation is false.
// The step on the parent time is skipped when final_copy is false.
inline bool taper_subcycling_step(bool const time_interpolation,
                                  bool const final_copy, int const iteration,
                                  int const do_every,
                                  int const parent_do_every) {
  if (skip_subcycling_prolongation(iteration, do_every, parent_do_every))
    return not time_interpolation;
  return not final_copy;
}

} // namespace Carpet

#endif
