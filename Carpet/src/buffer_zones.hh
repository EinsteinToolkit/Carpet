#ifndef CARPET_BUFFER_ZONES_HH
#define CARPET_BUFFER_ZONES_HH

#include <cassert>

namespace Carpet {

// Buffer zones on one face of one refinement level, excluding ghosts.
//
// A fine level takes time_substeps steps for each step of its parent.
// Each of those steps evaluates the RHS integrator_substeps times, and
// each evaluation reaches ghosts points in from the boundary.
// Prolongation fills the ghosts and the buffer only while the parent is
// active. The buffer has to absorb every later RHS, so its width is
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

} // namespace Carpet

#endif
