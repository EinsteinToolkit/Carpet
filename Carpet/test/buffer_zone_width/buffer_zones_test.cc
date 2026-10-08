#include "buffer_zones.hh"

#include <cstdio>
#include <cstdlib>

static void expect(int const time_substeps, int const integrator_substeps,
                   int const ghosts, int const additional, int const want) {
  int const got = Carpet::buffer_zone_width(time_substeps, integrator_substeps,
                                            ghosts, additional);
  if (got != want) {
    std::printf("buffer_zone_width(%d, %d, %d, %d) = %d, want %d\n",
                time_substeps, integrator_substeps, ghosts, additional, got,
                want);
    std::exit(1);
  }
}

int main() {
  // One RK4 step, 3 ghosts: the width Carpet used for every level.
  expect(1, 4, 3, 0, 9);
  // Two RK4 steps (2:1 subcycling): the width that contains the error.
  expect(2, 4, 3, 0, 21);
  // Four fine steps per parent step, as in time factors [1, 4].
  expect(4, 4, 3, 0, 45);
  // User-requested extra zones are scaled with the fine steps.
  expect(1, 4, 3, 2, 11);
  expect(2, 4, 3, 2, 25);
  // Buffer zones switched off: one stencil application, no extra zones.
  expect(1, 1, 3, 0, 0);
  // RK2 with 3 ghosts, one step and a 2:1 subcycling step.
  expect(1, 2, 3, 0, 3);
  expect(2, 2, 3, 0, 9);
  return 0;
}
