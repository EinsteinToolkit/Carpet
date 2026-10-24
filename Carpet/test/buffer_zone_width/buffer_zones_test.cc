#include "buffer_zones.hh"

#include <cstdio>
#include <cstdlib>

static void expect_skip(int const iteration, int const do_every,
                       int const parent_do_every, bool const want) {
  bool const got = Carpet::skip_subcycling_prolongation(iteration, do_every,
                                                        parent_do_every);
  if (got != want) {
    std::printf("skip_subcycling_prolongation(%d, %d, %d) = %d, want %d\n",
                iteration, do_every, parent_do_every, int(got), int(want));
    std::exit(1);
  }
}

static void expect_taper(bool const time_interpolation, bool const final_copy,
                         int const iteration, int const do_every,
                         int const parent_do_every, bool const want) {
  bool const got = Carpet::taper_subcycling_step(
      time_interpolation, final_copy, iteration, do_every, parent_do_every);
  if (got != want) {
    std::printf(
        "taper_subcycling_step(%d, %d, %d, %d, %d) = %d, want %d\n",
        int(time_interpolation), int(final_copy), iteration, do_every,
        parent_do_every, int(got), int(want));
    std::exit(1);
  }
}

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
  // CSAMR Hermite keeps a buffer only for the tapered fine step.
  expect(Carpet::hermite_csamr_buffer_substeps, 4, 3, 0, 9);
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

  // 2:1 on the finest level. Strides 1 and 2. Iterations 1 and 3 are the
  // halfway steps. Iterations 2 and 4 land on the parent time.
  expect_skip(1, 1, 2, true);
  expect_skip(2, 1, 2, false);
  expect_skip(3, 1, 2, true);
  expect_skip(4, 1, 2, false);
  // 4:1. Skip the three intermediate steps. Prolongate the matching step.
  expect_skip(1, 1, 4, true);
  expect_skip(2, 1, 4, true);
  expect_skip(3, 1, 4, true);
  expect_skip(4, 1, 4, false);
  expect_skip(5, 1, 4, true);
  expect_skip(8, 1, 4, false);
  // 2:1 on a level whose own stride is 2 and whose parent stride is 4.
  expect_skip(1, 2, 4, true);
  expect_skip(3, 2, 4, false);
  expect_skip(5, 2, 4, true);
  expect_skip(7, 2, 4, false);
  // Equal time factors. This level steps only when the parent does.
  expect_skip(1, 4, 4, false);
  expect_skip(5, 4, 4, false);

  // Defaults: time interpolation and the final copy. Nothing is skipped.
  expect_taper(true, true, 1, 1, 2, false);
  expect_taper(true, true, 2, 1, 2, false);
  // Current row: no interpolation, and the final copy.
  expect_taper(false, true, 1, 1, 2, true);
  expect_taper(false, true, 2, 1, 2, false);
  expect_taper(false, true, 1, 1, 4, true);
  expect_taper(false, true, 2, 1, 4, true);
  expect_taper(false, true, 3, 1, 4, true);
  expect_taper(false, true, 4, 1, 4, false);
  // Old row and code 4: time interpolation, and no final copy.
  expect_taper(true, false, 1, 1, 2, false);
  expect_taper(true, false, 2, 1, 2, true);
  // Neither refill.
  expect_taper(false, false, 1, 1, 2, true);
  expect_taper(false, false, 2, 1, 2, true);
  // Equal time factors. The only step is the copy.
  expect_taper(true, true, 1, 4, 4, false);
  expect_taper(false, true, 1, 4, 4, false);
  expect_taper(true, false, 1, 4, 4, true);
  return 0;
}
