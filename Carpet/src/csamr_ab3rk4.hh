#ifndef CARPET_CSAMR_AB3RK4_HH
#define CARPET_CSAMR_AB3RK4_HH

// Unequal CSAMR steps. The slot order stays the Figure 7 schedule in
// csamr.hh. On every refined level the even step of a parent interval is
// Adams-Bashforth 3 and the odd step is classical RK4. The coarse level
// takes one RK4 step. AB3 is the shorter piece.
//
// ab3_parent_fraction is one quarter, not a slight trim of an equal split.
// On a factor-of-two hierarchy the finest AB3 that sits inside the long
// parent step has Courant 4*f*(1-f) times the coarse Courant, and the
// finest RK4 inside that same step has Courant 4*(1-f)*(1-f) times the
// coarse Courant. At f = 1/4 those factors are 3/4 and 9/4. Coarse Courant
// 1/9 then puts the finest AB3 at 1/12 and the finest RK4 at 1/4.

#include "csamr.hh"

namespace Carpet {
namespace csamr {

double const ab3_parent_fraction = 0.25;

inline bool ab3rk4_fraction_ok(double const fraction) {
  return fraction > 0.0 and fraction < 0.5;
}

// Offsets are measured from the start of this coarse cycle. end_offset is
// in (0, coarse_dt]. A step that finishes the cycle stores coarse_dt
// exactly, so the next cycle's origin is the same floating-point value.
// dt is the offset span. The time the integrator must use is the
// difference of the absolute end and the absolute start: subtracting the
// offset span from the absolute end is not always that start.
struct Ab3Rk4Step {
  bool ab3;
  double dt;
  double start_offset;
  double end_offset;
};

inline Ab3Rk4Step ab3rk4_step(int const level, int const ordinal,
                              double const fraction, double const coarse_dt) {
  if (not ab3rk4_fraction_ok(fraction))
    fail("AB3-RK4 fraction must lie in (0, 0.5)");
  if (not(coarse_dt > 0.0))
    fail("AB3-RK4 coarse dt is not positive");
  if (level < 0 or level >= max_levels)
    fail("AB3-RK4 level is out of range");
  if (level == 0) {
    if (ordinal != 0)
      fail("AB3-RK4 coarse level takes one step per cycle");
    Ab3Rk4Step step;
    step.ab3 = false;
    step.dt = coarse_dt;
    step.start_offset = 0.0;
    step.end_offset = coarse_dt;
    return step;
  }

  int const nsteps = time_factor(level);
  if (ordinal < 0 or ordinal >= nsteps)
    fail("AB3-RK4 step ordinal is out of range");

  // Walk the parent intervals from coarse to fine. Bit (level-1) chooses
  // the piece of the coarse step. Bit 0 chooses this step: 0 is AB3.
  double start = 0.0;
  double span = coarse_dt;
  for (int bit = level - 1; bit >= 0; --bit) {
    bool const second = ((ordinal >> bit) & 1) != 0;
    if (second) {
      start += fraction * span;
      span = (1.0 - fraction) * span;
    } else {
      span = fraction * span;
    }
  }

  Ab3Rk4Step step;
  step.ab3 = (ordinal % 2) == 0;
  step.start_offset = start;
  if (ordinal + 1 == nsteps) {
    step.dt = coarse_dt - start;
    step.end_offset = coarse_dt;
  } else {
    step.dt = span;
    step.end_offset = start + span;
  }
  if (not(step.dt > 0.0))
    fail("AB3-RK4 step dt is not positive");
  return step;
}

// The integer schedule names this level's step in the slot. Its ordinal is
// how many of this level's steps in the cycle have finished, counting this
// one as the current step. Action::level is the subcycling-group index.
// This mode only runs when each group is one carpet level, so that index
// is the carpet level.
inline int ab3rk4_ordinal(int const n_levels, int const slot_index,
                          int const level) {
  if (n_levels < 1 or n_levels > max_levels)
    fail("AB3-RK4 level count is out of range");
  std::vector<Slot> const slots = coarse_cycle(n_levels);
  if (slot_index < 0 or slot_index >= static_cast<int>(slots.size()))
    fail("AB3-RK4 slot index is out of range");
  int const max_level = n_levels - 1;
  for (Action const &action : slots[static_cast<size_t>(slot_index)].steps) {
    if (action.kind != Kind::Step or action.level != level)
      continue;
    int const units = time_factor(max_level - action.level);
    if (units < 1 or action.end_time < units or action.end_time % units != 0)
      fail("AB3-RK4 step end is not a multiple of the level dt");
    return action.end_time / units - 1;
  }
  fail("AB3-RK4 level is not stepping in this slot");
  return -1;
}

// cctk_iteration selects the Figure 7 slot: the first evolution step is
// iteration 1 and slot 0. Level 0 is the coarse RK4. A refined level is
// Adams-Bashforth 3 on its even ordinal. The caller has to pass the
// driver's current iteration. An earlier iteration selects another slot.
inline bool ab3rk4_step_is_ab3(int const n_levels, int const iteration,
                               int const level) {
  if (level <= 0)
    return false;
  if (n_levels < 2 or n_levels > max_levels)
    fail("AB3-RK4 step query level count is out of range");
  if (iteration < 1)
    fail("AB3-RK4 step query ran before the first iteration");
  int const nfine = time_factor(n_levels - 1);
  int const slot = (iteration - 1) % nfine;
  return ab3rk4_ordinal(n_levels, slot, level) % 2 == 0;
}

inline double ab3rk4_absolute_start(double const initial, int const cycle,
                                    Ab3Rk4Step const &step,
                                    double const coarse_dt) {
  if (cycle < 0)
    fail("AB3-RK4 cycle index is negative");
  if (not(coarse_dt > 0.0))
    fail("AB3-RK4 coarse dt is not positive");
  if (step.start_offset == 0.0)
    return initial + static_cast<double>(cycle) * coarse_dt;
  return initial + static_cast<double>(cycle) * coarse_dt + step.start_offset;
}

inline double ab3rk4_absolute_end(double const initial, int const cycle,
                                  Ab3Rk4Step const &step,
                                  double const coarse_dt) {
  if (cycle < 0)
    fail("AB3-RK4 cycle index is negative");
  if (not(coarse_dt > 0.0))
    fail("AB3-RK4 coarse dt is not positive");
  if (step.end_offset == coarse_dt)
    return initial + static_cast<double>(cycle + 1) * coarse_dt;
  return initial + static_cast<double>(cycle) * coarse_dt + step.end_offset;
}

} // namespace csamr
} // namespace Carpet

#endif
