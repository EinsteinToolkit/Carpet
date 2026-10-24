// Standalone checks for the CSAMR schedule in csamr.hh.
// Build: g++ -std=c++11 -Wall -Wextra -Werror -o csamr_test csamr_test.cc

#include "csamr.hh"

#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <string>
#include <vector>

namespace {

int failures = 0;

void check(bool const ok, char const *const text, int const line) {
  if (ok)
    return;
  std::fprintf(stderr, "FAIL %s:%d: %s\n", "csamr_test.cc", line, text);
  ++failures;
}

#define CHECK(cond) check(static_cast<bool>(cond), #cond, __LINE__)

using Carpet::csamr::Action;
using Carpet::csamr::Kind;
using Carpet::csamr::Slot;

char const *kind_name(Kind const kind) {
  switch (kind) {
  case Kind::Step:
    return "step";
  case Kind::Restrict:
    return "restrict";
  case Kind::Prolong:
    return "prolong";
  case Kind::ProlongTimeInterpolate:
    return "prolong-ti";
  }
  return "unknown";
}

void print_slots(std::vector<Slot> const &slots) {
  for (int s = 0; s < static_cast<int>(slots.size()); ++s) {
    std::fprintf(stderr, "  slot %d\n", s);
    Slot const &slot = slots[s];
    for (Action const &action : slot.before)
      std::fprintf(stderr, "    before %s %d -> %d\n", kind_name(action.kind),
                   action.level, action.other);
    for (Action const &action : slot.steps)
      std::fprintf(stderr, "    step level %d end %d parallel %d\n",
                   action.level, action.end_time, action.other);
    for (Action const &action : slot.restricts)
      std::fprintf(stderr, "    restrict %d onto %d\n", action.other,
                   action.level);
  }
}

bool same_transfer(Action const &action, Kind const kind, int const source,
                   int const destination) {
  return action.kind == kind and action.level == source and
         action.other == destination;
}

bool same_step(Action const &action, int const level, int const end_time,
               int const partner) {
  return action.kind == Kind::Step and action.level == level and
         action.end_time == end_time and action.other == partner and
         action.in_parallel == (partner >= 0);
}

void check_fixed_schedules() {
  std::vector<Slot> const two = Carpet::csamr::coarse_cycle(2);
  CHECK(two.size() == 2);
  if (two.size() == 2) {
    CHECK(two[0].before.size() == 1);
    CHECK(same_transfer(two[0].before[0], Kind::Prolong, 0, 1));
    CHECK(two[0].steps.size() == 2);
    CHECK(same_step(two[0].steps[0], 1, 1, 0));
    CHECK(same_step(two[0].steps[1], 0, 2, 1));
    CHECK(two[0].restricts.empty());
    CHECK(two[1].before.empty());
    CHECK(two[1].steps.size() == 1);
    CHECK(same_step(two[1].steps[0], 1, 2, -1));
    CHECK(two[1].restricts.size() == 1);
    CHECK(same_transfer(two[1].restricts[0], Kind::Restrict, 0, 1));
  } else {
    print_slots(two);
  }

  // Figure 8, right panel, is two coarse cycles of a 3-level run.
  std::vector<Slot> const three = Carpet::csamr::coarse_cycle(3);
  CHECK(three.size() == 4);
  if (three.size() == 4) {
    CHECK(three[0].before.size() == 2);
    CHECK(same_transfer(three[0].before[0], Kind::Prolong, 0, 1));
    CHECK(same_transfer(three[0].before[1], Kind::Prolong, 1, 2));
    CHECK(same_step(three[0].steps[0], 2, 1, 1));
    CHECK(same_step(three[0].steps[1], 1, 2, 2));
    CHECK(three[0].restricts.empty());

    CHECK(three[1].before.empty());
    CHECK(same_step(three[1].steps[0], 2, 2, 0));
    CHECK(same_step(three[1].steps[1], 0, 4, 2));
    CHECK(three[1].restricts.size() == 1);
    CHECK(same_transfer(three[1].restricts[0], Kind::Restrict, 1, 2));

    CHECK(three[2].before.size() == 2);
    CHECK(same_transfer(three[2].before[0], Kind::ProlongTimeInterpolate, 0, 1));
    CHECK(same_transfer(three[2].before[1], Kind::Prolong, 1, 2));
    CHECK(same_step(three[2].steps[0], 2, 3, 1));
    CHECK(same_step(three[2].steps[1], 1, 4, 2));
    CHECK(three[2].restricts.empty());

    CHECK(three[3].before.empty());
    CHECK(same_step(three[3].steps[0], 2, 4, -1));
    CHECK(three[3].restricts.size() == 2);
    CHECK(same_transfer(three[3].restricts[0], Kind::Restrict, 1, 2));
    CHECK(same_transfer(three[3].restricts[1], Kind::Restrict, 0, 1));
  } else {
    print_slots(three);
  }
}

// Figure 8. Numbers are the time after the step, in finest-dt units.
// A zero on the coarse side means that process group is idle.
void check_figure_8() {
  std::vector<Slot> const two = Carpet::csamr::coarse_cycle(2);
  int const coarse_two[] = {2, 0, 4, 0};
  int const fine_two[] = {1, 2, 3, 4};
  for (int cycle = 0; cycle < 2; ++cycle) {
    for (int slot = 0; slot < 2; ++slot) {
      int const index = cycle * 2 + slot;
      int coarse_time = 0;
      int fine_time = 0;
      for (Action const &step : two[slot].steps) {
        int const time = cycle * 2 + step.end_time;
        if (step.level == 0)
          coarse_time = time;
        if (step.level == 1)
          fine_time = time;
      }
      CHECK(coarse_time == coarse_two[index]);
      CHECK(fine_time == fine_two[index]);
    }
  }

  std::vector<Slot> const three = Carpet::csamr::coarse_cycle(3);
  // Left column of the 3-level panel: the time of whichever coarser
  // level this slot steps. Black boxes in the figure are level 0.
  int const coarse_group_three[] = {2, 4, 4, 0, 6, 8, 8, 0};
  int const coarse_group_level[] = {1, 0, 1, -1, 1, 0, 1, -1};
  int const fine_three[] = {1, 2, 3, 4, 5, 6, 7, 8};
  for (int cycle = 0; cycle < 2; ++cycle) {
    for (int slot = 0; slot < 4; ++slot) {
      int const index = cycle * 4 + slot;
      int coarse_time = 0;
      int coarse_level = -1;
      int fine_time = 0;
      for (Action const &step : three[slot].steps) {
        int const time = cycle * 4 + step.end_time;
        if (step.level == 2)
          fine_time = time;
        else {
          coarse_time = time;
          coarse_level = step.level;
        }
      }
      CHECK(fine_time == fine_three[index]);
      CHECK(coarse_time == coarse_group_three[index]);
      CHECK(coarse_level == coarse_group_level[index]);
    }
  }
}

void check_invariants(int const n_levels) {
  std::vector<Slot> const slots = Carpet::csamr::coarse_cycle(n_levels);
  int const max_level = n_levels - 1;
  int const nfine = 1 << max_level;
  CHECK(static_cast<int>(slots.size()) == nfine);

  std::vector<int> counts(n_levels, 0);
  std::vector<int> time(n_levels, 0);
  int idle_slots = 0;
  for (int slot_index = 0; slot_index < nfine; ++slot_index) {
    Slot const &slot = slots[slot_index];
    bool coarse_stepped = false;
    for (Action const &step : slot.steps) {
      CHECK(step.level >= 0 and step.level < n_levels);
      int const dt = 1 << (max_level - step.level);
      time[step.level] += dt;
      CHECK(step.end_time == time[step.level]);
      ++counts[step.level];
      if (step.level != max_level)
        coarse_stepped = true;
      if (step.level > 0) {
        // Finest-first execution: the parent has not stepped in this
        // slot yet, so the recorded parent time is the pre-step value
        // only when the parent is later in the slot. Read it before
        // applying this step by using time[] after the increment of
        // earlier (finer) steps only. The parent is coarser, so it
        // has not been updated yet.
        CHECK(time[step.level - 1] >= time[step.level] - dt);
      }
    }
    if (not coarse_stepped)
      ++idle_slots;

    for (Action const &restriction : slot.restricts) {
      CHECK(restriction.other == restriction.level + 1);
      CHECK(time[restriction.other] == time[restriction.level]);
      // Carpet restricts onto level L when
      // iteration % 2^(max_level-L) == 0. iteration = slot_index+1
      // inside the first coarse cycle.
      int const do_every = 1 << (max_level - restriction.level);
      CHECK(((slot_index + 1) % do_every) == 0);
    }
  }

  for (int level = 0; level < n_levels; ++level) {
    // Level L advances by 2^(max_level-L) fine steps, so it takes
    // 2^L steps to cross one coarse step.
    CHECK(counts[level] == (1 << level));
    CHECK(time[level] == nfine);
  }
  // The coarse process group is idle for the trailing finest step
  // only. That is 1/2^(n_levels-1) of the slots.
  if (n_levels > 1)
    CHECK(idle_slots == 1);

  // Every Carpet restrict iteration is present, not only the other way.
  for (int level = 0; level < max_level; ++level) {
    int const do_every = 1 << (max_level - level);
    for (int slot_index = 0; slot_index < nfine; ++slot_index) {
      if (((slot_index + 1) % do_every) != 0)
        continue;
      bool found = false;
      for (Action const &restriction : slots[slot_index].restricts)
        if (restriction.level == level)
          found = true;
      CHECK(found);
    }
  }
}

// Adjacent prolongation has to see a parent time at or ahead of the
// child, and a time-interpolated prolong has to have both sides of
// the child time still stored.
void check_prolongation_history(int const n_levels, int const capacity) {
  std::vector<Slot> const slots = Carpet::csamr::coarse_cycle(n_levels);
  int const max_level = n_levels - 1;
  std::vector<int> time(n_levels, 0);
  std::vector<std::vector<int> > history(n_levels);
  for (int level = 0; level < n_levels; ++level)
    history[level].push_back(0);

  for (Slot const &slot : slots) {
    for (Action const &action : slot.before) {
      CHECK(action.kind == Kind::Prolong or
            action.kind == Kind::ProlongTimeInterpolate);
      int const source = action.level;
      int const destination = action.other;
      CHECK(destination == source + 1);
      CHECK(time[source] >= time[destination]);
      bool has_earlier_or_equal = false;
      bool has_later_or_equal = false;
      for (int stored : history[source]) {
        if (stored <= time[destination])
          has_earlier_or_equal = true;
        if (stored >= time[destination])
          has_later_or_equal = true;
      }
      CHECK(has_earlier_or_equal);
      CHECK(has_later_or_equal);
    }
    for (Action const &step : slot.steps) {
      int const dt = 1 << (max_level - step.level);
      time[step.level] += dt;
      history[step.level].insert(history[step.level].begin(), time[step.level]);
      if (static_cast<int>(history[step.level].size()) > capacity)
        history[step.level].pop_back();
    }
    for (Action const &restriction : slot.restricts)
      CHECK(time[restriction.other] == time[restriction.level]);
  }
}

} // namespace

void check_grouped_factors() {
  using Carpet::csamr::FactorError;
  std::vector<int> factors(7, 1);
  factors[5] = 2;
  factors[6] = 4;
  CHECK(Carpet::csamr::classify_factors(factors) == FactorError::Ok);
  std::vector<Carpet::csamr::Group> const groups =
      Carpet::csamr::subcycling_groups(factors);
  CHECK(groups.size() == 3);
  CHECK(groups[0].first == 0 and groups[0].last == 4 and groups[0].factor == 1);
  CHECK(groups[1].first == 5 and groups[1].last == 5 and groups[1].factor == 2);
  CHECK(groups[2].first == 6 and groups[2].last == 6 and groups[2].factor == 4);
  CHECK(Carpet::csamr::group_of(groups, 4) == 0);
  CHECK(Carpet::csamr::group_of(groups, 5) == 1);
  CHECK(Carpet::csamr::group_of(groups, 6) == 2);

  std::vector<int> jump(2);
  jump[0] = 1;
  jump[1] = 4;
  CHECK(Carpet::csamr::classify_factors(jump) == FactorError::Jump);
  std::vector<int> coarsest(1, 2);
  CHECK(Carpet::csamr::classify_factors(coarsest) == FactorError::Coarsest);

  // After slot 0 of a 1,2,4 hierarchy the fine level is at 1, the
  // middle level is at 2, and every equal-timestep coarse level is
  // still at 0. That is the charged-tapered time-factor pattern.
  std::vector<Carpet::csamr::Slot> const slots = Carpet::csamr::coarse_cycle(3);
  std::vector<int> const at =
      Carpet::csamr::fine_index_before_slot(groups, slots, 1, 0);
  CHECK(at.size() == 7);
  CHECK(at[0] == 0 and at[4] == 0);
  CHECK(at[5] == 2);
  CHECK(at[6] == 1);
}

// Brill factors [1, 2, 4], use_final_copy and no time interpolation.
// The first fine step of each parent interval finishes ahead of the
// parent and keeps the buffer. The second lands on the parent time
// and is the Hermite step. One coarse cycle has three of each.
void check_hybrid_refill() {
  std::vector<int> factors(3);
  factors[0] = 1;
  factors[1] = 2;
  factors[2] = 4;
  std::vector<Carpet::csamr::Group> const groups =
      Carpet::csamr::subcycling_groups(factors);
  std::vector<Slot> const slots = Carpet::csamr::coarse_cycle(3);
  CHECK(groups.size() == 3);
  CHECK(slots.size() == 4);
  char const *const want[] = {
      "0 L2 buffer", "0 L1 buffer", "1 L2 hermite",
      "2 L2 buffer", "2 L1 hermite", "3 L2 hermite",
  };
  int const nwant = 6;
  for (int cycle = 0; cycle < 2; ++cycle) {
    int seen = 0;
    int buffers = 0;
    int hermite = 0;
    int const nfine = groups.back().factor;
    for (int slot_index = 0; slot_index < static_cast<int>(slots.size());
         ++slot_index) {
      std::vector<int> at = Carpet::csamr::fine_index_before_slot(
          groups, slots, slot_index, cycle);
      for (Action const &step : slots[slot_index].steps) {
        Carpet::csamr::Group const &group = groups[step.level];
        int const fine_index = cycle * nfine + step.end_time;
        for (int rl = group.first; rl <= group.last; ++rl) {
          if (rl > 0) {
            bool const taper = Carpet::csamr::taper_step(
                fine_index, at[rl - 1], false, true);
            std::ostringstream label;
            label << slot_index << " L" << rl
                  << (taper ? " buffer" : " hermite");
            CHECK(seen < nwant);
            if (seen < nwant and label.str() != want[seen]) {
              std::fprintf(stderr,
                           "hybrid refill cycle %d got \"%s\", want \"%s\"\n",
                           cycle, label.str().c_str(), want[seen]);
            }
            CHECK(seen >= nwant or label.str() == want[seen]);
            if (taper)
              ++buffers;
            else
              ++hermite;
            ++seen;
          }
          at[rl] = fine_index;
        }
      }
    }
    CHECK(seen == nwant);
    CHECK(buffers == 3);
    CHECK(hermite == 3);
  }

  CHECK(Carpet::csamr::taper_step(3, 2, false, true));
  CHECK(not Carpet::csamr::taper_step(4, 4, false, true));
  CHECK(Carpet::csamr::taper_step(4, 4, false, false));
  CHECK(not Carpet::csamr::taper_step(2, 4, true, true));
  CHECK(Carpet::csamr::taper_step(2, 4, false, true));
}

void check_thread_teams() {
  Carpet::csamr::ThreadTeams const none = Carpet::csamr::thread_teams(0);
  CHECK(not none.parallel);
  CHECK(none.fine == 1 and none.coarse == 1);
  Carpet::csamr::ThreadTeams const negative = Carpet::csamr::thread_teams(-4);
  CHECK(not negative.parallel);
  CHECK(negative.fine == 1 and negative.coarse == 1);

  Carpet::csamr::ThreadTeams const one = Carpet::csamr::thread_teams(1);
  CHECK(not one.parallel);
  CHECK(one.fine == 1 and one.coarse == 1);

  Carpet::csamr::ThreadTeams const two = Carpet::csamr::thread_teams(2);
  CHECK(two.parallel);
  CHECK(two.fine == 1 and two.coarse == 1);

  // The leftover thread goes to the coarser update.
  Carpet::csamr::ThreadTeams const three = Carpet::csamr::thread_teams(3);
  CHECK(three.parallel);
  CHECK(three.fine == 1 and three.coarse == 2);

  Carpet::csamr::ThreadTeams const eight = Carpet::csamr::thread_teams(8);
  CHECK(eight.parallel);
  CHECK(eight.fine == 4 and eight.coarse == 4);
  CHECK(eight.fine + eight.coarse == 8);

  Carpet::csamr::ThreadTeams const nine = Carpet::csamr::thread_teams(9);
  CHECK(nine.parallel);
  CHECK(nine.fine == 4 and nine.coarse == 5);
  CHECK(nine.fine + nine.coarse == 9);
}

int main() {
  for (int level = 0; level < 8; ++level)
    CHECK(Carpet::csamr::time_factor(level) == (1 << level));
  check_grouped_factors();
  check_hybrid_refill();
  check_thread_teams();
  check_fixed_schedules();
  check_figure_8();
  for (int n_levels = 1; n_levels <= 8; ++n_levels)
    check_invariants(n_levels);
  for (int n_levels = 2; n_levels <= 8; ++n_levels) {
    // Two stored times are enough for an adjacent linear interpolant
    // because a time-interpolated prolong follows the parent step
    // that first reaches past the child.
    check_prolongation_history(n_levels, 2);
  }
  if (failures != 0) {
    std::fprintf(stderr, "%d checks failed\n", failures);
    return 1;
  }
  std::printf("csamr schedule checks passed\n");
  return 0;
}
