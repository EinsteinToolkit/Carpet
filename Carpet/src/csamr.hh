#ifndef CARPET_CSAMR_HH
#define CARPET_CSAMR_HH

// Concurrent structured AMR schedule.
//
// Löffler, Cao, Brandt, and Du, "A new parallelization scheme for
// adaptive mesh refinement", Journal of Computational Science 16
// (2016) 79–88, Figure 7. The refinement factor is 2. Level 0 is the
// coarsest level and level n_levels-1 is the finest.
//
// One coarse cycle is a list of fine-step slots. A slot with two
// steps is the parallel pair in Figure 7: the finest level and one
// coarser level, which do not depend on each other. Prolongations in
// "before" run before that slot's steps and after the previous
// slot's restrictions.

#include <cstdio>
#include <cstdlib>
#include <vector>

namespace Carpet {
namespace csamr {

// Deepest hierarchy this schedule will build. 2^(max_levels-1) fits
// in a signed 32-bit fine-step counter.
int const max_levels = 20;

enum class Kind { Step, Restrict, Prolong, ProlongTimeInterpolate };

struct Action {
  Kind kind;
  // Step: the level that advances. Restrict: the coarser destination.
  // Prolong: the coarser source.
  int level;
  // Step: the other level of a parallel pair, or -1. Restrict: the
  // finer source. Prolong: the finer destination.
  int other;
  // Step: time at the end of this step, in units of the finest dt,
  // measured from the start of the coarse cycle. Otherwise 0.
  int end_time;
  bool in_parallel;

  static Action step(int level, int partner) {
    Action action;
    action.kind = Kind::Step;
    action.level = level;
    action.other = partner;
    action.end_time = 0;
    action.in_parallel = partner >= 0;
    return action;
  }

  static Action restrict_onto(int destination, int source) {
    Action action;
    action.kind = Kind::Restrict;
    action.level = destination;
    action.other = source;
    action.end_time = 0;
    action.in_parallel = false;
    return action;
  }

  static Action prolong(int source, int destination, bool time_interpolate) {
    Action action;
    action.kind = time_interpolate ? Kind::ProlongTimeInterpolate : Kind::Prolong;
    action.level = source;
    action.other = destination;
    action.end_time = 0;
    action.in_parallel = false;
    return action;
  }
};

struct Slot {
  std::vector<Action> before;
  std::vector<Action> steps;
  std::vector<Action> restricts;
};

inline void fail(char const *message) {
  std::fprintf(stderr, "CSAMR: %s\n", message);
  std::abort();
}

// Factor by which this level's timestep divides the coarsest timestep.
// The paper's refinement factor is two, so the factor is 2^level. This
// is the only place that encodes that two.
inline int time_factor(int const level) {
  if (level < 0 or level >= max_levels)
    fail("time refinement level is out of range");
  return 1 << level;
}

// One run of equal timesteps. first and last are carpet refinement
// levels, inclusive. factor is the time refinement factor.
struct Group {
  int factor;
  int first;
  int last;
};

enum class FactorError { Ok, Empty, TooDeep, Coarsest, Jump };

// Carpet levels may share a timestep. Each new timestep is exactly
// twice the previous one, so the distinct factors are 1, 2, 4, ...
inline FactorError classify_factors(std::vector<int> const &factors) {
  if (factors.empty())
    return FactorError::Empty;
  if (factors[0] != 1)
    return FactorError::Coarsest;
  int groups = 1;
  int previous = factors[0];
  for (size_t i = 1; i < factors.size(); ++i) {
    int const factor = factors[static_cast<int>(i)];
    if (factor == previous)
      continue;
    if (factor != previous * 2)
      return FactorError::Jump;
    ++groups;
    previous = factor;
    if (groups > max_levels)
      return FactorError::TooDeep;
  }
  return FactorError::Ok;
}

inline std::vector<Group> subcycling_groups(std::vector<int> const &factors) {
  if (classify_factors(factors) != FactorError::Ok)
    fail("time refinement factors are not a 1,2,4,... hierarchy");
  std::vector<Group> groups;
  Group group;
  group.factor = factors[0];
  group.first = 0;
  group.last = 0;
  for (int rl = 1; rl < static_cast<int>(factors.size()); ++rl) {
    if (factors[rl] == group.factor) {
      group.last = rl;
      continue;
    }
    groups.push_back(group);
    group.factor = factors[rl];
    group.first = rl;
    group.last = rl;
  }
  groups.push_back(group);
  return groups;
}

inline int group_of(std::vector<Group> const &groups, int const level) {
  for (int i = 0; i < static_cast<int>(groups.size()); ++i)
    if (level >= groups[i].first and level <= groups[i].last)
      return i;
  fail("carpet level is not in a subcycling group");
  return -1;
}

// Figure 7, short_recursive_step. Appends the body, not the initial
// prolongation or the final restriction.
inline void short_recursive_step(int const level, int const max_level,
                                 std::vector<Action> &ops) {
  if (level < max_level) {
    short_recursive_step(level + 1, max_level, ops);
    // Lines 8–9 of Figure 7. Finest step first: it must not observe
    // the coarser step's update.
    ops.push_back(Action::step(max_level, level));
    ops.push_back(Action::step(level, max_level));
    if (max_level - 1 > level) {
      for (int lr = max_level - 1; lr >= level + 1; --lr)
        ops.push_back(Action::restrict_onto(lr, lr + 1));
      ops.push_back(Action::prolong(level, level + 1, true));
      for (int lp = level + 1; lp <= max_level - 1; ++lp)
        ops.push_back(Action::prolong(lp, lp + 1, false));
      short_recursive_step(level + 1, max_level, ops);
    }
  }
}

// One coarse step. The returned vector has length 2^(n_levels-1),
// one entry per finest-level step.
inline std::vector<Slot> coarse_cycle(int const n_levels) {
  if (n_levels < 1 or n_levels > max_levels)
    fail("refinement level count is out of range");

  int const max_level = n_levels - 1;
  std::vector<Action> ops;
  // Lines 19–20, before recursive_step(0).
  for (int lp = 0; lp <= max_level - 1; ++lp)
    ops.push_back(Action::prolong(lp, lp + 1, false));
  // recursive_step: short body, then the trailing finest step.
  short_recursive_step(0, max_level, ops);
  ops.push_back(Action::step(max_level, -1));
  // Lines 22–23.
  for (int lr = max_level - 1; lr >= 0; --lr)
    ops.push_back(Action::restrict_onto(lr, lr + 1));

  std::vector<int> time(n_levels, 0);
  for (Action &op : ops) {
    if (op.kind != Kind::Step)
      continue;
    int const dt = time_factor(max_level - op.level);
    time[op.level] += dt;
    op.end_time = time[op.level];
  }

  struct Raw {
    std::vector<Action> steps;
    std::vector<Action> trailing;
  };
  std::vector<Action> leading;
  std::vector<Raw> raws;
  for (Action const &op : ops) {
    bool const starts_slot = op.kind == Kind::Step and op.level == max_level;
    if (starts_slot) {
      if (raws.empty() or not raws.back().steps.empty())
        raws.push_back(Raw());
      raws.back().steps.push_back(op);
    } else if (op.kind == Kind::Step) {
      if (raws.empty() or raws.back().steps.empty())
        fail("coarser step is not paired with a finest step");
      raws.back().steps.push_back(op);
    } else if (raws.empty()) {
      leading.push_back(op);
    } else {
      raws.back().trailing.push_back(op);
    }
  }

  std::vector<Action> carry = leading;
  std::vector<Slot> slots;
  for (Raw const &raw : raws) {
    Slot slot;
    slot.before = carry;
    carry.clear();
    slot.steps = raw.steps;
    for (Action const &trailing : raw.trailing) {
      if (trailing.kind == Kind::Restrict)
        slot.restricts.push_back(trailing);
      else
        carry.push_back(trailing);
    }
    slots.push_back(slot);
  }
  if (not carry.empty())
    fail("prolongation follows the last step of the coarse cycle");

  int const nfine = time_factor(max_level);
  if (static_cast<int>(slots.size()) != nfine)
    fail("coarse cycle does not contain one slot per finest step");

  for (Slot const &slot : slots) {
    if (slot.steps.empty() or slot.steps.size() > 2)
      fail("a slot must take one step or one parallel pair");
    if (slot.steps[0].level != max_level)
      fail("every slot starts with the finest level");
    if (slot.steps.size() == 2) {
      if (not slot.steps[0].in_parallel or not slot.steps[1].in_parallel)
        fail("a pair of steps must be marked parallel");
      if (slot.steps[1].level >= max_level)
        fail("the parallel partner must be coarser than the finest level");
      if (slot.steps[0].other != slot.steps[1].level or
          slot.steps[1].other != slot.steps[0].level)
        fail("parallel partners do not name each other");
    } else if (n_levels > 1 and slot.steps[0].in_parallel) {
      fail("a single finest step must not be marked parallel");
    }
  }
  return slots;
}

// Whether a fine step skips boundary prolongation.
// end_index is this step's fine-step index. parent_index is the parent's
// index already visible. On a CSAMR slot that is the index from before
// the slot applies the parent's own step, and the finest level is
// planned first, so a child sees its parent as it was. A child that
// finishes ahead of the parent tapers. A child that lands on the parent
// copies when final_copy is set. A child that finishes behind the parent
// time-interpolates when that refill is set.
inline bool taper_step(int const end_index, int const parent_index,
                       bool const time_interpolation, bool const final_copy) {
  if (end_index > parent_index)
    return true;
  if (end_index == parent_index)
    return not final_copy;
  return not time_interpolation;
}

// Stock Adams-Bashforth 3 evaluates one right-hand side at the step
// start. The ghost zones cover that evaluation, so the step does not
// use buffer zones. The boundary sync runs at the step end, after MoL
// has restored that end time. parent_before is the parent's fine-step
// index at the start of the slot. parent_steps is true when this slot
// updates the parent. The sync may run during the step only when the
// parent stays idle and already stores that end. Otherwise it waits
// until the parent has reached the end.
inline bool ab3_prolongate_during_step(int const end_index,
                                      int const parent_before,
                                      bool const parent_steps) {
  if (parent_steps)
    return false;
  return parent_before >= end_index;
}

// partner_rl is the carpet level evolved in parallel with rl, or -1
// when rl steps alone. True when that partner is the parent.
inline bool ab3_fill_before_partner(int const rl, int const partner_rl) {
  return rl > 0 and partner_rl == rl - 1;
}

// Fine-step index of every carpet level at the start of this slot.
// Indices are counted from the initial time. Levels in one group share
// an index. slot_index counts from the start of this coarse cycle.
inline std::vector<int>
fine_index_before_slot(std::vector<Group> const &groups,
                       std::vector<Slot> const &slots, int const slot_index,
                       int const cycle) {
  if (groups.empty())
    fail("subcycling hierarchy is empty");
  int const nfine = groups.back().factor;
  int const ncarpet = groups.back().last + 1;
  if (slot_index < 0 or slot_index > static_cast<int>(slots.size()))
    fail("slot index is out of range");
  std::vector<int> at(ncarpet, cycle * nfine);
  for (int s = 0; s < slot_index; ++s) {
    for (Action const &step : slots[s].steps) {
      if (step.level < 0 or step.level >= static_cast<int>(groups.size()))
        fail("step names an unknown subcycling group");
      int const end = cycle * nfine + step.end_time;
      Group const &group = groups[step.level];
      for (int rl = group.first; rl <= group.last; ++rl)
        at[rl] = end;
    }
  }
  return at;
}

// How many OpenMP threads each update of a parallel pair receives.
// parallel is false when the process has fewer than two threads: the
// two updates then run one after the other on that single thread.
// An odd count gives the leftover thread to the coarser update, which
// may step more than one level while the finest update steps one.
struct ThreadTeams {
  bool parallel;
  int fine;
  int coarse;
};

inline ThreadTeams thread_teams(int const nthreads) {
  ThreadTeams teams;
  int const n = nthreads < 1 ? 1 : nthreads;
  teams.parallel = n >= 2;
  teams.fine = teams.parallel ? n / 2 : n;
  teams.coarse = teams.parallel ? n - teams.fine : n;
  return teams;
}

} // namespace csamr
} // namespace Carpet

#endif // CARPET_CSAMR_HH
