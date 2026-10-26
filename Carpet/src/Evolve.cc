#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>

#include <cctk.h>
#include <cctk_Parameter.h>
#include <cctk_Parameters.h>
#include <cctk_Termination.h>
#include <cctk_Timers.h>
#include <util_String.h>

#include <Requirements.hh>

#include <CactusTimerSet.hh>
#include <Timer.hh>

#include <dist.hh>
#include <th.hh>

#include <carpet.hh>
#include "csamr.hh"
#include "csamr_ab3rk4.hh"
#include "functions.hh"

// C linkage, so these names are MoL's globals. A weak declaration stays
// unresolved when MoL is not linked. Hermite.cc declares the same symbols.
extern "C" {
extern CCTK_INT *EvolvedVariableIndex __attribute__((weak));
extern CCTK_INT MoLNumEvolvedVariables __attribute__((weak));
extern CCTK_INT MoLNumEvolvedArrayVariables __attribute__((weak));
}

#ifdef _OPENMP
#include <omp.h>
#endif

#include "buffer_zones.hh"

namespace Carpet {

using namespace std;

static bool do_terminate(const cGH *cctkGH);

// Choose which fine steps prolongate.
// use_time_interpolation refills the ghost and buffer zones on fine steps
// that fall between parent times. use_final_copy refills them on the fine
// step that lands on the parent time. use_tapered_grids selects the buffer
// width and does not choose either refill.
static void SetTaperForSubcycling(cGH const *const cctkGH) {
  DECLARE_CCTK_PARAMETERS;
  do_taper = false;
  // Serial Hermite prolongates every fine stage, so this level never
  // tapers. CSAMR still tapers the fine step that runs before its parent.
  if (hermite_dense_output and not use_csamr)
    return;
  if (reflevel == 0)
    return;
  if (use_time_interpolation and use_final_copy)
    return;
  int const stride = ipow(mgfact, mglevel);
  int const do_every =
      stride * (maxtimereflevelfact / timereffacts.AT(reflevel));
  int const parent_do_every =
      stride * (maxtimereflevelfact / timereffacts.AT(reflevel - 1));
  if (do_every < 1 or parent_do_every < do_every or
      parent_do_every % do_every != 0) {
    CCTK_VERROR("Refinement level %d has stride %d and its parent has "
                "stride %d. Subcycling requires the parent stride to be "
                "a multiple of this level's stride.",
                reflevel, do_every, parent_do_every);
  }
  do_taper = taper_subcycling_step(use_time_interpolation, use_final_copy,
                                   cctkGH->cctk_iteration, do_every,
                                   parent_do_every);
}

static void AdvanceTime(cGH *cctkGH);
static void CallRegrid(cGH *cctkGH);
static void CallEvol(cGH *cctkGH);
static void CallEvolCsamr(cGH *cctkGH);
static void CallRestrict(cGH *cctkGH);
static void CallAnalysis(cGH *cctkGH);

static void print_internal_data();

static void ScheduleTraverse(char const *where, char const *name, cGH *cctkGH);
static void OutputGH(char const *where, cGH *cctkGH);

int Evolve(tFleshConfig *const fc) {
  DECLARE_CCTK_PARAMETERS;

  Waypoint("Starting evolution loop");

  int const convlev = 0;
  cGH *cctkGH = fc->GH[convlev];

  // Main loop
  BeginTimingEvolution(cctkGH);
  static Timers::Timer timer("Evolve");
  timer.start();
  while (not do_terminate(cctkGH)) {

    AdvanceTime(cctkGH);
    {
      int const do_every = maxtimereflevelfact / timereffacts.AT(reflevels - 1);
      if ((cctkGH->cctk_iteration - 1) % do_every == 0) {
        ENTER_GLOBAL_MODE(cctkGH, 0) {
          BEGIN_REFLEVEL_LOOP(cctkGH) { CallRegrid(cctkGH); }
          END_REFLEVEL_LOOP;
        }
        LEAVE_GLOBAL_MODE;
      }
    }
    // Serial Hermite also fills during restriction, where theta is 0 or 1.
    // CSAMR postrestrict keeps the existing non-Hermite path.
    hermite_prolongate = hermite_dense_output;
    CallEvol(cctkGH);
    if (use_csamr)
      hermite_prolongate = false;
    CallRestrict(cctkGH);
    hermite_prolongate = false;
    CallAnalysis(cctkGH);
    print_internal_data();

    // Print timer values
    {
      Timers::Timer timer("PrintTimers");
      timer.start();
      int const do_every = maxtimereflevelfact / timereffacts.AT(reflevels - 1);
      if (output_timers_every > 0 and
          cctkGH->cctk_iteration % output_timers_every == 0 and
          cctkGH->cctk_iteration % do_every == 0) {
        Timers::CactusTimerSet::writeData(cctkGH, timer_file);
      }

      if (output_timer_tree_every > 0 and
          cctkGH->cctk_iteration % output_timer_tree_every == 0 and
          cctkGH->cctk_iteration % do_every == 0) {
        Timers::Timer::outputTree("Evolve");
      }
      timer.stop();
    }

    // Ensure that all levels have consistent times
    {
      Timers::Timer timer("CheckLevelTimes");
      timer.start();
      CCTK_REAL const eps =
          pow(numeric_limits<CCTK_REAL>::epsilon(), CCTK_REAL(0.75));
      assert(fabs(cctkGH->cctk_time - global_time) <= eps * global_time);
      for (int ml = 0; ml < mglevels; ++ml) {
        for (int rl = 0; rl < reflevels; ++rl) {
          int const do_every =
              ipow(mgfact, ml) * (maxtimereflevelfact / timereffacts.AT(rl));
          if (cctkGH->cctk_iteration % do_every == 0) {
            // assert (fabs (leveltimes.AT(ml).AT(rl) - global_time) <=
            //         eps * global_time);
            assert(fabs(tt->get_time(ml, rl, 0) - global_time) <=
                   eps * global_time);
          }
        }
      }
      timer.stop();
    }

  } // end main loop
  timer.stop();

  Waypoint("Done with evolution loop");

  return 0;
}

// Rank 0 keeps the only stdout. An early stop has to be explained there.
static void csamr_explain_stop(cGH const *const cctkGH, int const local_term) {
  DECLARE_CCTK_PARAMETERS;

  int nprocs = 0;
  int rank = 0;
  if (MPI_Comm_size(MPI_COMM_WORLD, &nprocs) != MPI_SUCCESS or
      MPI_Comm_rank(MPI_COMM_WORLD, &rank) != MPI_SUCCESS) {
    CCTK_ERROR("CSAMR could not read the MPI world size at termination.");
  }
  if (nprocs < 1) {
    CCTK_ERROR("CSAMR termination check has an empty MPI world.");
  }

  double dsend[4] = {cctkGH->cctk_time, global_time, delta_time,
                     csamr_nan_counter()};
  int isend[4] = {CCTK_TerminationReached(cctkGH), cctkGH->cctk_iteration,
                  local_term, maxtimereflevelfact};
  vector<double> drecv(rank == 0 ? 4 * nprocs : 0);
  vector<int> irecv(rank == 0 ? 4 * nprocs : 0);
  if (MPI_Gather(dsend, 4, MPI_DOUBLE, rank == 0 ? &drecv.front() : NULL, 4,
                 MPI_DOUBLE, 0, MPI_COMM_WORLD) != MPI_SUCCESS or
      MPI_Gather(isend, 4, MPI_INT, rank == 0 ? &irecv.front() : NULL, 4,
                 MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS) {
    CCTK_ERROR("CSAMR could not gather the per-rank termination state.");
  }

  int fatal = 0;
  if (rank == 0) {
    ostringstream msg;
    msg << "CSAMR termination at iteration " << cctkGH->cctk_iteration
        << " with cctk_final_time " << cctk_final_time << ".";
    for (int proc = 0; proc < nprocs; ++proc) {
      msg << " rank " << proc << ": time " << drecv.AT(4 * proc) << " global "
          << drecv.AT(4 * proc + 1) << " delta " << drecv.AT(4 * proc + 2)
          << " nans " << drecv.AT(4 * proc + 3) << " reached "
          << irecv.AT(4 * proc) << " iter " << irecv.AT(4 * proc + 1)
          << " local " << irecv.AT(4 * proc + 2) << " factor "
          << irecv.AT(4 * proc + 3) << ".";
    }
    CCTK_REAL const tol = 1.0e-8 * fabs(delta_time);
    bool const early = cctkGH->cctk_time < cctk_final_time - tol;
    if (early) {
      fatal = 1;
      CCTK_VInfo(CCTK_THORNSTRING, "%s", msg.str().c_str());
      fflush(stdout);
    }
  }
  if (MPI_Bcast(&fatal, 1, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS) {
    CCTK_ERROR("CSAMR could not broadcast the termination decision.");
  }
  if (fatal) {
    CCTK_ERROR("CSAMR stopped before cctk_final_time. Rank 0 lists each "
               "rank's time, delta_time, NaN count, and termination flag.");
  }
}

static void csamr_require_same_clock() {
  double const mine_dt = delta_time;
  double root_dt = mine_dt;
  int const mine_factor = maxtimereflevelfact;
  int root_factor = mine_factor;
  if (MPI_Bcast(&root_dt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD) != MPI_SUCCESS or
      MPI_Bcast(&root_factor, 1, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS) {
    CCTK_ERROR("CSAMR could not compare the evolution clock across ranks.");
  }
  if (mine_dt != root_dt or mine_factor != root_factor) {
    CCTK_VERROR("CSAMR rank %d has delta_time %.17g and maxtimereflevelfact %d; "
                "rank 0 has delta_time %.17g and maxtimereflevelfact %d.",
                CCTK_MyProc(NULL), mine_dt, mine_factor, root_dt, root_factor);
  }
}

bool do_terminate(const cGH *cctkGH) {
  DECLARE_CCTK_PARAMETERS;

  static Timers::Timer timer("DoTerminate");
  timer.start();

  bool term;

  // Do not test on non-active reflevels to save the call to
  // MPI_Allreduce below
  int const do_every = maxtimereflevelfact / timereffacts.AT(reflevels - 1);
  if (cctkGH->cctk_iteration % do_every != 0) {

    term = false;

  } else {

    if (terminate_next or CCTK_TerminationReached(cctkGH)) {

      // Terminate if someone or something said so
      term = true;

    } else {

      // Test the various conditions
      bool const term_iter =
          cctk_itlast >= 0 and cctkGH->cctk_iteration >= cctk_itlast;
      bool const term_time =
          (delta_time > 0.0
               ? (cctkGH->cctk_time >=
                  cctk_final_time - 1.0e-8 * cctkGH->cctk_delta_time)
               : (cctkGH->cctk_time <=
                  cctk_final_time + 1.0e-8 * cctkGH->cctk_delta_time));
      bool const term_runtime =
          max_runtime > 0.0 and CCTK_RunTime() >= 60.0 * max_runtime;

      if (CCTK_Equals(terminate, "never")) {
        term = false;
      } else if (CCTK_Equals(terminate, "iteration")) {
        term = term_iter;
      } else if (CCTK_Equals(terminate, "time")) {
        term = term_time;
        if (term) {
          CCTK_VInfo(CCTK_THORNSTRING,
                     "Terminating due to cctk_final_time at t = %f",
                     cctkGH->cctk_time);
        }

      } else if (CCTK_Equals(terminate, "runtime")) {
        term = term_runtime;
      } else if (CCTK_Equals(terminate, "any")) {
        term = term_iter or term_time or term_runtime;
      } else if (CCTK_Equals(terminate, "all")) {
        term = term_iter and term_time and term_runtime;
      } else if (CCTK_Equals(terminate, "either")) {
        term = term_iter or term_time;
      } else if (CCTK_Equals(terminate, "both")) {
        term = term_iter and term_time;
      } else if (CCTK_Equals(terminate, "immediately")) {
        term = true;
      } else {
        CCTK_WARN(0, "Unsupported termination condition");
        abort(); // keep the compiler happy
      }
    }

    // Reduce termination condition
    int local, global;
    local = term;
    MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_LOR, dist::comm());
    term = global;
    if (use_csamr and term)
      csamr_explain_stop(cctkGH, local);
  }

  timer.stop();
  return term;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

void AdvanceTime(cGH *const cctkGH) {
  DECLARE_CCTK_PARAMETERS;

  static Timers::Timer timer("AdvanceTime");
  timer.start();

  Checkpoint("AdvanceTime");

  ++cctkGH->cctk_iteration;

  if (not adaptive_stepsize) {
    if (csamr_ab3_rk4) {
      // The finest step's closed-form end. A uniform fine dt would label
      // the unequal steps with the wrong time and fail CheckLevelTimes.
      int const n_levels = reflevels;
      int const nfine = maxtimereflevelfact;
      int const expect = n_levels >= 1 and n_levels <= csamr::max_levels
                             ? 1 << (n_levels - 1)
                             : 0;
      if (n_levels < 1 or nfine != expect or not(delta_time > 0.0)) {
        CCTK_ERROR("AB3-RK4 clock needs a positive coarse dt and "
                   "maxtimereflevelfact = 2^(levels-1).");
      }
      int const done = cctkGH->cctk_iteration - 1;
      if (done < 0)
        CCTK_ERROR("AB3-RK4 clock ran before the first iteration.");
      int const slot = done % nfine;
      int const cycle = done / nfine;
      int const finest = n_levels - 1;
      int const ordinal = csamr::ab3rk4_ordinal(n_levels, slot, finest);
      if (ordinal != slot) {
        CCTK_ERROR("AB3-RK4 finest step ordinal is not the slot index.");
      }
      csamr::Ab3Rk4Step const sized = csamr::ab3rk4_step(
          finest, ordinal, csamr_ab3_fraction, delta_time);
      global_time = csamr::ab3rk4_absolute_end(cctk_initial_time, cycle, sized,
                                               delta_time);
    } else {
      // Avoid accumulation of errors
      global_time = cctk_initial_time +
                    cctkGH->cctk_iteration * delta_time / maxtimereflevelfact;
    }
    cctkGH->cctk_time = global_time;
  } else {
    // Take varying step sizes into account
    cctkGH->cctk_time += cctkGH->cctk_delta_time;
    delta_time = cctkGH->cctk_delta_time;
    global_time = cctkGH->cctk_time;
  }

  if ((cctkGH->cctk_iteration - 1) %
          (maxtimereflevelfact / timereffacts.AT(reflevels - 1)) ==
      0) {
    Waypoint("Evolving iteration %d at t=%g", cctkGH->cctk_iteration,
             (double)cctkGH->cctk_time);
  }

  timer.stop();
}

void CallRegrid(cGH *const cctkGH) {
  DECLARE_CCTK_PARAMETERS;

  char const *const where = "CallRegrid";
  static Timers::Timer timer(where);
  timer.start();

  assert(is_level_mode());

  bool const old_do_global_mode = do_global_mode;
  bool const old_do_early_global_mode = do_early_global_mode;
  bool const old_do_late_global_mode = do_late_global_mode;
  bool const old_do_meta_mode = do_meta_mode;
  bool const old_do_early_meta_mode = do_early_meta_mode;
  bool const old_do_late_meta_mode = do_late_meta_mode;
  do_global_mode = true;
  do_early_global_mode = true;
  do_late_global_mode = true;
  do_meta_mode = true;
  do_early_meta_mode = true;
  do_late_meta_mode = true;

  // Preregrid
  Waypoint("Preregrid at iteration %d time %g%s%s", cctkGH->cctk_iteration,
           (double)cctkGH->cctk_time, (do_global_mode ? " (global)" : ""),
           (do_meta_mode ? " (meta)" : ""));
  ScheduleTraverse(where, "CCTK_PREREGRID", cctkGH);

  // Regrid
  Checkpoint("Regrid");
  int const oldreflevels = reflevels;
  bool const did_regrid = Regrid(cctkGH, false, true);
  bool const did_remove_level = reflevels < oldreflevels;
  assert(not did_remove_level or did_regrid);

  if (did_regrid) {
#ifdef REQUIREMENTS_HH
    Requirements::Regrid(reflevels);
#endif
    bool did_any_recompose = false;
    BEGIN_META_MODE(cctkGH) {

      bool have_done_global_mode = false;
      bool have_done_early_global_mode = false;
      bool have_done_late_global_mode = false;
      bool have_done_anything = false;

      for (int rl = 0; rl < reflevels; ++rl) {

        bool const did_recompose = Recompose(cctkGH, rl, true);
        did_any_recompose = did_any_recompose or did_recompose;
#ifdef REQUIREMENTS_HH
        Requirements::Recompose(cctkGH->cctk_iteration, rl,
                                not did_recompose
                                    ? Requirements::valid::everywhere
                                    : Requirements::valid::interior);
#endif

        // Carpet assumes that a regridding operation always changes
        // "level N and all finer levels" so we should call POSTREGRID
        // on all finer levels
        if (did_any_recompose or (did_remove_level and rl == reflevels - 1)) {
          BEGIN_MGLEVEL_LOOP(cctkGH) {
            ENTER_LEVEL_MODE(cctkGH, rl) {
              BeginTimingLevel(cctkGH);

              // Regridding should invalidate ghost regions
              int num_groups = CCTK_NumGroups();
              for (int gi=0; gi<num_groups; gi++) {
                int gtype = CCTK_GroupTypeI(gi);
                if (gtype == CCTK_GF) {
                  int num_vars = CCTK_NumVarsInGroupI(gi);
                  int tls = CCTK_ActiveTimeLevelsGI(cctkGH, gi);
                  int first_var = CCTK_FirstVarIndexI(gi);
                  for (int vi=0; vi<num_vars; vi++) {
                    int varindex = vi + first_var;
                    for (int tl=0; tl<tls; tl++) {
                      int valid = Driver_GetValidRegion(cctkGH, varindex, tl);
                      int new_valid = valid & ~CCTK_VALID_GHOSTS;
                      if (valid != new_valid) {
                        Driver_SetValidRegion(cctkGH, varindex, tl, new_valid);
                      }
                    }
                  }
                }
              }

              do_early_global_mode = not have_done_early_global_mode;
              do_late_global_mode = reflevel == reflevels - 1;
              do_early_meta_mode =
                  do_early_global_mode and mglevel == mglevels - 1;
              do_late_meta_mode = do_late_global_mode and mglevel == 0;
              do_global_mode = do_late_global_mode;
              do_meta_mode = do_late_meta_mode;
              assert(not(have_done_global_mode and do_global_mode));
              assert(not(have_done_early_global_mode and do_early_global_mode));
              assert(not(have_done_late_global_mode and do_late_global_mode));
              have_done_global_mode |= do_global_mode;
              have_done_early_global_mode |= do_early_global_mode;
              have_done_late_global_mode |= do_late_global_mode;
              have_done_anything = true;

              BEGIN_TIMELEVEL_LOOP(cctkGH) {

                Waypoint("Postregrid at iteration %d time %g timelevel %d%s%s",
                         cctkGH->cctk_iteration, (double)cctkGH->cctk_time,
                         timelevel, (do_global_mode ? " (global)" : ""),
                         (do_meta_mode ? " (meta)" : ""));

                // Postregrid
                ScheduleTraverse(where, "CCTK_POSTREGRID", cctkGH);
              }
              END_TIMELEVEL_LOOP;

              if (output_after_regridding) {
                // Output
                OutputGH(where, cctkGH);
              }

              EndTimingLevel(cctkGH);
            }
            LEAVE_LEVEL_MODE;
          }
          END_MGLEVEL_LOOP;

        } // if did_recompose

      } // for rl

      if (have_done_anything)
        assert(have_done_global_mode);
      if (have_done_anything)
        assert(have_done_early_global_mode);
      if (have_done_anything)
        assert(have_done_late_global_mode);
    }
    END_META_MODE;
#ifdef REQUIREMENTS_HH
    Requirements::RegridFree();
#endif
  } // if did_regrid

  RegridFree(cctkGH, true);

  do_global_mode = old_do_global_mode;
  do_early_global_mode = old_do_early_global_mode;
  do_late_global_mode = old_do_late_global_mode;
  do_meta_mode = old_do_meta_mode;
  do_early_meta_mode = old_do_early_meta_mode;
  do_late_meta_mode = old_do_late_meta_mode;

  timer.stop();
}

namespace {

struct csamr_mode_flags {
  bool early_global;
  bool late_global;
  bool early_meta;
  bool late_meta;
  bool global_mode;
  bool meta_mode;
};

// Same assignment as the sequential loop, keyed by level. ml is the
// multigrid level about to be entered. The caller is still in meta
// mode, where mglevel is unset, so the test uses ml.
csamr_mode_flags csamr_flags_for_level(int const level, int const coarsest,
                                       int const ml) {
  csamr_mode_flags flags;
  flags.early_global = level == coarsest;
  flags.late_global = level == reflevels - 1;
  flags.early_meta = flags.early_global and ml == mglevels - 1;
  flags.late_meta = flags.late_global and ml == 0;
  flags.global_mode = flags.early_global;
  flags.meta_mode = flags.early_meta;
  return flags;
}

// Parent time is the fine-step index already visible. A child that
// finishes ahead of that time cannot prolongate. A child that lands on
// the parent time copies when use_final_copy is set. A child that
// finishes behind the parent time-interpolates when that refill is set.
bool csamr_taper(int const rl, int const end_index, vector<int> const &at) {
  DECLARE_CCTK_PARAMETERS;
  if (rl == 0)
    return false;
  return csamr::taper_step(end_index, at.AT(rl - 1), use_time_interpolation,
                           use_final_copy);
}

struct csamr_level_plan {
  int rl;
  int partner;
  CCTK_REAL end_time;
  CCTK_REAL step_dt;
  bool unequal;
  bool taper;
  int end_index;
  // The step end is not yet stored on the parent, or this slot updates
  // the parent. Prolongation waits until the parent covers that end.
  bool ab3_prefill;
  csamr_mode_flags flags;
};

bool csamr_carpetlib_param(char const *const name) {
  int type = 0;
  void const *const ptr = CCTK_ParameterGet(name, "CarpetLib", &type);
  if (ptr == NULL) {
    CCTK_VERROR("CSAMR could not read CarpetLib::%s.", name);
  }
  return *static_cast<CCTK_INT const *>(ptr) != 0;
}

// Kind of the step this thread is inside. reflevel is threadprivate, so a
// slot that runs AB3 beside RK4 does not share a flag. Returns 0 when the
// mode is off or this level is coarse.
extern "C" CCTK_INT
Carpet_AB3RK4_ThisStepAB3(CCTK_POINTER_TO_CONST const cctkGH_) {
  DECLARE_CCTK_PARAMETERS;
  if (not csamr_ab3_rk4)
    return 0;
  if (cctkGH_ == nullptr)
    CCTK_ERROR("AB3-RK4 step query has no grid hierarchy.");
  cGH const *const cctkGH = static_cast<cGH const *>(cctkGH_);
  if (reflevel <= 0)
    return 0;
  int const n_levels = reflevels;
  int const nfine = maxtimereflevelfact;
  int const expect =
      n_levels >= 1 and n_levels <= csamr::max_levels ? 1 << (n_levels - 1) : 0;
  if (n_levels < 2 or nfine != expect)
    CCTK_ERROR("AB3-RK4 step query needs maxtimereflevelfact = 2^(levels-1).");
  if (cctkGH->cctk_iteration < 1)
    CCTK_ERROR("AB3-RK4 step query ran before the first iteration.");
  // The coarse update runs on a cGH allocated on the first paired slot.
  // That cGH has to carry this iteration. Slot 0's iteration makes the
  // later long refined step look like the short AB3 step.
  return csamr::ab3rk4_step_is_ab3(n_levels, cctkGH->cctk_iteration, reflevel)
             ? 1
             : 0;
}

// Rotate this level's grid functions and clock. Scalars stay put: the
// finest update has to see them as they were, and the coarser update
// cycles them later.
void csamr_prepare_level(cGH *const cctkGH, int const ml,
                         csamr_level_plan const &plan) {
  ENTER_GLOBAL_MODE(cctkGH, ml) {
    ENTER_LEVEL_MODE(cctkGH, plan.rl) {
      do_early_global_mode = plan.flags.early_global;
      do_late_global_mode = plan.flags.late_global;
      do_early_meta_mode = plan.flags.early_meta;
      do_late_meta_mode = plan.flags.late_meta;
      do_global_mode = plan.flags.global_mode;
      do_meta_mode = plan.flags.meta_mode;

      if (plan.unequal) {
        if (not(plan.step_dt > 0.0))
          CCTK_ERROR("AB3-RK4 step dt is not positive.");
        tt->set_delta(mglevel, reflevel, plan.step_dt);
      }
      CycleTimeLevels(cctkGH, false);
      cctkGH->cctk_time = plan.end_time;
      tt->set_time(mglevel, reflevel, timelevel, cctkGH->cctk_time);
    }
    LEAVE_LEVEL_MODE;
  }
  LEAVE_GLOBAL_MODE;
}

void csamr_traverse_level(cGH *const cctkGH, int const ml,
                          csamr_level_plan const &plan, bool const taper,
                          bool const cycle_replicated) {
  char const *const where = "CallEvolCsamr";
  bool const time_this = not csamr_overlap_active();
  ENTER_GLOBAL_MODE(cctkGH, ml) {
    ENTER_LEVEL_MODE(cctkGH, plan.rl) {
      if (time_this)
        BeginTimingLevel(cctkGH);

      do_early_global_mode = plan.flags.early_global;
      do_late_global_mode = plan.flags.late_global;
      do_early_meta_mode = plan.flags.early_meta;
      do_late_meta_mode = plan.flags.late_meta;
      do_global_mode = plan.flags.global_mode;
      do_meta_mode = plan.flags.meta_mode;
      do_taper = taper;

      // Prepare already rotated this level's grid functions.
      if (cycle_replicated and plan.flags.global_mode)
        CycleReplicatedLevels(cctkGH);
      cctkGH->cctk_time = plan.end_time;
      tt->set_time(mglevel, reflevel, timelevel, cctkGH->cctk_time);
      if (plan.unequal) {
        if (timereffacts.AT(reflevel) <= 0 or not(plan.step_dt > 0.0))
          CCTK_ERROR("AB3-RK4 cannot set the level timestep.");
        // CCTK_DELTA_TIME is cctk_delta_time / cctk_timefac. Level mode
        // has already stored the level's time factor in cctk_timefac.
        cctkGH->cctk_delta_time =
            plan.step_dt * static_cast<CCTK_REAL>(timereffacts.AT(reflevel));
      }

      Waypoint("Evolution I CSAMR at iteration %d time %g level %d "
               "partner %d%s%s%s%s",
               cctkGH->cctk_iteration, (double)cctkGH->cctk_time, reflevel,
               plan.partner, (do_global_mode ? " (global)" : ""),
               (do_meta_mode ? " (meta)" : ""),
               (plan.ab3_prefill ? " (AB3 ghosts deferred)"
                                 : (do_taper ? " (tapering)" : "")),
               (plan.unequal ? (reflevel > 0 and
                                        Carpet_AB3RK4_ThisStepAB3(cctkGH)
                                    ? " (AB3)"
                                    : " (RK4)")
                             : ""));

      CalculateChecksums(cctkGH, allbutcurrenttime);
      Poison(cctkGH, currenttimebutnotifonly);

      ScheduleTraverse(where, "CCTK_PRESTEP", cctkGH);
      ScheduleTraverse(where, "CCTK_EVOL", cctkGH);

      PoisonCheck(cctkGH, currenttime);
      if (time_this)
        StepTimingEvolution(cctkGH);

      do_taper = false;
      if (time_this)
        EndTimingLevel(cctkGH);
    }
    LEAVE_LEVEL_MODE;
  }
  LEAVE_GLOBAL_MODE;
}


void csamr_check_parameters(cGH const *const cctkGH) {
  DECLARE_CCTK_PARAMETERS;

  if (csamr_ab3_rk4) {
    if (use_tapered_grids or use_buffer_zones) {
      CCTK_ERROR("Carpet::csamr_ab3_rk4 does not use tapered grids or "
                 "buffer zones. The short step is one right-hand side at "
                 "the parent time.");
    }
    if (use_time_interpolation) {
      CCTK_ERROR("Carpet::csamr_ab3_rk4 requires "
                 "Carpet::use_time_interpolation = no.");
    }
    if (not use_final_copy or not hermite_dense_output) {
      CCTK_ERROR("Carpet::csamr_ab3_rk4 requires Carpet::use_final_copy and "
                 "Carpet::hermite_dense_output.");
    }
    if (not csamr::ab3rk4_fraction_ok(csamr_ab3_fraction)) {
      CCTK_ERROR("Carpet::csamr_ab3_fraction must lie in (0, 0.5).");
    }
    for (int rl = 0; rl < reflevels; ++rl) {
      if (timereffacts.AT(rl) != (1 << rl)) {
        CCTK_ERROR("Carpet::csamr_ab3_rk4 requires time factors 1,2,4,... "
                   "with one refinement level per factor.");
      }
    }
    int type = 0;
    void const *const method = CCTK_ParameterGet("ODE_Method", "MoL", &type);
    if (method == nullptr or type != PARAMETER_KEYWORD) {
      CCTK_ERROR("Carpet::csamr_ab3_rk4 could not read MoL::ODE_Method.");
    }
    char const *const ode = *static_cast<char const *const *>(method);
    if (not CCTK_EQUALS(ode, "AB3-RK4")) {
      CCTK_ERROR("Carpet::csamr_ab3_rk4 requires MoL::ODE_Method = "
                 "\"AB3-RK4\".");
    }
  } else {
    csamr_check_ab3_buffers(true);
  }
  if (adaptive_stepsize)
    CCTK_ERROR("Carpet::use_csamr does not support adaptive timesteps.");
  if (hermite_dense_output and not use_final_copy) {
    CCTK_ERROR("Carpet::hermite_dense_output with Carpet::use_csamr fills "
               "the fine step that lands on the parent time. That requires "
               "Carpet::use_final_copy.");
  }
  if (mglevels != 1)
    CCTK_ERROR("Carpet::use_csamr does not support multigrid.");
  if (cctkGH->cctk_iteration < 1)
    CCTK_ERROR("Carpet::use_csamr was called before the first iteration.");
  if (reflevels < 1) {
    CCTK_ERROR("Carpet::use_csamr was called with no refinement levels.");
  }
  vector<int> factors(reflevels);
  for (int rl = 0; rl < reflevels; ++rl)
    factors[rl] = timereffacts.AT(rl);
  csamr::FactorError const factors_ok = csamr::classify_factors(factors);
  if (factors_ok != csamr::FactorError::Ok) {
    CCTK_VERROR("Carpet::use_csamr requires time factors 1,2,4,... with "
                "equal-timestep levels grouped together. Level count %d is "
                "not that hierarchy.",
                reflevels);
  }
  if (maxtimereflevelfact != timereffacts.AT(reflevels - 1)) {
    CCTK_VERROR("Carpet::use_csamr requires maxtimereflevelfact %d, found %d.",
                timereffacts.AT(reflevels - 1), maxtimereflevelfact);
  }
}

// End times are counted from the start of a coarse cycle, and the
// initial prolongation lives on slot 0. A new level count is safe only
// when both the old cycle and the new one are on that boundary.
void csamr_note_reflevels(cGH const *const cctkGH, int const nfine) {
  static int active_reflevels = -1;
  static int active_nfine = -1;
  if (active_reflevels < 0) {
    active_reflevels = reflevels;
    active_nfine = nfine;
    return;
  }
  if (reflevels == active_reflevels)
    return;

  int const steps_done = cctkGH->cctk_iteration - 1;
  if (steps_done % active_nfine != 0 or steps_done % nfine != 0) {
    CCTK_VERROR("Carpet::use_csamr cannot change the number of refinement "
                "levels from %d to %d except at a shared coarse-cycle "
                "boundary (iteration %d).",
                active_reflevels, reflevels, cctkGH->cctk_iteration);
  }
  active_reflevels = reflevels;
  active_nfine = nfine;
}

void csamr_check_restricts(cGH const *const cctkGH, int const ml,
                           vector<csamr::Group> const &groups,
                           csamr::Slot const &slot) {
  vector<int> hits(groups.size(), 0);
  for (csamr::Action const &action : slot.restricts) {
    if (action.kind != csamr::Kind::Restrict)
      CCTK_ERROR("CSAMR restrict list contains a non-restrict action.");
    if (action.level < 0 or action.level >= static_cast<int>(groups.size()) or
        action.other != action.level + 1) {
      CCTK_VERROR("CSAMR restriction from group %d onto group %d is not an "
                  "adjacent restriction.",
                  action.other, action.level);
    }
    hits.AT(action.level) += 1;
  }

  for (int rl = reflevels - 2; rl >= 0; --rl) {
    int const group = csamr::group_of(groups, rl);
    int const do_every =
        ipow(mgfact, ml) * (maxtimereflevelfact / timereffacts.AT(rl));
    bool const carpet_restricts = cctkGH->cctk_iteration % do_every == 0;
    if (hits.AT(group) > 1) {
      CCTK_VERROR("CSAMR restricts onto group %d more than once in this step.",
                  group);
    }
    if (carpet_restricts != (hits.AT(group) > 0)) {
      CCTK_VERROR("CSAMR restriction onto level %d at iteration %d disagrees "
                  "with Carpet (schedule %d, Carpet %d).",
                  rl, cctkGH->cctk_iteration, hits.AT(group) > 0 ? 1 : 0,
                  carpet_restricts ? 1 : 0);
    }
  }
}

} // namespace

void csamr_fail(bool const fatal, char const *const message) {
  if (fatal)
    CCTK_ERROR(message);
  else
    CCTK_PARAMWARN(message);
}

bool csamr_mol_is(char const *const name, char const *const expected) {
  int type = 0;
  void const *const parameter = CCTK_ParameterGet(name, "MoL", &type);
  if (parameter == nullptr or type != PARAMETER_KEYWORD) {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                "Carpet could not read MoL::%s.", name);
  }
  char const *const value = *static_cast<char const *const *>(parameter);
  return CCTK_EQUALS(value, expected);
}

int csamr_mol_int(char const *const name) {
  int type = 0;
  void const *const parameter = CCTK_ParameterGet(name, "MoL", &type);
  if (parameter == nullptr or type != PARAMETER_INT) {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                "Carpet could not read MoL::%s.", name);
  }
  return static_cast<int>(*static_cast<CCTK_INT const *>(parameter));
}

bool csamr_mol_bool(char const *const name) {
  int type = 0;
  void const *const parameter = CCTK_ParameterGet(name, "MoL", &type);
  if (parameter == nullptr or type != PARAMETER_BOOLEAN) {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                "Carpet could not read MoL::%s.", name);
  }
  return *static_cast<CCTK_INT const *>(parameter) != 0;
}

bool csamr_ab3_ode() {
  DECLARE_CCTK_PARAMETERS;
  if (csamr_ab3_rk4 or not CCTK_IsThornActive("MoL"))
    return false;
  return csamr_mol_is("ODE_Method", "AB") and csamr_mol_is("AB_Type", "3");
}

bool csamr_stock_ab3() {
  if (not csamr_ab3_ode())
    return false;
  return csamr_mol_int("MoL_Intermediate_Steps") == 1;
}

bool csamr_ab3_unbuffered() {
  DECLARE_CCTK_PARAMETERS;
  return use_csamr and csamr_stock_ab3() and not use_buffer_zones and
         not use_tapered_grids;
}

void csamr_check_ab3_buffers(bool const fatal) {
  DECLARE_CCTK_PARAMETERS;
  if (not use_csamr or csamr_ab3_rk4)
    return;
  if (csamr_stock_ab3()) {
    if (use_buffer_zones or use_tapered_grids) {
      csamr_fail(fatal,
                 "Stock Adams-Bashforth 3 does not use buffer zones or "
                 "tapered grids. Its one right-hand side is evaluated at "
                 "the step start, and the ghost zones cover that "
                 "evaluation.");
    }
    if (not use_time_interpolation or not use_final_copy) {
      csamr_fail(fatal,
                 "Stock Adams-Bashforth 3 refills ghosts by time "
                 "interpolation between parent times and by a copy on the "
                 "parent time.");
    }
    if (hermite_dense_output) {
      csamr_fail(fatal, "Stock Adams-Bashforth 3 does not use Hermite "
                        "dense output.");
    }
    if (not csamr_mol_bool("AB3_rhs_storage")) {
      csamr_fail(fatal, "Stock Adams-Bashforth 3 needs "
                        "MoL::AB3_rhs_storage = yes.");
    }
    if (csamr_mol_bool("skip_initial_copy")) {
      csamr_fail(fatal,
                 "Stock Adams-Bashforth 3 CSAMR needs "
                 "MoL::skip_initial_copy = no. The ghost fill is rotated "
                 "with the time levels and copied back onto the "
                 "right-hand side.");
    }
    return;
  }
  if (csamr_ab3_ode() and (not use_buffer_zones or not use_tapered_grids)) {
    csamr_fail(fatal, "Adams-Bashforth 3 without buffer zones needs "
                      "MoL::MoL_Intermediate_Steps = 1.");
  }
  if (not use_tapered_grids) {
    csamr_fail(fatal,
               "Carpet::use_csamr requires Carpet::use_tapered_grids. "
               "Tapered grids and buffer zones stop an RK stage from "
               "prolongating out of the parent during the step.");
  }
  if (not use_buffer_zones) {
    csamr_fail(fatal,
               "Carpet::use_csamr requires Carpet::use_buffer_zones. "
               "Without buffer zones a refined grid has only one "
               "ghost-width, which does not cover a multi-stage "
               "Runge-Kutta step.");
  }
}


// Fine-step index whose ghosts still have to be written on this level.
// -1 means the current time level is already filled.
int &csamr_ab3_pending(int const rl) {
  static int pending[csamr::max_levels];
  static bool init = false;
  if (not init) {
    for (int i = 0; i < csamr::max_levels; ++i)
      pending[i] = -1;
    init = true;
  }
  if (rl < 0 or rl >= csamr::max_levels)
    CCTK_ERROR("Stock AB3 ghost level is out of range.");
  return pending[rl];
}

// Prolongate this level from its parent at the time currently stored on
// the level. The parent already covers that time. MoL's initial copy
// later puts these ghosts on the next right-hand side.
void csamr_ab3_prefill_ghosts(cGH *const cctkGH, int const ml,
                              csamr_level_plan const &plan) {
  if (plan.rl <= 0)
    CCTK_ERROR("Stock AB3 tried to pre-fill the coarse level.");
  if (&EvolvedVariableIndex == nullptr or EvolvedVariableIndex == nullptr or
      &MoLNumEvolvedVariables == nullptr or
      &MoLNumEvolvedArrayVariables == nullptr) {
    CCTK_ERROR("Stock AB3 ghost fill needs MoL's evolved-variable table.");
  }
  if (MoLNumEvolvedArrayVariables > 0)
    CCTK_ERROR("Stock AB3 CSAMR does not support evolved grid arrays.");
  if (MoLNumEvolvedVariables <= 0)
    CCTK_ERROR("Stock AB3 ghost fill found no evolved variables.");

  vector<int> groups;
  for (int i = 0; i < MoLNumEvolvedVariables; ++i) {
    int const vi = EvolvedVariableIndex[i];
    int const group = CCTK_GroupIndexFromVarI(vi);
    if (group < 0) {
      CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                  "Stock AB3 cannot find the group of evolved variable %d.",
                  vi);
    }
    bool seen = false;
    for (size_t g = 0; g < groups.size(); ++g) {
      if (groups.AT(g) == group)
        seen = true;
    }
    if (not seen)
      groups.push_back(group);
  }
  if (groups.empty())
    CCTK_ERROR("Stock AB3 ghost fill found no evolved groups.");

  ENTER_GLOBAL_MODE(cctkGH, ml) {
    ENTER_LEVEL_MODE(cctkGH, plan.rl) {
      do_taper = false;
      do_prolongate = true;
      cctkGH->cctk_time = tt->get_time(mglevel, reflevel, 0);
      int const status =
          SyncGroupsByDirI(cctkGH, static_cast<int>(groups.size()),
                           &groups.front(), NULL);
      if (status < 0) {
        CCTK_ERROR("Stock AB3 could not prolongate ghosts before the "
                   "parent step.");
      }
      do_taper = false;
    }
    LEAVE_LEVEL_MODE;
  }
  LEAVE_GLOBAL_MODE;
}

// A previous step ended before its parent had stored that time. Fill
// those ghosts once the parent has caught up, before this level steps.
void csamr_ab3_catch_up_ghosts(cGH *const cctkGH, int const ml,
                               vector<int> const &at_before) {
  DECLARE_CCTK_PARAMETERS;
  if (cctkGH->cctk_iteration <= 1) {
    for (int rl = 0; rl < reflevels and rl < csamr::max_levels; ++rl)
      csamr_ab3_pending(rl) = -1;
  }
  CCTK_REAL const fine_dt =
      delta_time / static_cast<CCTK_REAL>(maxtimereflevelfact);
  for (int rl = 1; rl < reflevels; ++rl) {
    int const pending = csamr_ab3_pending(rl);
    if (pending < 0)
      continue;
    if (rl - 1 >= static_cast<int>(at_before.size()) or
        at_before.AT(rl - 1) < pending)
      continue;
    CCTK_REAL const expected =
        cctk_initial_time + fine_dt * static_cast<CCTK_REAL>(pending);
    CCTK_REAL const actual = tt->get_time(ml, rl, 0);
    CCTK_REAL const scale = 1.0 + std::fabs(expected);
    if (not(std::fabs(actual - expected) <= 1.0e-12 * scale))
      CCTK_ERROR("Stock AB3 deferred ghost time does not match the "
                 "level clock.");
    csamr_level_plan plan;
    plan.rl = rl;
    csamr_ab3_prefill_ghosts(cctkGH, ml, plan);
    csamr_ab3_pending(rl) = -1;
  }
}

// Prolongate the steps that skipped their end sync. A parent that
// stepped in this slot now holds its new time. A parent that is still
// behind stays pending until a later slot.
void csamr_ab3_close_slot_ghosts(cGH *const cctkGH, int const ml,
                                 vector<csamr_level_plan> const &plans,
                                 vector<int> const &at_after) {
  DECLARE_CCTK_PARAMETERS;
  CCTK_REAL const fine_dt =
      delta_time / static_cast<CCTK_REAL>(maxtimereflevelfact);
  for (csamr_level_plan const &plan : plans) {
    if (not plan.ab3_prefill)
      continue;
    if (plan.rl <= 0)
      CCTK_ERROR("Stock AB3 deferred the coarse level.");
    if (plan.rl - 1 >= static_cast<int>(at_after.size()))
      CCTK_ERROR("Stock AB3 deferred a level with no parent.");
    int const parent_after = at_after.AT(plan.rl - 1);
    if (parent_after >= plan.end_index) {
      CCTK_REAL const expected =
          cctk_initial_time + fine_dt * static_cast<CCTK_REAL>(plan.end_index);
      CCTK_REAL const actual = tt->get_time(ml, plan.rl, 0);
      CCTK_REAL const scale = 1.0 + std::fabs(expected);
      if (not(std::fabs(actual - expected) <= 1.0e-12 * scale))
        CCTK_ERROR("Stock AB3 end ghost time does not match the level "
                   "clock.");
      if (not(std::fabs(plan.end_time - expected) <= 1.0e-12 * scale))
        CCTK_ERROR("Stock AB3 end ghost time does not match the step.");
      csamr_ab3_prefill_ghosts(cctkGH, ml, plan);
      csamr_ab3_pending(plan.rl) = -1;
    } else {
      csamr_ab3_pending(plan.rl) = plan.end_index;
    }
  }
}

void CallEvolCsamr(cGH *const cctkGH) {
  DECLARE_CCTK_PARAMETERS;

  char const *const where = "CallEvolCsamr";
  static Timers::Timer timer(where);
  timer.start();

  csamr_check_parameters(cctkGH);
  vector<int> factors(reflevels);
  for (int rl = 0; rl < reflevels; ++rl)
    factors[rl] = timereffacts.AT(rl);
  vector<csamr::Group> const groups = csamr::subcycling_groups(factors);
  int const nfine = maxtimereflevelfact;
  csamr_note_reflevels(cctkGH, nfine);

  vector<csamr::Slot> const slots =
      csamr::coarse_cycle(static_cast<int>(groups.size()));
  if (static_cast<int>(slots.size()) != nfine)
    CCTK_ERROR("CSAMR coarse cycle does not contain one slot per finest step.");

  int const steps_done = cctkGH->cctk_iteration - 1;
  int const slot_index = steps_done % nfine;
  int const cycle = steps_done / nfine;
  csamr::Slot const &slot = slots.AT(slot_index);
  int const ml = 0;

  csamr_check_restricts(cctkGH, ml, groups, slot);

  if (slot.steps.empty() or slot.steps.size() > 2)
    CCTK_ERROR("CSAMR slot must take one step or one parallel pair.");
  if (slot.steps.AT(0).level != static_cast<int>(groups.size()) - 1)
    CCTK_ERROR("CSAMR slot does not start with the finest level.");

  int coarsest_rl = reflevels;
  for (csamr::Action const &step : slot.steps) {
    if (step.level < 0 or step.level >= static_cast<int>(groups.size()))
      CCTK_ERROR("CSAMR step names an unknown subcycling group.");
    if (groups.AT(step.level).first < coarsest_rl)
      coarsest_rl = groups.AT(step.level).first;
  }

  vector<int> at =
      csamr::fine_index_before_slot(groups, slots, slot_index, cycle);
  vector<int> const at_before = at;
  CCTK_REAL const fine_dt =
      delta_time / static_cast<CCTK_REAL>(maxtimereflevelfact);
  vector<csamr_level_plan> plans;
  bool have_done_global_mode = false;
  bool const ab3_unbuffered = csamr_ab3_unbuffered();
  if (ab3_unbuffered) {
    for (csamr::Group const &group : groups) {
      if (group.first != group.last) {
        CCTK_ERROR("Stock Adams-Bashforth 3 CSAMR requires one level per "
                   "time factor.");
      }
    }
  }
  // The finest group is listed first. Equal-timestep levels inside a
  // group are then listed coarse to fine, so a final copy sees its parent.
  for (csamr::Action const &step : slot.steps) {
    if (step.kind != csamr::Kind::Step)
      CCTK_ERROR("CSAMR slot step is not a step.");
    csamr::Group const &group = groups.AT(step.level);
    int const fine_index = cycle * nfine + step.end_time;
    CCTK_REAL const uniform_end =
        cctk_initial_time + fine_dt * static_cast<CCTK_REAL>(fine_index);
    int const partner = step.other < 0 ? -1 : groups.AT(step.other).first;
    csamr::Ab3Rk4Step sized;
    sized.ab3 = false;
    sized.dt = 0.0;
    sized.start_offset = 0.0;
    sized.end_offset = 0.0;
    if (csamr_ab3_rk4) {
      if (group.first != group.last) {
        CCTK_ERROR("Carpet::csamr_ab3_rk4 requires one level per "
                   "subcycling group.");
      }
      int const ordinal = csamr::ab3rk4_ordinal(reflevels, slot_index, step.level);
      sized = csamr::ab3rk4_step(group.first, ordinal, csamr_ab3_fraction,
                                 delta_time);
    }
    for (int rl = group.first; rl <= group.last; ++rl) {
      csamr_level_plan plan;
      plan.rl = rl;
      plan.partner = partner;
      if (csamr_ab3_rk4) {
        plan.end_time = csamr::ab3rk4_absolute_end(cctk_initial_time, cycle,
                                                   sized, delta_time);
        // MoL subtracts this from the step end. The offset span does not
        // always invert to the absolute start, and that start is the
        // parent timestamp the AB3 prolongate has to match.
        CCTK_REAL const start_time = csamr::ab3rk4_absolute_start(
            cctk_initial_time, cycle, sized, delta_time);
        plan.step_dt = plan.end_time - start_time;
        plan.unequal = true;
        plan.taper = false;
        plan.end_index = fine_index;
        plan.ab3_prefill = false;
      } else if (ab3_unbuffered) {
        plan.end_time = uniform_end;
        plan.step_dt = 0.0;
        plan.unequal = false;
        plan.end_index = fine_index;
        if (rl == 0) {
          plan.ab3_prefill = false;
          plan.taper = false;
        } else {
          bool parent_steps = false;
          for (csamr::Action const &other : slot.steps) {
            csamr::Group const &other_group = groups.AT(other.level);
            if (rl - 1 >= other_group.first and rl - 1 <= other_group.last)
              parent_steps = true;
          }
          bool const during = csamr::ab3_prolongate_during_step(
              fine_index, at_before.AT(rl - 1), parent_steps);
          plan.ab3_prefill = not during;
          plan.taper = plan.ab3_prefill;
        }
      } else {
        plan.end_time = uniform_end;
        plan.step_dt = 0.0;
        plan.unequal = false;
        plan.end_index = fine_index;
        plan.taper = csamr_taper(rl, fine_index, at);
        plan.ab3_prefill = false;
      }
      plan.flags = csamr_flags_for_level(rl, coarsest_rl, ml);
      assert(not(have_done_global_mode and plan.flags.global_mode));
      have_done_global_mode = have_done_global_mode or plan.flags.global_mode;
      plans.push_back(plan);
      at.AT(rl) = fine_index;
    }
  }
  if (not have_done_global_mode)
    CCTK_ERROR("CSAMR slot did not assign a global-mode level.");

  // Every rank owns every level. A pair overlaps only when this
  // process has two OpenMP threads to give out.
  csamr::ThreadTeams const teams = csamr::thread_teams(dist::num_threads());
  bool const paired = slot.steps.size() == 2 and reflevels >= 2;
  int fine_n = 0;
  if (paired) {
    csamr::Group const &fine_group = groups.AT(slot.steps.AT(0).level);
    fine_n = fine_group.last - fine_group.first + 1;
    if (fine_n < 1 or fine_n >= static_cast<int>(plans.size()))
      CCTK_ERROR("CSAMR finest update is not a proper prefix of this slot.");
    for (int i = 0; i < fine_n; ++i) {
      if (plans.AT(i).rl != fine_group.first + i)
        CCTK_ERROR("CSAMR finest update is not a proper prefix of this slot.");
    }
    // A rank that forks while another runs the updates one after the
    // other would deadlock inside MPI_Comm_dup.
    int local_parallel = teams.parallel ? 1 : 0;
    int parallel_ranks = 0;
    if (MPI_Allreduce(&local_parallel, &parallel_ranks, 1, MPI_INT, MPI_SUM,
                      MPI_COMM_WORLD) != MPI_SUCCESS) {
      CCTK_ERROR("CSAMR could not compare OpenMP thread counts.");
    }
    if (parallel_ranks != 0 and parallel_ranks != CCTK_nProcs(cctkGH))
      CCTK_ERROR("CSAMR OpenMP thread counts differ across ranks.");
  }

  if (reflevels >= 2) {
    static bool announced = false;
    if (not announced and CCTK_MyProc(cctkGH) == 0) {
      if (teams.parallel) {
        CCTK_VInfo(CCTK_THORNSTRING,
                   "CSAMR gives the finest update %d OpenMP threads and the "
                   "coarser update %d.",
                   teams.fine, teams.coarse);
      } else {
        CCTK_VInfo(CCTK_THORNSTRING,
                   "CSAMR runs a slot's two updates one after the other on "
                   "%d OpenMP thread.",
                   teams.fine);
      }
      announced = true;
    }
    static bool clock_checked = false;
    if (not clock_checked) {
      csamr_require_same_clock();
      clock_checked = true;
    }
    csamr_reset_nan_counter();
  }

  // One level, one thread, or the solo finest slot: prepare and step
  // each level in plan order. That is the schedule that matched.
  if (ab3_unbuffered) {
    csamr_ab3_catch_up_ghosts(cctkGH, ml, at_before);
    for (csamr_level_plan const &plan : plans) {
      if (csamr_ab3_pending(plan.rl) >= 0)
        CCTK_ERROR("Stock AB3 is stepping a level whose ghosts are still "
                   "deferred.");
    }
  }
  if (not paired or not teams.parallel) {
    for (csamr_level_plan const &plan : plans) {
      csamr_prepare_level(cctkGH, ml, plan);
      csamr_traverse_level(cctkGH, ml, plan, plan.taper, true);
    }
    if (ab3_unbuffered)
      csamr_ab3_close_slot_ghosts(cctkGH, ml, plans, at);
    timer.stop();
    return;
  }

  if (barriers or schedule_barriers or sync_barriers) {
    CCTK_ERROR("Carpet::use_csamr cannot overlap while barriers, "
               "schedule_barriers, or sync_barriers is yes. Those barriers "
               "join every rank.");
  }
  if (csamr_carpetlib_param("barrier_between_stages") or
      csamr_carpetlib_param("check_communication_schedule")) {
    CCTK_ERROR("Carpet::use_csamr cannot overlap while CarpetLib inserts a "
               "world barrier or an alltoall into ghost exchange.");
  }
  if (CCTK_IsFunctionAliased("Accelerator_Cycle")) {
    CCTK_ERROR("Carpet::use_csamr cannot overlap while Accelerator_Cycle is "
               "registered.");
  }
  csamr_reject_evolved_arrays();
  csamr_require_thread_multiple();
  csamr_create_thread_comms();

  for (csamr_level_plan const &plan : plans)
    csamr_prepare_level(cctkGH, ml, plan);

  // The finest update reads replicated variables from before they are
  // cycled. The coarser update, which owns global mode, reads the
  // cycled storage.
  csamr_snapshot_replicated();
  int global_rl = -1;
  for (int i = 0; i < static_cast<int>(plans.size()); ++i) {
    if (not plans.AT(i).flags.global_mode)
      continue;
    if (i < fine_n)
      CCTK_ERROR("CSAMR global mode is on the finest update.");
    if (global_rl >= 0)
      CCTK_ERROR("CSAMR slot assigned global mode twice.");
    global_rl = plans.AT(i).rl;
  }
  if (global_rl < 0)
    CCTK_ERROR("CSAMR pair has no global-mode level.");
  ENTER_GLOBAL_MODE(cctkGH, ml) {
    ENTER_LEVEL_MODE(cctkGH, global_rl) { CycleReplicatedLevels(cctkGH); }
    LEAVE_LEVEL_MODE;
  }
  LEAVE_GLOBAL_MODE;

  cGH *const coarse_gh = csamr_clone_gh(cctkGH);
  csamr_share_groupdata();
  Timers::Timer::suppress_tree(true);
  CCTK_TimersSuspend(1);
  csamr_suppress_lib_timers(true);

#ifdef _OPENMP
  int const saved_nested = omp_get_nested();
  int const saved_levels = omp_get_max_active_levels();
  int const saved_threads = omp_get_max_threads();
  int const saved_dynamic = omp_get_dynamic();
  omp_set_nested(1);
  if (saved_levels < 2)
    omp_set_max_active_levels(2);
  omp_set_dynamic(0);

  // Inner RHS teams read the cGH passed into the traverse. They do
  // not copy these threadprivate variables.
#pragma omp parallel num_threads(2) copyin(CARPET_MODE_THREADPRIVATE_SCALARS)
  {
    if (omp_get_num_threads() != 2)
      CCTK_ERROR("CSAMR overlap did not get two OpenMP threads.");
    int const tid = omp_get_thread_num();
    if (tid == 0) {
      dist::set_thread_comm(csamr_fine_comm());
      omp_set_num_threads(teams.fine);
      csamr_use_replicated_snapshot(true);
      csamr_set_overlap_worker(false);
      csamr_set_overlap_active(true);
    } else if (tid == 1) {
      dist::set_thread_comm(csamr_coarse_comm());
      omp_set_num_threads(teams.coarse);
      csamr_use_replicated_snapshot(false);
      csamr_set_overlap_worker(true);
      csamr_set_overlap_active(true);
      csamr_privatize_group_info();
    } else {
      CCTK_ERROR("CSAMR overlap launched more than two threads.");
    }
#pragma omp barrier
    if (tid == 0) {
      for (int i = 0; i < fine_n; ++i) {
        csamr_traverse_level(cctkGH, ml, plans.AT(i), plans.AT(i).taper,
                             false);
      }
    } else {
      for (int i = fine_n; i < static_cast<int>(plans.size()); ++i) {
        csamr_traverse_level(coarse_gh, ml, plans.AT(i), plans.AT(i).taper,
                             false);
      }
      csamr_release_group_info();
    }
    dist::set_thread_comm(MPI_COMM_NULL);
    csamr_use_replicated_snapshot(false);
    csamr_set_overlap_worker(false);
    csamr_set_overlap_active(false);
    omp_set_num_threads(saved_threads);
  }

  omp_set_nested(saved_nested);
  omp_set_max_active_levels(saved_levels);
  omp_set_num_threads(saved_threads);
  omp_set_dynamic(saved_dynamic);
#else
  CCTK_ERROR("CSAMR thread overlap needs OpenMP.");
#endif

  Timers::Timer::suppress_tree(false);
  CCTK_TimersSuspend(0);
  csamr_suppress_lib_timers(false);
  csamr_free_gh(coarse_gh);
  csamr_release_replicated();
  csamr_release_shared_groupdata();
  if (ab3_unbuffered)
    csamr_ab3_close_slot_ghosts(cctkGH, ml, plans, at);
  timer.stop();
}

void CallEvol(cGH *const cctkGH) {
  DECLARE_CCTK_PARAMETERS;

  if (use_csamr) {
    CallEvolCsamr(cctkGH);
    return;
  }

  char const *const where = "CallEvol";
  static Timers::Timer timer(where);
  timer.start();

  for (int ml = mglevels - 1; ml >= 0; --ml) {

    bool have_done_global_mode = false;
    bool have_done_anything = false;

    for (int rl = 0; rl < reflevels; ++rl) {
      int const do_every =
          ipow(mgfact, ml) * (maxtimereflevelfact / timereffacts.AT(rl));
      if ((cctkGH->cctk_iteration - 1) % do_every == 0) {
        ENTER_GLOBAL_MODE(cctkGH, ml) {
          ENTER_LEVEL_MODE(cctkGH, rl) {
            BeginTimingLevel(cctkGH);

            do_early_global_mode = not have_done_global_mode;
            do_late_global_mode = reflevel == reflevels - 1;
            do_early_meta_mode =
                do_early_global_mode and mglevel == mglevels - 1;
            do_late_meta_mode = do_late_global_mode and mglevel == 0;
            do_global_mode = do_early_global_mode;
            do_meta_mode = do_early_meta_mode;
            assert(not(have_done_global_mode and do_global_mode));
            have_done_global_mode |= do_global_mode;
            have_done_anything = true;

            SetTaperForSubcycling(cctkGH);

            // Advance times
            CycleTimeLevels(cctkGH);
            if (not adaptive_stepsize) {
              cctkGH->cctk_time =
                  (global_time - delta_time / maxtimereflevelfact +
                   delta_time * mglevelfact / timereflevelfact);
            }
            tt->set_time(mglevel, reflevel, timelevel, cctkGH->cctk_time);

            Waypoint("Evolution I at iteration %d time %g%s%s%s",
                     cctkGH->cctk_iteration, (double)cctkGH->cctk_time,
                     (do_global_mode ? " (global)" : ""),
                     (do_meta_mode ? " (meta)" : ""),
                     (do_taper ? " (tapering)" : ""));

            // Checking
            CalculateChecksums(cctkGH, allbutcurrenttime);
            Poison(cctkGH, currenttimebutnotifonly);

            // Evolve
            ScheduleTraverse(where, "CCTK_PRESTEP", cctkGH);
            ScheduleTraverse(where, "CCTK_EVOL", cctkGH);

            // Checking
            PoisonCheck(cctkGH, currenttime);

            // Timing statistics
            StepTimingEvolution(cctkGH);

            do_taper = false;

            EndTimingLevel(cctkGH);
          }
          LEAVE_LEVEL_MODE;
        }
        LEAVE_GLOBAL_MODE;
      } // if do_every
    }   // for rl

    if (have_done_anything)
      assert(have_done_global_mode);

  } // for ml

  timer.stop();
}

void CallRestrict(cGH *const cctkGH) {
  DECLARE_CCTK_PARAMETERS;

  char const *const where = "Evolve::CallRestrict";
  static Timers::Timer timer("CallRestrict");
  timer.start();

  for (int ml = mglevels - 1; ml >= 0; --ml) {

    bool did_restrict = false;

    for (int rl = reflevels - 2; rl >= 0; --rl) {
      int const do_every =
          ipow(mgfact, ml) * (maxtimereflevelfact / timereffacts.AT(rl));
      if (cctkGH->cctk_iteration % do_every == 0) {
        ENTER_GLOBAL_MODE(cctkGH, ml) {
          ENTER_LEVEL_MODE(cctkGH, rl) {
            BeginTimingLevel(cctkGH);

            Waypoint("Evolution/Restrict at iteration %d time %g",
                     cctkGH->cctk_iteration, (double)cctkGH->cctk_time);

            Restrict(cctkGH);
            did_restrict = true;

            if (use_higher_order_restriction) {
              do_early_global_mode = false;
              do_late_global_mode = false;
              do_early_meta_mode = false;
              do_late_meta_mode = false;
              do_global_mode = false;
              do_meta_mode = false;

              SetTaperForSubcycling(cctkGH);

              // TODO: disable prolongation (but not
              // synchronization)
              Waypoint("Evolution/PostRestrict (intermediate) at iteration %d "
                       "time %g",
                       cctkGH->cctk_iteration, (double)cctkGH->cctk_time);

              ScheduleTraverse(where, "CCTK_POSTRESTRICT", cctkGH);

              do_taper = false;
            }

            EndTimingLevel(cctkGH);
          }
          LEAVE_LEVEL_MODE;
        }
        LEAVE_GLOBAL_MODE;
      } // if do_every
    }   // for rl

    if (did_restrict) {

      bool have_done_global_mode = false;
      bool have_done_early_global_mode = false;
      bool have_done_late_global_mode = false;
      bool have_done_anything = false;

      for (int rl = 0; rl < reflevels; ++rl) {
        int const do_every =
            ipow(mgfact, ml) * (maxtimereflevelfact / timereffacts.AT(rl));
        if (cctkGH->cctk_iteration % do_every == 0) {
          ENTER_GLOBAL_MODE(cctkGH, ml) {
            ENTER_LEVEL_MODE(cctkGH, rl) {
              BeginTimingLevel(cctkGH);

              do_early_global_mode = not have_done_early_global_mode;
              do_late_global_mode = reflevel == reflevels - 1;
              do_early_meta_mode =
                  do_early_global_mode and mglevel == mglevels - 1;
              do_late_meta_mode = do_late_global_mode and mglevel == 0;
              do_global_mode = do_late_global_mode;
              do_meta_mode = do_global_mode and do_late_meta_mode;
              assert(not(have_done_global_mode and do_global_mode));
              assert(not(have_done_early_global_mode and do_early_global_mode));
              assert(not(have_done_late_global_mode and do_late_global_mode));
              have_done_global_mode |= do_global_mode;
              have_done_early_global_mode |= do_early_global_mode;
              have_done_late_global_mode |= do_late_global_mode;
              have_done_anything = true;

              SetTaperForSubcycling(cctkGH);

              Waypoint("Evolution/PostRestrict at iteration %d time %g",
                       cctkGH->cctk_iteration, (double)cctkGH->cctk_time);

              ScheduleTraverse(where, "CCTK_POSTRESTRICT", cctkGH);

              do_taper = false;

              EndTimingLevel(cctkGH);
            }
            LEAVE_LEVEL_MODE;
          }
          LEAVE_GLOBAL_MODE;
        } // if do_every
      }   // for rl

      if (have_done_anything)
        assert(have_done_global_mode);
      if (have_done_anything)
        assert(have_done_early_global_mode);
      if (have_done_anything)
        assert(have_done_late_global_mode);

    } // if did_restrict

  } // for ml

  timer.stop();
}

void CallAnalysis(cGH *const cctkGH) {
  DECLARE_CCTK_PARAMETERS;

  char const *const where = "CallAnalysis";
  static Timers::Timer timer(where);
  timer.start();

  for (int ml = mglevels - 1; ml >= 0; --ml) {

    bool have_done_global_mode = false;
    bool have_done_early_global_mode = false;
    bool have_done_late_global_mode = false;
    bool have_done_anything = false;

    for (int rl = 0; rl < reflevels; ++rl) {
      int const do_every =
          ipow(mgfact, ml) * (maxtimereflevelfact / timereffacts.AT(rl));
      if (cctkGH->cctk_iteration % do_every == 0) {
        ENTER_GLOBAL_MODE(cctkGH, ml) {
          ENTER_LEVEL_MODE(cctkGH, rl) {
            BeginTimingLevel(cctkGH);

            do_early_global_mode = not have_done_early_global_mode;
            do_late_global_mode = reflevel == reflevels - 1;
            do_early_meta_mode =
                do_early_global_mode and mglevel == mglevels - 1;
            do_late_meta_mode = do_late_global_mode and mglevel == 0;
            do_global_mode = do_late_global_mode;
            do_meta_mode = do_global_mode and do_late_meta_mode;
            assert(not(have_done_global_mode and do_global_mode));
            assert(not(have_done_early_global_mode and do_early_global_mode));
            assert(not(have_done_late_global_mode and do_late_global_mode));
            have_done_global_mode |= do_global_mode;
            have_done_early_global_mode |= do_early_global_mode;
            have_done_late_global_mode |= do_late_global_mode;
            have_done_anything = true;

            SetTaperForSubcycling(cctkGH);

            Waypoint("Evolution II at iteration %d time %g%s%s%s",
                     cctkGH->cctk_iteration, (double)cctkGH->cctk_time,
                     (do_global_mode ? " (global)" : ""),
                     (do_meta_mode ? " (meta)" : ""),
                     (do_taper ? " (tapering)" : ""));

#if 0
              if (reflevel < reflevels-1) {
                ScheduleTraverse (where, "CCTK_POSTRESTRICT", cctkGH);
              }
#endif

            // Poststep
            ScheduleTraverse(where, "CCTK_POSTSTEP", cctkGH);

            // Checking
            PoisonCheck(cctkGH, currenttime);
            CalculateChecksums(cctkGH, currenttime);

            // Checkpoint
            ScheduleTraverse(where, "CCTK_CHECKPOINT", cctkGH);

            // Analysis
            in_analysis_bin = true;
            ScheduleTraverse(where, "CCTK_ANALYSIS", cctkGH);
            in_analysis_bin = false;

            if (do_late_global_mode) {
              // Timing statistics
              UpdateTimingStats(cctkGH);
            }

            // Output
            OutputGH(where, cctkGH);

            // Checking
            CheckChecksums(cctkGH, alltimes);

            if (do_late_global_mode) {
              // Timing statistics
              PrintTimingStats(cctkGH);
            }

            do_taper = false;

            EndTimingLevel(cctkGH);
          }
          LEAVE_LEVEL_MODE;
        }
        LEAVE_GLOBAL_MODE;
      } // if do_every
    }   // for rl

    if (have_done_anything)
      assert(have_done_global_mode);
    if (have_done_anything)
      assert(have_done_early_global_mode);
    if (have_done_anything)
      assert(have_done_late_global_mode);

  } // for ml

  timer.stop();
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

void print_internal_data() {
  DECLARE_CCTK_PARAMETERS;

  if (output_internal_data) {
    CCTK_INFO("Internal data dump:");
    streamsize const oldprecision = cout.precision();
    cout.precision(17);
    cout << "   global_time: " << global_time << endl
         // << "   leveltimes: " << leveltimes << endl
         << "   delta_time: " << delta_time << endl;
    cout.precision(oldprecision);
  }
}

void ScheduleTraverse(char const *const where, char const *const name,
                      cGH *const cctkGH) {
  Timers::Timer timer(name);
  timer.start();
  ostringstream infobuf;
  infobuf << "Scheduling " << name;
  string const info = infobuf.str();
  Checkpoint(info.c_str());
  CCTK_ScheduleTraverse(name, cctkGH, CallFunction);
  timer.stop();
}

void OutputGH(char const *const where, cGH *const cctkGH) {
  static Timers::Timer timer("OutputGH");
  timer.start();
  CCTK_OutputGH(cctkGH);
  timer.stop();
}

} // namespace Carpet
