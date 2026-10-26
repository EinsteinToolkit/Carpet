#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <vector>

#include <cctk.h>
#include <cctk_Parameters.h>

#include <Requirements.hh>

#include <Timer.hh>

#include <ggf.hh>
#include <gh.hh>

#include <carpet.hh>
#include "modes.hh"
#include "PreSync.hh"
#include "hermite_fill.hh"

extern "C" CCTK_INT Carpet_AB3RK4_ThisStepAB3(CCTK_POINTER_TO_CONST const);

namespace Carpet {

using namespace std;

static void ProlongateGroupBoundaries(const cGH *cctkGH,
                                      const vector<int> &groups);

// Carpet's overload function for CCTK_SyncGroupsByDirI()
// which synchronises a set of groups given by their indices.
// Synchronisation of individual directions is not (yet?) implemented.
//
// returns the number of groups successfully synchronised
//        -1 if a group in the set doesn't have storage assigned
int SyncGroupsByDirI(const cGH *cctkGH, int num_groups, const int *groups,
                     const int *directions) {
  int group, retval = 0;
  vector<int> groups_set;
  groups_set.reserve(num_groups);

  // individual directions aren't supported (yet?)
  if (directions != NULL) {
    CCTK_ERROR(
        "Carpet doesn't support synchronisation of individual directions");
  }

  for (group = 0; group < num_groups; group++) {
    if (CCTK_NumVarsInGroupI(groups[group]) > 0) {
      groups_set.push_back(groups[group]);
    }
  }

  if (groups_set.size() > 0) {
    retval = SyncProlongateGroups(cctkGH, groups_set);
    if (retval == 0) {
      retval = groups_set.size();
    }
  }

  return retval;
}


// synchronises ghostzones and prolongates boundaries of a set of groups
//
// returns 0 for success and -1 if the set contains a group with no storage
int SyncProlongateGroups(const cGH *cctkGH, const vector<int> &groups,
                         cFunctionData const *function_data) {
  int retval = 0;
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  assert(groups.size() > 0);

  // check consistency of all groups:
  // create a new set with empty and no-storage groups removed
  vector<int> goodgroups;
  goodgroups.reserve(groups.size());
  for (size_t group = 0; group < groups.size(); group++) {
    const int g = groups.AT(group);
    const int grouptype = CCTK_GroupTypeI(g);
    char *const groupname = CCTK_GroupName(g);
    Checkpoint("SyncGroup \"%s\" iteration=%d time=%g", groupname,
               cctkGH->cctk_iteration, (double)cctkGH->cctk_time);

    if (grouptype == CCTK_GF) {
      if (reflevel == -1) {
        CCTK_VERROR("Cannot synchronise in global mode "
                   "(Tried to synchronise group \"%s\")",
                   groupname);
      }
      if (map != -1 and component == -1) {
        if (maps == 1) {
          CCTK_VWARN(2, "Synchronising group \"%s\" in singlemap mode",
                     groupname);
        } else {
          CCTK_VERROR("Cannot synchronise in singlemap mode "
                     "(Tried to synchronise group \"%s\")",
                     groupname);
        }
      }
      if (component != -1) {
        if (maps == 1 and vhh.AT(map)->local_components(reflevel) == 1) {
          CCTK_VWARN(2, "Synchronising group \"%s\" in local mode", groupname);
        } else {
          CCTK_VERROR("Cannot synchronise in local mode "
                     "(Tried to synchronise group \"%s\")",
                     groupname);
        }
      }
    }

    if (not CCTK_QueryGroupStorageI(cctkGH, g)) {
      CCTK_VWARN(4, "Cannot synchronise group \"%s\" because it has no storage",
                 groupname);
      retval = -1;
    } else if (CCTK_NumVarsInGroupI(g) > 0) {
      goodgroups.push_back(g);
    }

    free(groupname);
  } // for g

  if (goodgroups.size() > 0) {

  // TODO: move this into Timing.cc
  for (size_t gp = 0; gp < goodgroups.size(); gp++) {
    const int g = goodgroups.AT(gp);
    *syncs_count += CCTK_NumVarsInGroupI(g);
  }
#ifdef REQUIREMENTS_HH
    Requirements::Sync(function_data, goodgroups, cctkGH->cctk_iteration,
                       reflevel, timelevel);
#endif

    // prolongate boundaries
    bool const local_do_prolongate = do_prolongate and not do_taper;
    if (local_do_prolongate) {
      if (sync_barriers) {
        static Timers::Timer barrier_timer("PreProlongateBarrier");
        barrier_timer.start();
        CCTK_Barrier(cctkGH);
        barrier_timer.stop();
      }

      static Timers::Timer timer("Prolongate");
      timer.start();
      ProlongateGroupBoundaries(cctkGH, goodgroups);
      NoteTaperHaloRefilled();
      timer.stop();

      if (sync_barriers) {
        static Timers::Timer barrier_timer("PostProlongateBarrier");
        barrier_timer.start();
        CCTK_Barrier(cctkGH);
        barrier_timer.stop();
      }
    }

// This was found to be necessary on Hopper, otherwise memory
// seems to be overwritten while prolongating/synchronizing. It
// looks like this might be an MPI implementation issue, but
// this is not clear. A barrier at this point seems to be a
// sufficient workaround, and is now used on Hopper. For more
// information about this ask Frank Loeffler
// <knarf@cct.lsu.edu>.
#ifdef CARPET_MPI_BARRIER_PROLONGATE_SYNC
    Carpet::NamedBarrier(cctkGH, 8472211063,
                         "CARPET_MPI_BARRIER_PROLONGATE_SYNC");
#endif

    // synchronise ghostzones
    if (sync_during_time_integration or local_do_prolongate) {
      if (sync_barriers) {
        static Timers::Timer barrier_timer("PreSyncBarrier");
        barrier_timer.start();
        CCTK_Barrier(cctkGH);
        barrier_timer.stop();
      }

      static Timers::Timer timer("Sync");
      timer.start();
      SyncGroups(cctkGH, goodgroups);
      timer.stop();

      if (sync_barriers) {
        static Timers::Timer barrier_timer("PostSyncBarrier");
        barrier_timer.start();
        CCTK_Barrier(cctkGH);
        barrier_timer.stop();
      }
    }
    if(not CCTK_EQUALS(presync_mode, "off") and
       not CCTK_EQUALS(presync_mode, "warn-only")) {
      for(auto group : goodgroups) {
        ApplyPhysicalBCsForGroupI(cctkGH,group);
      }
    }
  }

  return retval;
}

static bool prolongates_time(operator_type const op) {
  return op != op_none and op != op_sync and op != op_restrict;
}

// Parent time level whose timestamp is this AB3 right-hand side. Both
// clocks use the same closed form, so the comparison is exact. A single
// stored time, or every stored time equal to the requested time, is time
// level 0: at startup the older levels can share that timestamp and still
// hold poison. Several matches with any other timestamp are fatal. No
// match at the step end is the AB3 update itself: the parent has no value
// there, and the right-hand side already ran at the shared start. No match
// at any other time is fatal. Returns -1 when prolongation should be
// skipped.
static int ab3_parent_timelevel(const cGH *const cctkGH, int const g,
                                ggf *const gv) {
  if (tt == NULL)
    CCTK_ERROR("AB3-RK4 prolongation has no Carpet time hierarchy.");
  int const parent = reflevel - 1;
  if (parent < 0)
    CCTK_ERROR("AB3-RK4 prolongation has no parent level.");
  int const parent_active =
      groupdata.AT(g).activetimelevels.AT(mglevel).AT(parent);
  if (parent_active <= 1)
    return 0;
  if (not prolongates_time(gv->transport_operator))
    return 0;
  if (parent_active > tt->timelevels) {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                "AB3-RK4 group %d has %d active parent time levels and the "
                "time hierarchy has %d.",
                g, parent_active, tt->timelevels);
    return -1;
  }

  CCTK_REAL const time = cctkGH->cctk_time;
  int nmatch = 0;
  int matched = -1;
  bool all_equal = true;
  for (int tl = 0; tl < parent_active; ++tl) {
    CCTK_REAL const stored = tt->get_time(mglevel, parent, tl);
    if (stored != time)
      all_equal = false;
    else {
      if (nmatch == 0)
        matched = tl;
      ++nmatch;
    }
  }
  if (all_equal)
    return 0;
  if (nmatch == 1) {
    if (matched < 0 or matched >= parent_active) {
      CCTK_ERROR("AB3-RK4 matched a parent time level this grid function "
                 "does not store.");
    }
    return matched;
  }

  ostringstream times;
  times.setf(ios::scientific);
  times.precision(17);
  for (int tl = 0; tl < parent_active; ++tl) {
    if (tl > 0)
      times << ", ";
    times << "tl" << tl << "=" << tt->get_time(mglevel, parent, tl);
  }
  char *const fullname = CCTK_FullName(gv->varindex);
  string const name = fullname ? fullname : "<unknown variable>";
  free(fullname);
  if (nmatch > 1) {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                "AB3-RK4 found %d parent times equal to %.17g for %s on "
                "level %d, and the stored times are not all that time (%s).",
                nmatch, static_cast<double>(time), name.c_str(), parent,
                times.str().c_str());
    return -1;
  }
  CCTK_REAL const level_now = tt->get_time(mglevel, reflevel, 0);
  if (time == level_now)
    return -1;
  CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
              "AB3-RK4 found no parent time equal to %.17g for %s on level "
              "%d. This level's current time is %.17g. Stored parent times: "
              "%s.",
              static_cast<double>(time), name.c_str(), parent,
              static_cast<double>(level_now), times.str().c_str());
  return -1;
}

// Prolongate the boundaries of all CCTK_GF groups in the given set
static void ProlongateGroupBoundaries(const cGH *cctkGH,
                                      const vector<int> &groups) {
  DECLARE_CCTK_PARAMETERS;

  if (reflevel == 0)
    return;

  Checkpoint("ProlongateGroups");

  assert(groups.size() > 0);

  // TODO: Transmit only the regions that are actually needed
  if (CCTK_IsFunctionAliased("Accelerator_RequireValidData")) {
    // TODO: Require only as many time levels as needed
    // TODO: Handle time levels correctly
    vector<CCTK_INT> vis, rls, tls;
    for (int group = 0; group < (int)groups.size(); ++group) {
      const int g = groups.AT(group);
      const int grouptype = CCTK_GroupTypeI(g);
      if (grouptype != CCTK_GF) {
        continue;
      }

      const int alloc_num_tl = CCTK_ActiveTimeLevelsGI(cctkGH, g);
      const int num_tl = prolongation_order_time + 1;
      if(num_tl > alloc_num_tl) {
        char *groupname = CCTK_GroupName(g);
        CCTK_VERROR("Cannot prolongate %s because only %d time levels have been allocated and %d are needed.",
                    groupname, alloc_num_tl, num_tl);
        free(groupname);
      }
      const int var0 = CCTK_FirstVarIndexI(g);
      const int varn = CCTK_NumVarsInGroupI(g);
      for (int vi = var0; vi < var0 + varn; ++vi) {
        vis.push_back(vi);
        rls.push_back(reflevel);
        tls.push_back(timelevel);
        for (int tl = 0; tl < num_tl; ++tl) {
          vis.push_back(vi);
          rls.push_back(reflevel - 1);
          tls.push_back(tl);
        }
      }
    }
    const CCTK_INT on_device = 0;
    Accelerator_RequireValidData(cctkGH, &vis.front(), &rls.front(),
                                 &tls.front(), vis.size(), on_device);
  }

  // use the current time here (which may be modified by the user)
  const CCTK_REAL time = cctkGH->cctk_time;

  static vector<Timers::Timer *> timers;
  if (timers.empty()) {
    timers.push_back(new Timers::Timer("comm_state[0].create"));
    for (astate state = static_cast<astate>(0); state != state_done;
         state = static_cast<astate>(static_cast<int>(state) + 1)) {
      ostringstream name1;
      name1 << "comm_state[" << timers.size() << "]"
            << "." << tostring(state) << ".user";
      timers.push_back(new Timers::Timer(name1.str()));
      ostringstream name2;
      name2 << "comm_state[" << timers.size() << "]"
            << "." << tostring(state) << ".step";
      timers.push_back(new Timers::Timer(name2.str()));
    }
  }

  vector<Timers::Timer *>::iterator ti = timers.begin();
  (*ti)->start();
  // The short step copies one parent time level. HermiteFill would read
  // the parent's current and previous slots, and a zero Lagrange weight
  // on a poisoned slot is NaN. The timer walk stays the same length.
  bool const ab3_step = csamr_ab3_rk4 and hermite_prolongate and
                        Carpet_AB3RK4_ThisStepAB3(cctkGH);
  if (ab3_step) {
    for (comm_state state; not state.done(); state.step()) {
      (*ti)->stop();
      ++ti;
      (*ti)->start();
      for (int group = 0; group < (int)groups.size(); ++group) {
        const int g = groups.AT(group);
        const int grouptype = CCTK_GroupTypeI(g);
        if (grouptype != CCTK_GF) {
          continue;
        }
        assert(reflevel >= 0 and reflevel < reflevels);
        const int active_tl =
            groupdata.AT(g).activetimelevels.AT(mglevel).AT(reflevel);
        assert(active_tl >= 0);
        const int tl = active_tl > 1 ? timelevel : 0;

        for (int m = 0; m < (int)arrdata.AT(g).size(); ++m) {
          for (int v = 0; v < (int)arrdata.AT(g).AT(m).data.size(); ++v) {
            ggf *const gv = arrdata.AT(g).AT(m).data.AT(v);
            int const src_tl = ab3_parent_timelevel(cctkGH, g, gv);
            if (src_tl < 0)
              continue;
            gv->ref_bnd_prolongate_spatial(state, tl, reflevel, mglevel,
                                           src_tl, time);
          }
        }
      }
      (*ti)->stop();
      ++ti;
      (*ti)->start();
    }
  } else {
    HermiteFill hermite_fill(cctkGH, groups, time);
    for (comm_state state; not state.done(); state.step()) {
      (*ti)->stop();
      ++ti;
      (*ti)->start();
      for (int group = 0; group < (int)groups.size(); ++group) {
        const int g = groups.AT(group);
        const int grouptype = CCTK_GroupTypeI(g);
        if (grouptype != CCTK_GF) {
          continue;
        }
        assert(reflevel >= 0 and reflevel < reflevels);
        const int active_tl =
            groupdata.AT(g).activetimelevels.AT(mglevel).AT(reflevel);
        assert(active_tl >= 0);
        const int tl = active_tl > 1 ? timelevel : 0;

        for (int m = 0; m < (int)arrdata.AT(g).size(); ++m) {
          for (int v = 0; v < (int)arrdata.AT(g).AT(m).data.size(); ++v) {
            ggf *const gv = arrdata.AT(g).AT(m).data.AT(v);
            int src_tl = -1;
            if (hermite_fill.source_timelevel(gv, src_tl)) {
              gv->ref_bnd_prolongate_spatial(state, tl, reflevel, mglevel,
                                             src_tl, time);
            } else {
              gv->ref_bnd_prolongate_all(state, tl, reflevel, mglevel, time);
            }
          }
        }
      }
      (*ti)->stop();
      ++ti;
      (*ti)->start();
    }
  }
  (*ti)->stop();
  ++ti;
  assert(ti == timers.end());

  // TODO: Transmit only the regions that were actually written
  // (i.e. the buffer zones)
  if (CCTK_IsFunctionAliased("Accelerator_NotifyDataModified")) {
    // TODO: Handle time levels correctly
    vector<CCTK_INT> vis, rls, tls;
    for (int group = 0; group < (int)groups.size(); ++group) {
      const int g = groups.AT(group);
      const int grouptype = CCTK_GroupTypeI(g);
      if (grouptype != CCTK_GF) {
        continue;
      }
      const int var0 = CCTK_FirstVarIndexI(g);
      const int varn = CCTK_NumVarsInGroupI(g);
      for (int vi = var0; vi < var0 + varn; ++vi) {
        vis.push_back(vi);
        rls.push_back(reflevel);
        tls.push_back(timelevel);
      }
    }
    const CCTK_INT on_device = 0;
    Accelerator_NotifyDataModified(cctkGH, &vis.front(), &rls.front(),
                                   &tls.front(), vis.size(), on_device);
  }
}

// synchronises a set of groups
void SyncGroups(const cGH *cctkGH, const vector<int> &groups) {
  DECLARE_CCTK_PARAMETERS;

  Checkpoint("SyncGroups");

  assert(groups.size() > 0);

  if (CCTK_IsFunctionAliased("Accelerator_PreSync")) {
    vector<CCTK_INT> groups_(groups.size());
    for (size_t i = 0; i < groups.size(); ++i)
      groups_[i] = groups[i];
    Accelerator_PreSync(cctkGH, &groups_.front(), groups_.size());
  }

  static vector<Timers::Timer *> timers;
  if (timers.empty()) {
    timers.push_back(new Timers::Timer("comm_state[0].create"));
    for (astate state = static_cast<astate>(0); state != state_done;
         state = static_cast<astate>(static_cast<int>(state) + 1)) {
      ostringstream name1;
      name1 << "comm_state[" << timers.size() << "]"
            << "." << tostring(state) << ".user";
      timers.push_back(new Timers::Timer(name1.str()));
      ostringstream name2;
      name2 << "comm_state[" << timers.size() << "]"
            << "." << tostring(state) << ".step";
      timers.push_back(new Timers::Timer(name2.str()));
    }
  }

  vector<Timers::Timer *>::iterator ti = timers.begin();
  (*ti)->start();
  for (comm_state state; not state.done(); state.step()) {
    (*ti)->stop();
    ++ti;
    (*ti)->start();
    for (int group = 0; group < (int)groups.size(); ++group) {
      const int g = groups.AT(group);
      const int grouptype = CCTK_GroupTypeI(g);
      const int ml = grouptype == CCTK_GF ? mglevel : 0;
      const int rl = grouptype == CCTK_GF ? reflevel : 0;
      const int active_tl = groupdata.AT(g).activetimelevels.AT(ml).AT(rl);
      assert(active_tl >= 0);
      const int tl = active_tl > 1 ? timelevel : 0;
      for (int m = 0; m < (int)arrdata.AT(g).size(); ++m) {
        for (int v = 0; v < (int)arrdata.AT(g).AT(m).data.size(); ++v) {
          arrdesc &array = arrdata.AT(g).AT(m);
          array.data.AT(v)->sync_all(state, tl, rl, ml);
        }
      }
    }
    (*ti)->stop();
    ++ti;
    (*ti)->start();
  }
  (*ti)->stop();
  ++ti;
  assert(ti == timers.end());

  if (CCTK_IsFunctionAliased("Accelerator_PostSync")) {
    vector<CCTK_INT> groups_(groups.size());
    for (size_t i = 0; i < groups.size(); ++i)
      groups_[i] = groups[i];
    Accelerator_PostSync(cctkGH, &groups_.front(), groups_.size());
  }
}

int EnableGroupComm(const cGH *cctkGH, const char *groupname) {
  // Communication is always enabled
  return 0;
}

int DisableGroupComm(const cGH *cctkGH, const char *groupname) {
  // Communication is always enabled
  return -1;
}

/* The caller may be inside one component. Leave that mode, sync every
   component of this refinement level, then restore the caller. */
extern "C" CCTK_INT Carpet_SyncVarIInLevelMode(CCTK_POINTER_TO_CONST cctkGH_,
                                               CCTK_INT varindex) {
  cGH const *const cctkGH = static_cast<cGH const *>(cctkGH_);
  if (!cctkGH) {
    CCTK_ERROR("Driver_SyncVarIInLevelMode has no grid hierarchy");
  }
  if (varindex < 0) {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                "Driver_SyncVarIInLevelMode has variable index %d",
                (int)varindex);
  }
  if (is_meta_mode() or is_global_mode()) {
    CCTK_ERROR("Driver_SyncVarIInLevelMode needs a refinement level");
  }
  CCTK_INT ierr = -1;
  BEGIN_LEVEL_MODE(cctkGH) {
    ierr = ::CCTK_SyncGroupWithVarI(cctkGH, varindex);
  }
  END_LEVEL_MODE;
  return ierr;
}

} // namespace Carpet
