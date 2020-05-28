#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <cctk_Functions.h>

#include <cmath>
#include <cassert>
#include <cstring>

#include <iostream>
#include <locale>
#include <map>
#include <set>
#include <tuple>
#include <vector>
#include <sstream>

#include "variables.hh"
#include "modes.hh"
#include "carpet.hh"
#include "PreSync.hh"
#include "Carpet_Prototypes.h"

#include <Timer.hh>

#undef PRESYNC_DEBUG

namespace Carpet {

// local functions and types
namespace {

// boundary condition handling
struct Bounds {
  CCTK_INT selected_faces;
  struct Bound {
    std::string bc_name;
    CCTK_INT faces;
    CCTK_INT width;
    CCTK_INT table_handle;
    Bound() : faces(0), width(-1), table_handle(-1) {}
  };
  std::vector<Bound> bounds;
  Bounds() : selected_faces(0), bounds() {}
};

// these keep track of selected boundary condistions
std::map<int,Bounds> boundary_conditions;

// set to true while inside of Carpet's own BC function, used to report / no
// report to thorn Boundary which BC it should ignore
bool do_applyphysicalbcs = false;

/**
 *  helpers to handle valid states
 */

/**
 * check flags
 */
inline bool is_set(int flags,int flag) {
  return (flags & flag) == flag;
}

/**
 * Provide a string representation of a variable's name
 */
std::string format_var_tuple(const int vi, const int rl, const int tl) {
  std::ostringstream os;
  os << CCTK_FullVarName(vi);
  for(int i=0;i<tl;i++)
    os << "_p";
  os << " (rl=" << rl << ")";
  return os.str();
}

/**
 * Provide a string representation of the code for the boundary.
 */
std::string format_where(const int wh) {
  if(wh == CCTK_VALID_EVERYWHERE) {
    return "Everywhere";
  }
  if(wh == CCTK_VALID_NOWHERE) {
    return "Nowhere";
  }
  std::string s;
  if(is_set(wh,CCTK_VALID_INTERIOR))
    s += "Interior";
  if(is_set(wh,CCTK_VALID_BOUNDARY)) {
    if(s.size() > 0) s+= "With";
    s += "Boundary";
  }
  if(is_set(wh,CCTK_VALID_GHOSTS)) {
    if(s.size() > 0) s+= "With";
    s += "Ghosts";
  }
  return s;
}

/**
 * Provide a string representation of the curent information I have about a
 * variable
 */
std::string format_valid(const int vi) {
  std::ostringstream os;
  int const gi = CCTK_GroupIndexFromVarI(vi);
  assert(gi >= 0);
  int const var = vi - CCTK_FirstVarIndexI(gi);
  assert(var >= 0);

  int const map0 = 0;
  ggf *const ff = arrdata.AT(gi).AT(map0).data.AT(var);
  assert(ff);

  os << "\nValid entries:";
  for (int rl=0; rl<ff->h.reflevels(); ++rl) {
    for (int tl=0; tl<ff->timelevels(mglevel, rl); ++tl) {
      os << " " << format_var_tuple(vi, rl, tl);
      os << " " << format_where(ff->valid(mglevel, rl, tl));
    } // tl
  } // rl

  return os.str();
}

/**
 * Schedule helper
 */
int ScheduleTraverse(char const *const where, char const *const name,
                     cGH *const cctkGH) {
  Timers::Timer timer(name);
  timer.start();
  ostringstream infobuf;
  infobuf << "Scheduling " << name;
  string const info = infobuf.str();
  Checkpoint(info.c_str());
  int const ierr = CCTK_ScheduleTraverse(name, cctkGH, CallFunction);
  timer.stop();
  return ierr;
}

/**
 * translate a relative timelevel from a READS statement to and absolute timelevel of a ggf
 */
int get_timelevel(int const active_tl, int const tl) {
  assert(active_tl >= 0);
  int tl_offset;
  if (active_tl == 0) {
    // group has no storage
    tl_offset = 0;
  } else if (do_allow_past_timelevels) {
    // regular case; all timelevels are accessible
    tl_offset = 0;
  } else {
    // only one timelevel is accessible
    if (timelevel < active_tl) {
      // timelevel "timelevel" exists
      tl_offset = timelevel;
    } else {
      // timelevel "timelevel" does not exist
      tl_offset = active_tl - 1;
    }
  }
  return tl_offset + tl;
}

} // namespace

/**********************************************************************
 ***********SCHEDULE READS / WRITES ***********************************
 **********************************************************************/

/**
 * Computes which groups need to be presync'd.
 */
void PreCheckValid(cFunctionData *attribute,cGH *cctkGH, std::vector<int>& pre_groups) {
  DECLARE_CCTK_PARAMETERS;

  assert(attribute);
  std::set<int> tmpgroups;

  assert(not CCTK_EQUALS(presync_mode, "off"));

  bool const warn = CCTK_EQUALS(presync_mode, "warn-only") or
                    CCTK_EQUALS(presync_mode, "mixed-warn");
  bool const error = CCTK_EQUALS(presync_mode, "mixed-error") or
                     CCTK_EQUALS(presync_mode, "presync-only");
  assert(warn + error == 1);
  bool const may_sync = not CCTK_EQUALS(presync_mode, "warn-only");

  for(int i=0;i<attribute->n_RDWR;i++) {
    const RDWR_entry& entry = attribute->RDWR[i];

    assert(0 <= entry.varindex and entry.varindex < CCTK_NumVars());

    if(entry.where_rd == CCTK_VALID_NOWHERE) { // only READS cn trigger sync or errors
      continue;
    }

    if(not do_allow_past_timelevels and entry.timelevel > 0) {
      /* must not abort or do anything since Carpet passed NULL pointers to
       * thorns who can detect invalid past level */
      continue;
    }

    int const group_index = CCTK_GroupIndexFromVarI(entry.varindex);
    assert(group_index >= 0);
    const int type = CCTK_GroupTypeI(group_index);
    const int rl = type == CCTK_GF ? reflevel : 0;

    // some thorns only read from variables based on extra parameters, eg
    // TmunuBase's stress_energy_state grid scalar. This check mimics the
    // behaviour of calling SYNC on a variable without storage by outputting a
    // high level warning. A thorn that actually tries to read from the
    // variable will encounter a SEGFAULT.
    int const active_tl =
      groupdata.AT(group_index).activetimelevels.AT(mglevel).AT(rl);
    int const tl = get_timelevel(active_tl, entry.timelevel);

    if(tl >= active_tl) {
      CCTK_VWARN(CCTK_WARN_DEBUG,
                 "Declared access to '%s' on time level %d which has no storage",
                 CCTK_FullVarName(entry.varindex), tl);
      continue;
    }

    int const var = entry.varindex - CCTK_FirstVarIndexI(group_index);
    int const map0 = 0;
    ggf *const ff = arrdata.AT(group_index).AT(map0).data.AT(var);
    assert(ff);

    int const valid = ff->valid(mglevel, rl, tl);

    if(not is_set(valid, entry.where_rd)) {
      // we have: interior required-for ghosts required-for boundary

      if((is_set(entry.where_rd, CCTK_VALID_GHOSTS) and
          not is_set(valid, CCTK_VALID_GHOSTS)) or
         (is_set(entry.where_rd, CCTK_VALID_BOUNDARY) and
          not is_set(valid, CCTK_VALID_BOUNDARY))) {

        // give warning / error if we cannot make things ok by SYNCing
        if(entry.timelevel > 0 or not is_set(valid, CCTK_VALID_INTERIOR) or
           not may_sync) {
          std::ostringstream msg;
          msg << "Required read for "
              << format_var_tuple(entry.varindex, rl, entry.timelevel)
              << " not satisfied. Have " << format_where(valid)
              << " and require " << format_where(entry.where_rd) << " missing "
              << format_where(~valid & entry.where_rd)
              << " at the start of routine "
              << attribute->thorn << "::" << attribute->routine << ".";
          if(entry.timelevel > 0) {
            msg << " Cannot SYNC past timelevels.";
          }
          if(error) {
            msg << " Current valid state: " << format_valid(entry.varindex) << ".";
          }
          CCTK_WARN(warn ? CCTK_WARN_ALERT : CCTK_WARN_ABORT, msg.str().c_str());
        } else {
          assert(is_set(valid, CCTK_VALID_INTERIOR));
          assert(not is_set(valid, CCTK_VALID_GHOSTS) or
                 not is_set(valid, CCTK_VALID_BOUNDARY));
          tmpgroups.insert(group_index);
        }
      } // invalid needed GHOSTS or BOUNDARY
    } // no everything valie
  } // for RDWR

  pre_groups.assign(tmpgroups.begin(), tmpgroups.end());
}

/**
 * Do the actual presync of the groups.
 **/
void PreSyncGroups(cFunctionData *attribute,cGH *cctkGH,const std::vector<int>& pre_groups) {
  DECLARE_CCTK_PARAMETERS;

  assert(not CCTK_EQUALS(presync_mode, "off") and
         not CCTK_EQUALS(presync_mode, "warn-only"));

  if(pre_groups.empty())
    return;

  if(reflevel > 0) {
    // recurse to check that all coarsers levels are properly SYNCed
    CCTK_REAL previous_time = cctkGH->cctk_time;
    const int parent_reflevel = reflevel - 1;
    BEGIN_GLOBAL_MODE(cctkGH) {
      ENTER_LEVEL_MODE(cctkGH, parent_reflevel) {
        cctkGH->cctk_time = tt->get_time(mglevel, reflevel, timelevel);
        PreSyncGroups(attribute, cctkGH, pre_groups);
      } LEAVE_LEVEL_MODE;
    } END_GLOBAL_MODE;
    cctkGH->cctk_time = previous_time;
  }

  // ask Carpet to do the SYNC, this will apply BC as well
  SyncProlongateGroups(cctkGH, pre_groups, attribute);
}

/**
 * after a scheduled routined finished, update the valid states
 */
void PostCheckValid(cFunctionData *attribute, cGH *cctkGH) {
  DECLARE_CCTK_PARAMETERS;

  assert(not CCTK_EQUALS(presync_mode, "off"));

  for (int i=0;i<attribute->n_RDWR;i++) {
    const RDWR_entry& entry = attribute->RDWR[i];

    assert(0 <= entry.varindex and entry.varindex < CCTK_NumVars());

    if(entry.where_wr == CCTK_VALID_NOWHERE) { // nothing to do
      continue;
    }

    if(not do_allow_past_timelevels and entry.timelevel > 0) {
      /* must not abort or do anything since Carpet passed NULL pointers to
       * thorns who can detect invalid past level */
      continue;
    }

    int const group_index = CCTK_GroupIndexFromVarI(entry.varindex);
    assert(group_index >= 0);
    const int type = CCTK_GroupTypeI(group_index);
    const int rl = type == CCTK_GF ? reflevel : 0;

    // some thorns only read from variables based on extra parameters, eg
    // TmunuBase's stress_energy_state grid scalar. This check mimics the
    // behaviour of calling SYNC on a variable without storage by outputting a
    // high level warning. A thorn that actually tries to read from the
    // variable will encounter a SEGFAULT.
    int const active_tl =
      groupdata.AT(group_index).activetimelevels.AT(mglevel).AT(rl);
    int const tl = get_timelevel(active_tl, entry.timelevel);

    if(tl >= active_tl) {
      CCTK_VWARN(CCTK_WARN_DEBUG,
                 "Declared access to '%s' on time level %d which has no storage",
                 CCTK_FullVarName(entry.varindex), tl);
      continue;
    }

    int const var = entry.varindex - CCTK_FirstVarIndexI(group_index);
    int const map0 = 0;
    ggf *const ff = arrdata.AT(group_index).AT(map0).data.AT(var);
    assert(ff);

    if(entry.where_wr == CCTK_VALID_INTERIOR) {
      ff->set_valid(mglevel, rl, tl, CCTK_VALID_INTERIOR);
    } else {
      int const old_valid = ff->valid(mglevel, rl, tl);
      ff->set_valid(mglevel, rl, tl, old_valid | entry.where_wr);
    }
  }
}

/**********************************************************************
 **************Programmatic READS / WRITES ****************************
 **********************************************************************/

/**
 * Apply manual synchronization to some variables. If these variables is
 * already valid in the requested regions , this routine does nothing. When the
 * routine finishes, it will be valid in the requested regions.
 */
namespace {
CCTK_INT RequireValidData(const cGH* cctkGH,
   CCTK_INT const * const variables, CCTK_INT const * const tls,
   CCTK_INT const nvariables, CCTK_INT const * wheres) {
  DECLARE_CCTK_PARAMETERS;

  // TODO: this is almost the same as the READS checks only, consider factoring
  // out common code in particular the warning logic.

  assert(not CCTK_EQUALS(presync_mode, "off"));

  assert(variables or nvariables == 0);
  assert(tls or nvariables == 0);
  assert(wheres or nvariables == 0);

  assert(is_level_mode());

  if(reflevel > 0) {
    cGH* nonconstGH = const_cast<cGH*>(cctkGH);
    // recurse to check that all coarsers levels are properly SYNCed
    CCTK_REAL previous_time = nonconstGH->cctk_time;
    const int parent_reflevel = reflevel - 1;
    BEGIN_GLOBAL_MODE(nonconstGH) {
      ENTER_LEVEL_MODE(nonconstGH, parent_reflevel) {
        nonconstGH->cctk_time = tt->get_time(mglevel, reflevel, timelevel);
        RequireValidData(nonconstGH, variables, tls, nvariables, wheres);
      } LEAVE_LEVEL_MODE;
    } END_GLOBAL_MODE;
    nonconstGH->cctk_time = previous_time;
  }

  bool const warn = CCTK_EQUALS(presync_mode, "warn-only") or
                    CCTK_EQUALS(presync_mode, "mixed-warn");
  bool const error = CCTK_EQUALS(presync_mode, "mixed-error") or
                     CCTK_EQUALS(presync_mode, "presync-only");
  assert(warn + error == 1);
  bool const presync_only = CCTK_EQUALS(presync_mode, "presync-only");
  bool const may_sync = not CCTK_EQUALS(presync_mode, "warn-only");

  std::set<std::tuple<int,int>> processed_groups; // [gi] [tl]
  for(int i = 0; i < nvariables; i++) {
    int const vi = variables[i];
    int const where = wheres[i];

    assert(0 <= vi and vi < CCTK_NumVars());

    if(vi < 0 or vi >= CCTK_NumVars()) {
      CCTK_VERROR("Invalid variable index %d", vi);
    }

    if(not do_allow_past_timelevels and tls[i] > 0) {
      /* must not abort or do anything since Carpet passed NULL pointers to
       * thorns who can detect invalid past level */
      continue;
    }

    int const group_index = CCTK_GroupIndexFromVarI(vi);
    assert(group_index >= 0);
    const int type = CCTK_GroupTypeI(group_index);
    const int rl = type == CCTK_GF ? reflevel : 0;

    // some thorns only read from variables based on extra parameters, eg
    // TmunuBase's stress_energy_state grid scalar. This check mimics the
    // behaviour of calling SYNC on a variable without storage by outputting a
    // high level warning. A thorn that actually tries to read from the
    // variable will encounter a SEGFAULT.
    int const active_tl =
      groupdata.AT(group_index).activetimelevels.AT(mglevel).AT(rl);
    int const tl = get_timelevel(active_tl, tls[i]);

    if(tl >= active_tl) {
      CCTK_VWARN(CCTK_WARN_DEBUG,
                 "Declared access to '%s' on time level %d which has no storage",
                 CCTK_FullVarName(vi), tl);
      continue;
    }

    int const var = vi - CCTK_FirstVarIndexI(group_index);
    int const map0 = 0;
    ggf *const ff = arrdata.AT(group_index).AT(map0).data.AT(var);
    assert(ff);

    int const valid = ff->valid(mglevel, rl, tl);

    if(not is_set(valid, where)) {
      // we have: interior required-for ghosts required-for boundary

      bool const bc_selected = QueryDriverBCForVarI(cctkGH, vi);

      if((is_set(where, CCTK_VALID_GHOSTS) and
          not is_set(valid, CCTK_VALID_GHOSTS)) or
         (is_set(where, CCTK_VALID_BOUNDARY) and
          not is_set(valid, CCTK_VALID_BOUNDARY))) {

        // give warning / error if we cannot make things ok by SYNCing
        if(tls[i] > 0 or not is_set(valid, CCTK_VALID_INTERIOR) or not may_sync) {
          /* only warn / error about functions that the client thorn told me
           * about one way or the other */
          if(not presync_only and not bc_selected)
            continue;

          std::ostringstream msg;
          msg << "Required read for "
              << format_var_tuple(vi, rl, tls[i])
              << " not satisfied. Have " << format_where(valid) << " and require "
              << format_where(where) << " missing "
              << format_where(~valid & where) << ".";
          if(tl > 0) {
            msg << " Cannot SYNC past timelevels.";
          }
          if(error) {
            msg << " Current valid state: " << format_valid(vi) << ".";
          }
          CCTK_WARN(warn ? CCTK_WARN_ALERT : CCTK_WARN_ABORT, msg.str().c_str());
        } else if(bc_selected and may_sync) {
          // we need to ad can SYNC

          assert(is_set(valid, CCTK_VALID_INTERIOR));
          assert(not is_set(valid, CCTK_VALID_GHOSTS) or
                 not is_set(valid, CCTK_VALID_BOUNDARY));

          int const gi = CCTK_GroupIndexFromVarI(vi);

          // since SYNC processes whole groups, if any variable is unset then we must
          // not have encountered any other variable from that group either
          auto res = processed_groups.insert(std::tuple<int,int>(gi,tl));
          assert(res.second); // tuple was not marked processed before

          // actually apply SYNC
          std::vector<int> sync_groups;
          sync_groups.push_back(gi);
          const cFunctionData *attribute = CCTK_ScheduleQueryCurrentFunction(cctkGH);
          // copied out of modes.hh
          int const old_timelevel_offset = timelevel_offset;
          int const old_timelevel = timelevel;
          CCTK_REAL old_cctk_time = cctkGH->cctk_time;
          auto const old_do_allow_past_timelevels = do_allow_past_timelevels;
          do_allow_past_timelevels = tl == 0;
          timelevel_offset = timelevel = tl;
          const_cast<cGH*>(cctkGH)->cctk_time = tt->get_time(mglevel, rl, timelevel);
          int ierr = SyncProlongateGroups(cctkGH, sync_groups, attribute);
          assert(not ierr);
          const_cast<cGH*>(cctkGH)->cctk_time = old_cctk_time;
          timelevel_offset = old_timelevel_offset;
          timelevel = old_timelevel;
          do_allow_past_timelevels = old_do_allow_past_timelevels;

          int const new_valid = ff->valid(mglevel, rl, tl);
          assert(is_set(new_valid, where));
        } // bs_selected and may_sync
      } // invalid needed GHOSTS or BOUNDARY
    } // no already valid
  } // for nvariables

  return 0;
}
}

extern "C"
CCTK_INT Carpet_RequireValidData(CCTK_POINTER_TO_CONST cctkGH_,
   CCTK_INT const * const variables, CCTK_INT const * const tls,
   CCTK_INT const nvariables, CCTK_INT const * wheres) {
  DECLARE_CCTK_PARAMETERS;

  const cGH *cctkGH = static_cast<const cGH*>(cctkGH_);

  if(CCTK_EQUALS(presync_mode, "off") or CCTK_EQUALS(presync_mode, "warn-only"))
    return 0;

  // Take action by mode
  if(is_singlemap_mode() or is_local_mode()) {
    BEGIN_LEVEL_MODE(cctkGH) {
      RequireValidData(cctkGH, variables, tls, nvariables, wheres);
    }
    END_LEVEL_MODE;
  } else if(is_level_mode()) {
    RequireValidData(cctkGH, variables, tls, nvariables, wheres);
  } else if(is_global_mode()) {
    BEGIN_REFLEVEL_LOOP(cctkGH) {
      RequireValidData(cctkGH, variables, tls, nvariables, wheres);
    }
    END_REFLEVEL_LOOP;
  } else if(is_meta_mode()) {
    BEGIN_MGLEVEL_LOOP(cctkGH) {
      BEGIN_REFLEVEL_LOOP(cctkGH) {
        RequireValidData(cctkGH, variables, tls, nvariables, wheres);
      }
      END_REFLEVEL_LOOP;
    }
    END_MGLEVEL_LOOP;
  } else {
    CCTK_VERROR("%s must be called in local, singlemap, level, global, or meta mode", __func__);
  }

  return 0;
}

namespace {
CCTK_INT NotifyDataModified(const cGH * cctkGH,
   CCTK_INT const * const variables, CCTK_INT const * const tls,
   CCTK_INT const nvariables, CCTK_INT const * wheres) {
  DECLARE_CCTK_PARAMETERS;

  // TODO: pretty much the same as the WRITE handling. Consder factoring out
  // common code in particular the warnings.

  assert(not CCTK_EQUALS(presync_mode, "off"));

  assert(variables or nvariables == 0);
  assert(tls or nvariables == 0);

  assert(is_level_mode());

  for(int i = 0; i < nvariables; i++) {
    int const vi = variables[i];
    int const where = wheres[i];

    if(vi < 0 or vi >= CCTK_NumVars()) {
      CCTK_VERROR("Invalid variable index %d", vi);
    }

    if(not do_allow_past_timelevels and tls[i] > 0) {
      /* must not abort or do anything since Carpet passed NULL pointers to
       * thorns who can detect invalid past level */
      continue;
    }

    int const group_index = CCTK_GroupIndexFromVarI(vi);
    assert(group_index >= 0);
    const int type = CCTK_GroupTypeI(group_index);
    const int rl = type == CCTK_GF ? reflevel : 0;

    // some thorns only read from variables based on extra parameters, eg
    // TmunuBase's stress_energy_state grid scalar. This check mimics the
    // behaviour of calling SYNC on a variable without storage by outputting a
    // high level warning. A thorn that actually tries to read from the
    // variable will encounter a SEGFAULT.
    int const active_tl =
      groupdata.AT(group_index).activetimelevels.AT(mglevel).AT(rl);
    int const tl = get_timelevel(active_tl, tls[i]);

    if(tl >= active_tl) {
      CCTK_VWARN(CCTK_WARN_DEBUG,
                 "Declared access to '%s' on time level %d which has no storage",
                 CCTK_FullVarName(vi), tl);
      continue;
    }

    int const var = vi - CCTK_FirstVarIndexI(group_index);
    int const map0 = 0;
    ggf *const ff = arrdata.AT(group_index).AT(map0).data.AT(var);
    assert(ff);

    if(where == CCTK_VALID_INTERIOR) {
      ff->set_valid(mglevel, rl, tl, CCTK_VALID_INTERIOR);
    } else {
      int const old_valid = ff->valid(mglevel, rl, tl);
      ff->set_valid(mglevel, rl, tl, old_valid | where);
    }
  }

  return 0;
}
}

extern "C"
CCTK_INT Carpet_NotifyDataModified(CCTK_POINTER_TO_CONST cctkGH_,
   CCTK_INT const * const variables, CCTK_INT const * const tls,
   CCTK_INT const nvariables, CCTK_INT const * wheres) {
  DECLARE_CCTK_PARAMETERS;

  const cGH *cctkGH = static_cast<const cGH*>(cctkGH_);

  if(CCTK_EQUALS(presync_mode, "off"))
    return 0;

  // Take action by mode
  if(is_singlemap_mode() or is_local_mode()) {
    BEGIN_LEVEL_MODE(cctkGH) {
      NotifyDataModified(cctkGH, variables, tls, nvariables, wheres);
    }
    END_LEVEL_MODE;
  } else if(is_level_mode()) {
    NotifyDataModified(cctkGH, variables, tls, nvariables, wheres);
  } else if(is_global_mode()) {
    BEGIN_REFLEVEL_LOOP(cctkGH) {
      NotifyDataModified(cctkGH, variables, tls, nvariables, wheres);
    }
    END_REFLEVEL_LOOP;
  } else if(is_meta_mode()) {
    BEGIN_MGLEVEL_LOOP(cctkGH) {
      BEGIN_REFLEVEL_LOOP(cctkGH) {
        NotifyDataModified(cctkGH, variables, tls, nvariables, wheres);
      }
      END_REFLEVEL_LOOP;
    }
    END_MGLEVEL_LOOP;
  } else {
    CCTK_VERROR("%s must be called in local, singlemap, level, global, or meta mode", __func__);
  }

  return 0;
}

/**********************************************************************
 ************** Low level access to valid state ************************
 **********************************************************************/

/**
 * Given a variable and a timelevel, set the region
 * of the grid where that variable is valid (i.e. the where_spec).
 */
extern "C"
void Carpet_SetValidRegion(CCTK_POINTER_TO_CONST /*cctkGH_*/, CCTK_INT vi,
                           CCTK_INT tl_in, CCTK_INT wh) {
  DECLARE_CCTK_PARAMETERS;

  if(vi < 0 or vi >= CCTK_NumVars()) {
    CCTK_VERROR("Invalid variable index %d", vi);
  }

  if(not do_allow_past_timelevels and tl_in > 0) {
      /* must not abort or do anything since Carpet passed NULL pointers to
       * thorns who can detect invalid past level */
    return;
  }

  int const group_index = CCTK_GroupIndexFromVarI(vi);
  assert(group_index >= 0);
  const int type = CCTK_GroupTypeI(group_index);
  const int rl = type == CCTK_GF ? reflevel : 0;

  // some thorns only read from variables based on extra parameters, eg
  // TmunuBase's stress_energy_state grid scalar. This check mimics the
  // behaviour of calling SYNC on a variable without storage by outputting a
  // high level warning. A thorn that actually tries to read from the
  // variable will encounter a SEGFAULT.
  int const active_tl =
    groupdata.AT(group_index).activetimelevels.AT(mglevel).AT(rl);
  int const tl = get_timelevel(active_tl, tl_in);

  if(tl >= active_tl) {
    CCTK_VWARN(CCTK_WARN_DEBUG,
               "Declared access to '%s' on time level %d which has no storage",
               CCTK_FullVarName(vi), tl);
    return;
  }

  int const var = vi - CCTK_FirstVarIndexI(group_index);
  int const map0 = 0;
  ggf *const ff = arrdata.AT(group_index).AT(map0).data.AT(var);
  assert(ff);

  ff->set_valid(mglevel, rl, tl, wh);
}

/**
 * Given a variable and a timelevel, return the region
 * of the grid where that variable is valid (i.e. the where_spec).
 */
extern "C"
CCTK_INT Carpet_GetValidRegion(CCTK_POINTER_TO_CONST /*cctkGH_*/, CCTK_INT vi,
                               CCTK_INT tl_in) {
  DECLARE_CCTK_PARAMETERS;

  if(vi < 0 or vi >= CCTK_NumVars()) {
    CCTK_VERROR("Invalid variable index %d", vi);
  }

  if(not do_allow_past_timelevels and tl_in > 0) {
      /* must not abort or do anything since Carpet passed NULL pointers to
       * thorns who can detect invalid past level */
    return CCTK_VALID_NOWHERE;
  }

  int const group_index = CCTK_GroupIndexFromVarI(vi);
  assert(group_index >= 0);
  const int type = CCTK_GroupTypeI(group_index);
  const int rl = type == CCTK_GF ? reflevel : 0;

  // some thorns only read from variables based on extra parameters, eg
  // TmunuBase's stress_energy_state grid scalar. This check mimics the
  // behaviour of calling SYNC on a variable without storage by outputting a
  // high level warning. A thorn that actually tries to read from the
  // variable will encounter a SEGFAULT.
  int const active_tl =
    groupdata.AT(group_index).activetimelevels.AT(mglevel).AT(rl);
  int const tl = get_timelevel(active_tl, tl_in);

  if(tl >= active_tl) {
    CCTK_VWARN(CCTK_WARN_DEBUG,
               "Declared access to '%s' on time level %d which has no storage",
               CCTK_FullVarName(vi), tl);
    return CCTK_VALID_NOWHERE;
  }

  int const var = vi - CCTK_FirstVarIndexI(group_index);
  int const map0 = 0;
  ggf *const ff = arrdata.AT(group_index).AT(map0).data.AT(var);
  assert(ff);

  return ff->valid(mglevel, rl, tl);
}

/**********************************************************************
 ************** Boundary condtions controlled by driver ****************
 **********************************************************************/

namespace {
// return codes are those of Boundary_SelectVarForBCI
CCTK_INT SelectVarForBCI(
    const cGH *cctkGH,
    const CCTK_INT faces,
    const CCTK_INT width,
    const CCTK_INT table_handle,
    const CCTK_INT var_index,
    const CCTK_STRING bc_name) {
  if(not Boundary_QueryRegisteredPhysicalBC(cctkGH, bc_name)) {
    CCTK_VWARN(CCTK_WARN_ALERT,
               "There is no function implementing the physical boundary condition %s",
               bc_name);
    return -2;
  }
  if(var_index < 0 or var_index >= CCTK_NumVars()) {
    CCTK_VWARN(CCTK_WARN_ALERT, "Invalid variable index");
    return -7;
  }

  Bounds& bounds = boundary_conditions[var_index];
  if(bounds.selected_faces & faces) {
    // in order to support DriverSelectVarForBC inside of an old-style SelectBC
    // scheduled function, accept multiple identical registrations
    for(auto bound: bounds.bounds) {
      if(bound.faces == faces and bound.width == width and
         bound.table_handle == table_handle and
         CCTK_EQUALS(bound.bc_name.c_str(), bc_name)) {
        return 0;
      }
    }
    std::ostringstream msg;
    msg << "'" << CCTK_FullVarName(var_index) << "'"
        << " has already been selected for a bc. Selected BCs:";
    for(auto bound: bounds.bounds) {
      msg << " " << bound.bc_name;
    }
    CCTK_VWARN(CCTK_WARN_ALERT, msg.str().c_str());
    return -3;
  }

  Bounds::Bound b;
  b.faces = faces;
  b.width = width;
  b.table_handle = table_handle;
  b.bc_name = bc_name;

  bounds.selected_faces |= faces;
  bounds.bounds.push_back(b);
  return 0;
}
}

/**
 * Choose which boundary condtion to apply when needed. These mimic Boundary's
 * functions.
 */
extern "C"
CCTK_INT Carpet_SelectVarForBC(
    const CCTK_POINTER_TO_CONST cctkGH_,
    const CCTK_INT faces,
    const CCTK_INT width,
    const CCTK_INT table_handle,
    const CCTK_STRING var_name,
    const CCTK_STRING bc_name) {
  const cGH *cctkGH = static_cast<const cGH*>(cctkGH_);
  const CCTK_INT vi = CCTK_VarIndex(var_name);
  return SelectVarForBCI(cctkGH,faces,width,table_handle,vi,bc_name);
}

extern "C"
CCTK_INT Carpet_SelectGroupForBC(
    const CCTK_POINTER_TO_CONST cctkGH_,
    const CCTK_INT faces,
    const CCTK_INT width,
    const CCTK_INT table_handle,
    const CCTK_STRING group_name,
    const CCTK_STRING bc_name) {
  const cGH *cctkGH = static_cast<const cGH*>(cctkGH_);
  const CCTK_INT group = CCTK_GroupIndex(group_name);
  if(group < 0)
    return -6; // see documenation from Boundary
  const CCTK_INT vstart = CCTK_FirstVarIndexI(group);
  const CCTK_INT vnum   = CCTK_NumVarsInGroupI(group);
  CCTK_INT ierr = 0;
  for(CCTK_INT vi=vstart;vi<vstart+vnum;vi++) {
    const CCTK_INT myierr =
      SelectVarForBCI(cctkGH,faces,width,table_handle,vi,bc_name);
    if(ierr == 0 and myierr != 0)
      ierr = myierr;
  }
  return ierr;
}

/**
 * A hook for Boundary to call. Carpet reports whether it would like Boundary
 * to now apply the selected boundary conditons. Usually false for grid
 * functions with a driver boundary conditions selected but true when Carpet
 * itself wants the boundary condtiosn applied.
 */
extern "C"
CCTK_INT Carpet_FilterOutVarForBCI(
    const CCTK_POINTER_TO_CONST cctkGH_,
    const CCTK_INT var_index) {
  DECLARE_CCTK_PARAMETERS;

  // do apply physical BC if we are being called from Carpet's own Driver BC
  // routine
  static bool const presync_only = CCTK_EQUALS(presync_mode, "presync-only");
  static bool const no_psync = CCTK_EQUALS(presync_mode, "off") or
                               CCTK_EQUALS(presync_mode, "warn-only");

  if(do_applyphysicalbcs)
    return false;
  else if(no_psync)
    return false;
  else if(presync_only)
    return true;

  return QueryDriverBCForVarI(static_cast<const cGH*>(cctkGH_), var_index);
}

/**
 * return true or false depending on whether a driver boundary condtion is
 * registered for this variable
 */
int QueryDriverBCForVarI(const cGH *cgh, const int varindex) {
  auto it = boundary_conditions.find(varindex);
  return it != boundary_conditions.end();
}

/**
 * Apply boundary conditions for a group. Called from inside SyncProlongateGroups.
 */
void ApplyPhysicalBCsForGroupI(const cGH *cctkGH, const int group_index) {
  DECLARE_CCTK_PARAMETERS;

  assert(not CCTK_EQUALS(presync_mode, "off") and
         not CCTK_EQUALS(presync_mode, "warn-only"));

  assert(0 <= group_index and group_index < CCTK_NumGroups());

  // some thorns only read from variables based on extra parameters, eg
  // TmunuBase's stress_energy_state grid scalar. This check mimics the
  // behaviour of calling SYNC on a variable without storage by outputting a
  // high level warning. A thorn that actually tries to read from the
  // variable will encounter a SEGFAULT.
  const int type = CCTK_GroupTypeI(group_index);
  const int rl = type == CCTK_GF ? reflevel : 0;

  int const active_tl =
    groupdata.AT(group_index).activetimelevels.AT(mglevel).AT(rl);
  const int tl = active_tl > 1 ? timelevel : 0;

  char const *const where = "ApplyPhysicalBCsForVarI";
  static Timers::Timer timer(where);
  timer.start();

  bool any_driver_bc = false;
  int const vstart = CCTK_FirstVarIndexI(group_index);
  int const vnum   = CCTK_NumVarsInGroupI(group_index);
  for(int var = 0; var < vnum; var++) {
    int const var_index = vstart + var;

    if(not boundary_conditions.count(var_index)) {
#ifdef PRESYNC_DEBUG
      std::cout << "ApplyBC: No boundary_conditions for " << CCTK_GroupName(var_index) << std::endl;
#endif
      continue;
    }

    int const map0 = 0;
    ggf *const ff = arrdata.AT(group_index).AT(map0).data.AT(var);
    assert(ff);

    // TODO: keep track of which faces are valid?
    assert(is_set(ff->valid(mglevel, rl, tl),
                  CCTK_VALID_INTERIOR | CCTK_VALID_GHOSTS));

    const Bounds& bounds = boundary_conditions[var_index];
    for(auto b: bounds.bounds) {
      int const ierr = Boundary_SelectVarForBCI(cctkGH, b.faces, b.width,
        b.table_handle, var_index, b.bc_name.c_str());
      if(ierr < 0) {
        CCTK_VWARN(CCTK_WARN_ALERT, "Failed to select boundary condition '%s' for variable '%s': %d",
        b.bc_name.c_str(), CCTK_FullVarName(var_index), ierr);
      }
      any_driver_bc = true;
    }
  }
  if(any_driver_bc) {
    do_applyphysicalbcs = true;
    int const ierr = ScheduleTraverse(where, "Driver_ApplyBCs",
                                      const_cast<cGH*>(cctkGH));
    if(ierr)
      CCTK_VWARN(CCTK_WARN_ALERT, "Failed ot traverse group Driver_ApplyBCs: %d\n", ierr);
    do_applyphysicalbcs = false;

    // sanity check on thorn Boundary
    const int type = CCTK_GroupTypeI(group_index);
    const int rl = type == CCTK_GF ? reflevel : 0;

    bool any_driver_bc_failed = false;
    for(int var = 0; var < vnum; var++) {
      int const var_index = vstart + var;

      int const map0 = 0;
      ggf *const ff = arrdata.AT(group_index).AT(map0).data.AT(var);
      assert(ff);

      if(boundary_conditions.count(var_index) and
         not is_set(ff->valid(mglevel, rl, tl), CCTK_VALID_EVERYWHERE)) {
        CCTK_VWARN(CCTK_WARN_ALERT,
                   "Internal error: thorn Boundary did not mark boundary of '%s' as valid when applying boundary conditions",
                   CCTK_FullVarName(var_index));
        any_driver_bc_failed = true;
      }
    }
    if(any_driver_bc_failed) {
        CCTK_VERROR("Internal error: thorn Boundary did not mark boundary of '%s' as valid when applying boundary conditions",
                   CCTK_GroupName(group_index));
    }
  }

  timer.stop();
}

} // namespace Carpet
