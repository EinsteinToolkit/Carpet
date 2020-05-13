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
#include "PreSyncCarpet.hh"
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

std::string format_var_tuple(const int vi, const int rl, const int tl) {
  std::ostringstream o;
  o << CCTK_FullVarName(vi);
  for(int i=0;i<tl;i++)
    o << "_p";
  o << " (rl=" << rl << ")";
  return o.str();
}

inline bool is_set(int flags,int flag) {
  return (flags & flag) == flag;
}

/**
 * Provide a string representation
 * of the code for the boundary.
 */
inline std::string wstr(int wh) {
  std::string s;
  if(wh == WH_EVERYWHERE) {
    s = "Everywhere";
    return s;
  }
  if(wh == 0) {
    s = "Nowhere";
    return s;
  }
  if(is_set(wh,WH_INTERIOR))
    s += "Interior";
  if(is_set(wh,WH_BOUNDARY)) {
    if(s.size() > 0) s+= "With";
    s += "Boundary";
  }
  if(is_set(wh,WH_GHOSTS)) {
    if(s.size() > 0) s+= "With";
    s += "Ghosts";
  }
  return s;
}

ostream& dumpValid(ostream& os, const int vi) {
  assert(vi < CCTK_NumVars());
  int const gi = CCTK_GroupIndexFromVarI(vi);
  assert(gi >= 0);
  int const var = vi - CCTK_FirstVarIndexI(gi);

  int const m = 0; // FIXME: this assumes that validity is the same on all maps

  ggf *const ff = arrdata.AT(gi).AT(m).data.AT(var);
  assert(ff);

  os << "\nValid entries:";
  for (int rl=0; rl<ff->h.reflevels(); ++rl) {
    for (int tl=0; tl<ff->timelevels(mglevel, rl); ++tl) {
      os << " " << format_var_tuple(vi, rl, tl);
      os << " " << wstr(ff->valid(mglevel, rl, tl));
    } // tl
  } // rl

  return os;
}

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

} // namespace

void PostCheckValid(cFunctionData *attribute, cGH *cctkGH, vector<int> const &sync_groups) {
  DECLARE_CCTK_PARAMETERS;

  if (not use_psync)
    return;

  for (int i=0;i<attribute->n_RDWR;i++) {
    const RDWR_entry& entry = attribute->RDWR[i];

    assert(entry.var_id < CCTK_NumVars());
    int const gi = CCTK_GroupIndexFromVarI(entry.var_id);
    assert(gi >= 0);
    int const var = entry.var_id - CCTK_FirstVarIndexI(gi);
    int const m = 0; // FIXME: this assumes that validity is the same on all maps
    ggf *const ff = arrdata.AT(gi).AT(m).data.AT(var);
    assert(ff);

    int type = CCTK_GroupTypeFromVarI(entry.var_id);
    int const rl = type == CCTK_GF ? reflevel : 0;

    // skip checking variables without storage
    if(not(entry.time_level < ff->timelevels(mglevel, rl)))
      continue;

    if(entry.where_wr == WH_INTERIOR) {
      ff->set_valid(mglevel, rl, entry.time_level, WH_INTERIOR);
    } else {
      int const old_valid = ff->valid(mglevel, rl, entry.time_level);
      ff->set_valid(mglevel, rl, entry.time_level, old_valid | entry.where_wr);
    }
  }
}

/**
 * Do the actual presync of the groups.
 **/
void PreSyncGroups(cFunctionData *attribute,cGH *cctkGH,const std::set<int>& pregroups) {
  DECLARE_CCTK_PARAMETERS;
  std::vector<int> sync_groups;

  if(use_psync and reflevel > 0) {
    // recurse to check that all coarsers levels are properly SYNCed
    CCTK_REAL previous_time = cctkGH->cctk_time;
    const int parent_reflevel = reflevel - 1;
    BEGIN_GLOBAL_MODE(cctkGH) {
      ENTER_LEVEL_MODE(cctkGH, parent_reflevel) {
        cctkGH->cctk_time = tt->get_time(mglevel, reflevel, timelevel);
        PreSyncGroups(attribute, cctkGH, pregroups);
      } LEAVE_LEVEL_MODE;
    } END_GLOBAL_MODE;
    cctkGH->cctk_time = previous_time;
  }

  for(auto i=pregroups.begin();i != pregroups.end();++i) {
    int gi = *i;
    int i0 = CCTK_FirstVarIndexI(gi);
    int iN = i0+CCTK_NumVarsInGroupI(gi);
    bool push = true;
    if(not use_psync && psync_error) {
      for(int vi=i0;vi<iN;vi++) {
        int const m = 0; // FIXME: this assumes that validity is the same on all maps
        ggf *const ff = arrdata.AT(gi).AT(m).data.AT(vi - i0);
        assert(ff);
        int type = CCTK_GroupTypeFromVarI(vi);
        int const rl = type == CCTK_GF ? reflevel : 0;
        int const wh = ff->valid(mglevel, rl, 0);
        if(not is_set(wh,WH_GHOSTS)) {
          std::ostringstream msg;
          msg << "  Presync: missing sync for ";
          msg << attribute->thorn << "::" << attribute->routine;
          msg << " in/at " << attribute->where << " variable " << CCTK_FullVarName(vi);
          CCTK_WARN(0,msg.str().c_str());
        }
        if(not is_set(wh,WH_BOUNDARY)) {
          std::ostringstream msg;
          msg << "  Presync: BC not valid for ";
          msg << attribute->thorn << "::" << attribute->routine;
          msg << " in/at " << attribute->where << " variable " << CCTK_FullVarName(vi);
          CCTK_WARN(0,msg.str().c_str());
        }
      }
    }
    for(int vi=i0;vi<iN;vi++) {
      int type = CCTK_GroupTypeFromVarI(vi);
      if(type == CCTK_GF && CCTK_VarTypeSize(CCTK_VarTypeI(vi)) == sizeof(CCTK_REAL)) {
        int const m = 0; // FIXME: this assumes that validity is the same on all maps
        ggf *const ff = arrdata.AT(gi).AT(m).data.AT(vi - i0);
        assert(ff);
        int wh = ff->valid(mglevel, reflevel, 0);
        if(is_set(wh,WH_BOUNDARY|WH_GHOSTS)) {
          continue;
        }
        if(!is_set(wh,WH_INTERIOR)) {
          std::ostringstream msg;
          msg << "SYNC of variable with invalid interior. Name: "
            << CCTK_FullVarName(vi) << " in routine "
            << attribute->thorn << "::" << attribute->routine;
          int level = psync_error ? 0 : 1;
          static bool have_warned = false;
          if(not have_warned) {
            CCTK_WARN(level,msg.str().c_str());
            have_warned = true;
          }
        }
        if(push) {
          sync_groups.push_back(gi);
          push = false;
        }
      }
    }
  }
  if(sync_groups.size()>0) {
    if(use_psync) {
      SyncProlongateGroups(cctkGH, sync_groups, attribute);
      for (size_t sgi=0;sgi<sync_groups.size();sgi++) {
        int i0 = CCTK_FirstVarIndexI(sync_groups[sgi]);
        int iN = i0+CCTK_NumVarsInGroupI(sync_groups[sgi]);
        for (int vi=i0;vi<iN;vi++) {
          int const m = 0; // FIXME: this assumes that validity is the same on all maps
          ggf *const ff = arrdata.AT(sync_groups[sgi]).AT(m).data.AT(vi - i0);
          assert(ff);
          int type = CCTK_GroupTypeFromVarI(vi);
          int const rl = type == CCTK_GF ? reflevel : 0;
          assert(rl >= 0);
          if(not is_set(ff->valid(mglevel, rl, 0), WH_EVERYWHERE)) {
            std::ostringstream msg;
            msg << "Required: Valid Everywhere, Observed: Valid " << wstr(ff->valid(mglevel, rl, 0)) << " " << CCTK_FullVarName(vi);
            msg << " Routine: " << attribute->thorn << "::" << attribute->routine;
            dumpValid(msg,vi);
            msg << std::endl;
            CCTK_WARN(0,msg.str().c_str());
          }
        }
      }
    }
  }
}

/**
 * Computes which groups need to be presync'd.
 */
void PreCheckValid(cFunctionData *attribute,cGH *cctkGH,std::set<int>& pregroups) {
  DECLARE_CCTK_PARAMETERS;
  if(cctkGH == 0) return;
  if(attribute == 0) return;
  if(not use_psync) return;

  for(int i=0;i<attribute->n_RDWR;i++) {
    const RDWR_entry& entry = attribute->RDWR[i];
    if(entry.where_rd == WH_NOWHERE) { // only READS cn trigger sync or errors
      continue;
    }
    // clauses to not refer to reflevel but the valid states have them so
    // inject them here
    const int vi = entry.var_id;
    if(CCTK_GroupTypeFromVarI(vi) != CCTK_GF) {
      continue;
    }
    // some thorns only read from variables based on extra parameters, eg
    // TmunuBase's stress_energy_state grid scalar. This check mimics the
    // behaviour of calling SYNC on a variable without storage by outputting a
    // high level warning. A thorn that actually tries to read from the
    // variable will encounter a SEGFAULT.
    if(CCTK_ActiveTimeLevelsVI(cctkGH, vi) == 0) {
      CCTK_VWARN(CCTK_WARN_DEBUG,
                 "Declared access to '%s' which has no storage",
                 CCTK_FullVarName(vi));
      continue;
    }

    assert(entry.var_id < CCTK_NumVars());
    int const gi = CCTK_GroupIndexFromVarI(entry.var_id);
    assert(gi >= 0);
    int const var = entry.var_id - CCTK_FirstVarIndexI(gi);
    int const m = 0; // FIXME: this assumes that validity is the same on all maps
    ggf *const ff = arrdata.AT(gi).AT(m).data.AT(var);
    assert(ff);

    int type = CCTK_GroupTypeFromVarI(vi);
    int const rl = type == CCTK_GF ? reflevel : 0;

    // skip checking variables without storage
    if(not(entry.time_level < ff->timelevels(mglevel, rl)))
      continue;

    int const valid = ff->valid(mglevel, rl, entry.time_level);

    // TODO: only need to check that what is READ is valid
    if(not is_set(valid,WH_INTERIOR))
    {
      // If the read spec is everywhere and we only have
      // interior, that's ok. The system will sync.
      std::ostringstream msg; 
      msg << "Required read for " << format_var_tuple(vi, rl, entry.time_level)
          << " not satisfied. Invalid interior"
          << " at the start of routine "
          << attribute->thorn << "::" << attribute->routine;
      dumpValid(msg, vi);
      int level = psync_error ? 0 : 1;
      static bool have_warned = false;
      if(not have_warned) {
        CCTK_WARN(level,msg.str().c_str()); 
        have_warned = true;
      }
    } else if(entry.time_level > 0) {
      if(not is_set(valid,entry.where_rd)) {
        // If the read spec is everywhere, that's not
        // okay for previous time levels as they won't sync.
        std::ostringstream msg; 
        msg << "Required read for " << format_var_tuple(vi, rl, entry.time_level)
          << " not satisfied. Wanted '" << wstr(entry.where_rd)
          << "' found '" << wstr(valid)
          << "' at the start of routine "
          << attribute->thorn << "::" << attribute->routine;
        dumpValid(msg, vi);
        CCTK_Error(__LINE__,__FILE__,CCTK_THORNSTRING,msg.str().c_str());
      }
    }
    if(entry.where_rd == WH_EVERYWHERE) {
      if(is_set(valid,WH_INTERIOR) && not is_set(valid, WH_EVERYWHERE))
      {
        if(entry.time_level != 0) {
          std::ostringstream msg;
          msg << "Attempt to sync previous time level (tl=" << entry.time_level << ")"
              << " for " << CCTK_FullVarName(vi);
          CCTK_Error(__LINE__,__FILE__,CCTK_THORNSTRING,msg.str().c_str());
        }

        int g = CCTK_GroupIndexFromVarI(vi);
        pregroups.insert(g);
      } else if(not is_set(valid,WH_INTERIOR)) {
        std::ostringstream msg; 
        msg << "Cannot sync " << CCTK_FullVarName(vi)
            << " because it is not valid in the interior.";
        dumpValid(msg, vi);
        int level = psync_error ? CCTK_WARN_ABORT : CCTK_WARN_ALERT;
        static bool have_warned = false;
        if(not have_warned) {
          CCTK_WARN(level,msg.str().c_str()); 
          have_warned = true;
        }
      }
    }
  }
}

/**
 * Given a variable and a timelevel, set the region
 * of the grid where that variable is valid (i.e. the where_spec).
 */
extern "C"
void Carpet_SetValidRegion(CCTK_POINTER_TO_CONST /*cctkGH_*/, CCTK_INT vi,
                           CCTK_INT tl, CCTK_INT wh) {
  if(vi < 0 or vi >= CCTK_NumVars()) {
    CCTK_VERROR("Invalid variable index %d", vi);
  }

  int const gi = CCTK_GroupIndexFromVarI(vi);
  assert(gi >= 0);
  int const var = vi - CCTK_FirstVarIndexI(gi);
  int const m = 0; // FIXME: this assumes that validity is the same on all maps
  ggf *const ff = arrdata.AT(gi).AT(m).data.AT(var);
  assert(ff);
  int type = CCTK_GroupTypeFromVarI(vi);
  int const rl = type == CCTK_GF ? reflevel : 0;

  if(tl >= ff->timelevels(mglevel, rl)) {
    CCTK_VERROR("Invalid time level %d for variable %s.", tl,
                CCTK_FullVarName(vi));
  }

  ff->set_valid(mglevel, rl, tl, wh);
}

/**
 * Given a variable and a timelevel, return the region
 * of the grid where that variable is valid (i.e. the where_spec).
 */
extern "C"
CCTK_INT Carpet_GetValidRegion(CCTK_POINTER_TO_CONST /*cctkGH_*/, CCTK_INT vi,
                               CCTK_INT tl) {
  if(vi < 0 or vi >= CCTK_NumVars()) {
    CCTK_VERROR("Invalid variable index %d", vi);
  }
  int const gi = CCTK_GroupIndexFromVarI(vi);
  assert(gi >= 0);
  int const var = vi - CCTK_FirstVarIndexI(gi);
  int const m = 0; // FIXME: this assumes that validity is the same on all maps
  ggf *const ff = arrdata.AT(gi).AT(m).data.AT(var);
  assert(ff);
  int type = CCTK_GroupTypeFromVarI(vi);
  int const rl = type == CCTK_GF ? reflevel : 0;

  if(tl >= ff->timelevels(mglevel, rl)) {
    CCTK_VERROR("Invalid time level %d for variable %s.", tl,
                CCTK_FullVarName(vi));
  }

  return ff->valid(mglevel, rl, tl);
}

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

  assert(variables or nvariables == 0);
  assert(tls or nvariables == 0);

  assert(is_level_mode());

  if(use_psync and reflevel > 0) {
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

  std::set<std::tuple<int,int>> processed_groups; // [gi] [tl]
  for(int i = 0; i < nvariables; i++) {
    int const vi = variables[i];
    int const tl = tls[i];
    int const where = wheres[i];

    // Do nothing if this is not a grid function
    if(CCTK_GroupTypeFromVarI(vi) != CCTK_GF) {
      continue;
    }

    if(vi < 0 or vi >= CCTK_NumVars()) {
      CCTK_VERROR("Invalid variable index %d", vi);
    }
    int const gi = CCTK_GroupIndexFromVarI(vi);
    assert(gi >= 0);
    int const var = vi - CCTK_FirstVarIndexI(gi);
    int const m = 0; // FIXME: this assumes that validity is the same on all maps
    ggf *const ff = arrdata.AT(gi).AT(m).data.AT(var);
    assert(ff);
    if(tl >= ff->timelevels(mglevel, reflevel)) {
      CCTK_VERROR("Invalid time level %d for variable %s.", tl,
                  CCTK_FullVarName(vi));
    }
    int const valid = ff->valid(mglevel, reflevel, tl);

    // Check if anything needs to be done
    if(is_set(valid, where)) {
      continue;
    }

    if(not is_set(valid, WH_INTERIOR)) {
      static std::set<int> have_warned_about;
      if(not have_warned_about.count(vi)) {
        CCTK_VWARN(psync_error ? CCTK_WARN_ABORT : CCTK_WARN_ALERT,
                   "SYNC requires valid data in interior %s rl=%d tl=%d but have only %s",
                   CCTK_FullVarName(vi), reflevel, tl, wstr(valid).c_str());
        have_warned_about.insert(vi);
      }
      continue; // cannot SYNC, do not update valid state
    }

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
    const_cast<cGH*>(cctkGH)->cctk_time = tt->get_time(mglevel, reflevel, timelevel);
    int ierr = SyncProlongateGroups(cctkGH, sync_groups, attribute);
    assert(not ierr);
    const_cast<cGH*>(cctkGH)->cctk_time = old_cctk_time;
    timelevel_offset = old_timelevel_offset;
    timelevel = old_timelevel;
    do_allow_past_timelevels = old_do_allow_past_timelevels;
  }

  return 0;
}
}

extern "C"
CCTK_INT Carpet_RequireValidData(CCTK_POINTER_TO_CONST cctkGH_,
   CCTK_INT const * const variables, CCTK_INT const * const tls,
   CCTK_INT const nvariables, CCTK_INT const * wheres) {

  const cGH *cctkGH = static_cast<const cGH*>(cctkGH_);

  if(not is_level_mode() and not is_global_mode() and not is_meta_mode())
    CCTK_VERROR("%s must be called in level, global, or meta mode", __func__);

  // Take action by mode
  if(is_level_mode()) {
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
    CCTK_BUILTIN_UNREACHABLE();
  }

  return 0;
}

extern "C"
CCTK_INT Carpet_NotifyDataModified(CCTK_POINTER_TO_CONST /*cctkGH_*/,
   CCTK_INT const * const variables, CCTK_INT const * const tls,
   CCTK_INT const nvariables, CCTK_INT const * wheres) {
  DECLARE_CCTK_PARAMETERS;

  // technically LEVEL mode is the only allowed one since I keep track of
  // validity on a level granularity. However that makes the routine awkward to
  // use in locally scheduled routines.
  if(not is_level_mode() and not is_local_mode())
    CCTK_VERROR("%s must be called in level, of local mode", __func__);

  assert(variables or nvariables == 0);
  assert(tls or nvariables == 0);

  for(int i = 0; i < nvariables; i++) {
    int const vi = variables[i];
    int const tl = tls[i];
    int const where = wheres[i];

    // Do nothing if this is not a grid function
    if(CCTK_GroupTypeFromVarI(vi) != CCTK_GF) {
      continue;
    }

    if(vi < 0 or vi >= CCTK_NumVars()) {
      CCTK_VERROR("Invalid variable index %d", vi);
    }

    int const gi = CCTK_GroupIndexFromVarI(vi);
    assert(gi >= 0);
    int const var = vi - CCTK_FirstVarIndexI(gi);
    int const m = 0; // FIXME: this assumes that validity is the same on all maps
    ggf *const ff = arrdata.AT(gi).AT(m).data.AT(var);
    assert(ff);
    if(tl >= ff->timelevels(mglevel, reflevel)) {
      CCTK_VERROR("Invalid time level %d for variable %s.", tl,
                  CCTK_FullVarName(vi));
    }
    int const valid = ff->valid(mglevel, reflevel, tl);
    ff->set_valid(mglevel, reflevel, tl, valid | where);
  }

  return 0;
}

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
  const CCTK_INT vstart = CCTK_FirstVarIndexI(group);
  const CCTK_INT vnum   = CCTK_NumVarsInGroupI(group);
  CCTK_INT ierr = 0;
  for(CCTK_INT vi=vstart;vi<vstart+vnum;vi++) {
    const CCTK_INT myierr =
      SelectVarForBCI(cctkGH,faces,width,table_handle,vi,bc_name);
    if(ierr and not myierr)
      ierr = myierr;
  }
  return ierr;
}

extern "C"
CCTK_INT Carpet_FilterOutVarForBCI(
    const CCTK_POINTER_TO_CONST cctkGH_,
    const CCTK_INT var_index) {
  // do apply physical BC if we are being called from Carpet's own Driver BC
  // routine
  if(do_applyphysicalbcs)
    return false;
  auto it = boundary_conditions.find(var_index);
  return it != boundary_conditions.end();
}

/**
 * Apply boundary conditions for a group. Called from inside SyncProlongateGroups.
 */
void ApplyPhysicalBCsForGroupI(const cGH *cctkGH, const int group_index) {
  DECLARE_CCTK_PARAMETERS;

  if(not use_psync)
    return;

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

    int const m = 0; // FIXME: this assumes that validity is the same on all maps
    ggf *const ff = arrdata.AT(group_index).AT(m).data.AT(var);
    assert(ff);
    int const type = CCTK_GroupTypeI(group_index);
    int const rl = type == CCTK_GF ? reflevel : 0;

    // TODO: keep track of which faces are valid?
    assert(is_set(ff->valid(mglevel, rl, 0), WH_INTERIOR | WH_GHOSTS));

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
    bool any_driver_bc_failed = false;
    for(int var = 0; var < vnum; var++) {
      int const var_index = vstart + var;

      int const m = 0; // FIXME: this assumes that validity is the same on all maps
      ggf *const ff = arrdata.AT(group_index).AT(m).data.AT(var);
      assert(ff);
      int const type = CCTK_GroupTypeI(group_index);
      int const rl = type == CCTK_GF ? reflevel : 0;

      if(boundary_conditions.count(var_index) and
         not is_set(ff->valid(mglevel, rl, 0), WH_EVERYWHERE)) {
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
