#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cmath>
#include <cassert>
#include <cstring>

#include <iostream>
#include <locale>
#include <map>
#include <set>
#include <vector>
#include <sstream>

#include "variables.hh"
#include "modes.hh"
#include "carpet.hh"
#include "PreSyncCarpet.hh"
#include "Carpet_Prototypes.h"

#undef PRESYNC_DEBUG

namespace Carpet {

// local functions and types
namespace {

// boundary and symmstry condition handling
typedef CCTK_INT (*boundary_function)(
  const cGH *cctkGH,
  int num_vars,
  const int *var_indices,
  const int *faces,
  const int *widths,
  const int *table_handles);

typedef CCTK_INT (*iface_boundary_function)(
  CCTK_POINTER_TO_CONST cctkGH,
  CCTK_INT num_vars,
  const CCTK_INT *var_indices,
  const CCTK_INT *faces,
  const CCTK_INT *boundary_widths,
  const CCTK_INT *table_handles);

typedef CCTK_INT (*sym_boundary_function)(
  const cGH *cctkGH, const CCTK_INT var_index);

typedef CCTK_INT (*sym_iface_boundary_function)(
  CCTK_POINTER_TO_CONST cctkGH, const CCTK_INT var_index);

struct Bound {
  std::string bc_name;
  CCTK_INT faces;
  CCTK_INT width;
  CCTK_INT table_handle;
};

struct Func {
  boundary_function func;
};

struct SymFunc {
  sym_boundary_function func;
  SymFunc() {
    func = nullptr;
  }
};

// these keep track of registered boundary and symmetry routines as well as
// selected boundary condistions
std::map<std::string,Func> boundary_functions;
std::map<std::string,SymFunc> symmetry_functions;
std::map<int,std::vector<Bound>> boundary_conditions;

std::string format_var_tuple(const int vi, const int rl, const int tl) {
  std::ostringstream o;
  o << CCTK_FullVarName(vi);
  for(int i=0;i<tl;i++)
    o << "_p";
  if(rl != -1) o << " (rl=" << rl << ")";
  return o.str();
}

inline bool on(int flags,int flag) {
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
  if(on(wh,WH_INTERIOR))
    s += "Interior";
  if(on(wh,WH_BOUNDARY)) {
    if(s.size() > 0) s+= "With";
    s += "Boundary";
  }
  if(on(wh,WH_GHOSTS)) {
    if(s.size() > 0) s+= "With";
    s += "Ghosts";
  }
  return s;
}

inline void tolower(std::string& s) {
  for(auto si=s.begin();si != s.end();++si) {
    *si = std::tolower(*si);
  }
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

/**
 * Called by ManualSyncGF.
 */
void Sync1(const cGH *cctkGH,int tl,int gi) {
  std::vector<int> sync_groups;
  sync_groups.push_back(gi);
  cFunctionData *attribute = 0;
  // copied out of modes.hh
  int const old_timelevel_offset = timelevel_offset;
  int const old_timelevel = timelevel;
  CCTK_REAL old_cctk_time = cctkGH->cctk_time;
  auto const old_do_allow_past_timelevels = do_allow_past_timelevels;
  do_allow_past_timelevels = tl == 0;
  timelevel_offset = timelevel = tl;
  int type = CCTK_GroupTypeI(gi);
  int const rl = type == CCTK_GF ? reflevel : 0;
  const_cast<cGH*>(cctkGH)->cctk_time = tt->get_time(mglevel, rl, timelevel);
  int ierr = SyncProlongateGroups(cctkGH, sync_groups, attribute);
  assert(!ierr);
  const_cast<cGH*>(cctkGH)->cctk_time = old_cctk_time;
  timelevel_offset = old_timelevel_offset;
  timelevel = old_timelevel;
  do_allow_past_timelevels = old_do_allow_past_timelevels;
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
    if(!use_psync && psync_error) {
      for(int vi=i0;vi<iN;vi++) {
        int const m = 0; // FIXME: this assumes that validity is the same on all maps
        ggf *const ff = arrdata.AT(gi).AT(m).data.AT(vi - i0);
        assert(ff);
        int type = CCTK_GroupTypeFromVarI(vi);
        int const rl = type == CCTK_GF ? reflevel : 0;
        int const wh = ff->valid(mglevel, rl, 0);
        if(!on(wh,WH_GHOSTS)) {
          std::ostringstream msg;
          msg << "  Presync: missing sync for ";
          msg << attribute->thorn << "::" << attribute->routine;
          msg << " in/at " << attribute->where << " variable " << CCTK_FullVarName(vi);
          CCTK_WARN(0,msg.str().c_str());
        }
        if(!on(wh,WH_BOUNDARY)) {
          std::ostringstream msg;
          msg << "  Presync: BC not valid for ";
          msg << attribute->thorn << "::" << attribute->routine;
          msg << " in/at " << attribute->where << " variable " << CCTK_FullVarName(vi);
          CCTK_WARN(0,msg.str().c_str());
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
          if(ff->valid(mglevel, reflevel, 0) != WH_EVERYWHERE) {
            std::ostringstream msg;
            msg << "Required: Valid Everywhere, Observed: Valid " << wstr(ff->valid(mglevel, reflevel, 0)) << " " << CCTK_FullVarName(vi);
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
  if(!use_psync) return;

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
    if(!on(valid,WH_INTERIOR)) // and !silent_psync)
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
      if(!on(valid,entry.where_rd)) {
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
      if(on(valid,WH_INTERIOR) && valid != WH_EVERYWHERE)
      {
        if(entry.time_level != 0) {
          std::ostringstream msg;
          msg << "Attempt to sync previous time level (tl=" << entry.time_level << ")"
              << " for " << CCTK_FullVarName(vi);
          CCTK_Error(__LINE__,__FILE__,CCTK_THORNSTRING,msg.str().c_str());
        }

        int g = CCTK_GroupIndexFromVarI(vi);
        pregroups.insert(g);
      } else if(!on(valid,WH_INTERIOR))
      {
        std::ostringstream msg; 
        msg << "Cannot sync " << CCTK_FullVarName(vi)
            << " because it is not valid in the interior.";
        dumpValid(msg, vi);
        int level = psync_error ? 0 : 1;
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
// TODO: expand to take a reflevel argument?
void SetValidRegion(int vi,int tl,int wh) {
  assert(vi < CCTK_NumVars());
  int const gi = CCTK_GroupIndexFromVarI(vi);
  assert(gi >= 0);
  int const var = vi - CCTK_FirstVarIndexI(gi);
  int const m = 0; // FIXME: this assumes that validity is the same on all maps
  ggf *const ff = arrdata.AT(gi).AT(m).data.AT(var);
  assert(ff);
  int type = CCTK_GroupTypeFromVarI(vi);
  int const rl = type == CCTK_GF ? reflevel : 0;

  ff->set_valid(mglevel, rl, tl, wh);
}

/**
 * Given a variable and a timelevel, return the region
 * of the grid where that variable is valid (i.e. the where_spec).
 */
// TODO: expand to take a reflevel argument?
int GetValidRegion(int vi,int tl) {
  assert(vi < CCTK_NumVars());
  int const gi = CCTK_GroupIndexFromVarI(vi);
  assert(gi >= 0);
  int const var = vi - CCTK_FirstVarIndexI(gi);
  int const m = 0; // FIXME: this assumes that validity is the same on all maps
  ggf *const ff = arrdata.AT(gi).AT(m).data.AT(var);
  assert(ff);
  int type = CCTK_GroupTypeFromVarI(vi);
  int const rl = type == CCTK_GF ? reflevel : 0;

  return ff->valid(mglevel, rl, tl);
}

/**
 * Apply manual synchronization to variable vi. If this variable
 * is already valid everywhere, this routine does nothing. When
 * the routine finishes, it will be valid everywhere.
 */
extern "C" void Carpet_ManualSyncGF(CCTK_POINTER_TO_CONST cctkGH_,const CCTK_INT tl,const CCTK_INT vi) {
  // Do nothing if this is not a GF
  if(CCTK_GroupTypeFromVarI(vi) != CCTK_GF) {
    return;
  }
  const cGH *cctkGH = static_cast<const cGH*>(cctkGH_);
  assert(vi < CCTK_NumVars());
  int const gi = CCTK_GroupIndexFromVarI(vi);
  assert(gi >= 0);
  int const var = vi - CCTK_FirstVarIndexI(gi);
  int const m = 0; // FIXME: this assumes that validity is the same on all maps
  ggf *const ff = arrdata.AT(gi).AT(m).data.AT(var);
  assert(ff);
  int type = CCTK_GroupTypeFromVarI(vi);
  int const rl = type == CCTK_GF ? reflevel : 0;
  int const valid = ff->valid(mglevel, rl, tl);
  // Check if anything needs to be done
  if(valid == WH_EVERYWHERE) {
    return;
  }
  if(on(valid,WH_INTERIOR)) {
    dumpValid(std::cerr, vi) << std::endl;
    CCTK_VERROR("SYNC requires valid data in interior %s rl=%d tl=%d", CCTK_FullVarName(vi), rl, tl);
  }
  assert(on(valid,WH_INTERIOR));

  // Take action by mode
  if(is_level_mode()) {
    Sync1(cctkGH,tl,gi);
  } else if(is_global_mode()) {
    BEGIN_REFLEVEL_LOOP(cctkGH) {
      Sync1(cctkGH,tl,gi);
    }
    END_REFLEVEL_LOOP;
  } else if(is_meta_mode()) {
    BEGIN_MGLEVEL_LOOP(cctkGH) {
      BEGIN_REFLEVEL_LOOP(cctkGH) {
        Sync1(cctkGH,tl,gi);
      }
      END_REFLEVEL_LOOP;
    }
    END_MGLEVEL_LOOP;
  } else {
    abort();
  }
}

extern "C"
void Carpet_SynchronizationRecovery(CCTK_ARGUMENTS) {
  for(int vi = 0; vi < CCTK_NumVars(); vi++) {
    SetValidRegion(vi,0,WH_INTERIOR);
    int tl = CCTK_ActiveTimeLevelsVI(cctkGH, vi);
    for(int time = 1; time < tl; time++) {
      SetValidRegion(vi,time,WH_EVERYWHERE);
    }
  }
}

extern "C"
CCTK_INT Carpet_RegisterPhysicalBC(
    CCTK_POINTER_TO_CONST /*cctkGH_*/,
    iface_boundary_function func_,
    CCTK_STRING bc_name) {
  boundary_function func = reinterpret_cast<boundary_function>(func_);
  if(NULL==func) {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,  
               "Physical Boundary condition '%s' points to NULL.", bc_name);
  }
  assert(bc_name != nullptr);
  std::string name{bc_name};
  tolower(name);
  Func& f = boundary_functions[name];
  f.func = func;
  return 0;
}

extern "C"
CCTK_INT Carpet_RegisterSymmetryBC(
    CCTK_POINTER_TO_CONST /*cctkGH_*/,
    sym_iface_boundary_function func_,
    CCTK_STRING bc_name) {
  sym_boundary_function func = reinterpret_cast<sym_boundary_function>(func_);
  if(NULL==func) {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,  
               "Symmetry Boundary condition '%s' points to NULL.", bc_name);
  }
  SymFunc& f = symmetry_functions[bc_name];
  f.func = func;
  return 0;
}

namespace {
CCTK_INT SelectVarForBCI(
    const cGH *cctkGH,
    const CCTK_INT faces,
    const CCTK_INT width,
    const CCTK_INT table_handle,
    const CCTK_INT var_index,
    const CCTK_STRING bc_name) {
  std::string name{bc_name};
  tolower(name);
  if(!boundary_functions.count(name.c_str())) {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,  
               "Requested BC '%s' not found.", bc_name);
  }
  assert(var_index != 0);
  auto i = boundary_conditions.find(var_index);
  for(; i != boundary_conditions.end();++i) {
    for(auto b = i->second.begin(); b != i->second.end(); ++b) {
      if(b->bc_name == name) {
#ifdef PRESYNC_DEBUG
        std::cout << "Variable " << CCTK_FullVarName(var_index) << " has been selected twice for '" << bc_name << "'" << std::endl;
#endif
        return 0;
      }
    }
  }
  std::vector<Bound>& bv = boundary_conditions[var_index];
  Bound b;
  b.faces = faces;
  b.width = width;
  b.table_handle = table_handle;
  b.bc_name = name;
  bv.push_back(b);
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

/**
 * Apply boundary conditions for a single variable.
 */
void ApplyPhysicalBCsForVarI(const cGH *cctkGH, const int var_index) {
  DECLARE_CCTK_PARAMETERS;
  if(!use_psync) return;
  auto bc = boundary_conditions;
  if(bc.find(var_index) == bc.end()) {
#ifdef PRESYNC_DEBUG
    std::cout << "ApplyBC: No bc's for " << CCTK_FullVarName(var_index) << std::endl;
#endif
    return;
  }

  assert(var_index < CCTK_NumVars());
  int const gi = CCTK_GroupIndexFromVarI(var_index);
  assert(gi >= 0);
  int const var = var_index - CCTK_FirstVarIndexI(gi);
  int const m = 0; // FIXME: this assumes that validity is the same on all maps
  ggf *const ff = arrdata.AT(gi).AT(m).data.AT(var);
  assert(ff);
  int type = CCTK_GroupTypeFromVarI(var_index);
  int const rl = type == CCTK_GF ? reflevel : 0;

  std::vector<Bound>& bv = bc[var_index];
  BEGIN_LOCAL_MAP_LOOP(cctkGH, CCTK_GF) {
    BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) {
      for(auto j=bv.begin(); j != bv.end(); ++j) {
        Bound& b = *j;
        Func& f = boundary_functions.at(b.bc_name);

        if(use_psync) {
          int ierr = (*f.func)(cctkGH,1,&var_index,&b.faces,&b.width,&b.table_handle);
          assert(!ierr);
          for (auto iter = symmetry_functions.begin(); iter != symmetry_functions.end(); iter++) {
            std::string name = iter->first;
            SymFunc& fsym = symmetry_functions.at(name);
            ierr = (*fsym.func)(cctkGH, var_index);
          }
	  int const valid = ff->valid(mglevel, rl, 0);
	  ff->set_valid(mglevel, rl, 0, valid | WH_BOUNDARY);
        }
      }
    }
    END_LOCAL_COMPONENT_LOOP;
  }
  END_LOCAL_MAP_LOOP;
}

/**
 * Apply boundary conditions for a group. Called from inside SyncProlongateGroups.
 */
void ApplyPhysicalBCsForGroupI(const cGH *cctkGH, const int group_index) {
  int vstart = CCTK_FirstVarIndexI(group_index);
  int vnum   = CCTK_NumVarsInGroupI(group_index);
  for(int var_index=vstart;var_index<vstart+vnum;var_index++) {
    ApplyPhysicalBCsForVarI(cctkGH,var_index);
  }
}

} // namespace Carpet
