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
#include "PreSync.h"

#undef PRESYNC_DEBUG

extern "C" void ShowValid();

namespace Carpet {
int bnd_vi = -1;

struct var_tuple {
  const int vi; // var index
  const int rl; // refinement level
  const int tl; // time level;
  var_tuple() : vi(-1), rl(-1), tl(-1) {}
  var_tuple(int vi_,int rl_,int tl_) : vi(vi_), rl(CCTK_GroupTypeFromVarI(vi_) == CCTK_GF ? rl_ : -1), tl(tl_) {}
};

std::ostream& operator<<(std::ostream& o,const var_tuple& vt) {
  o << CCTK_FullVarName(vt.vi);
  for(int i=0;i<vt.tl;i++)
    o << "_p";
  if(vt.rl != -1) o << " (rl=" << vt.rl << ")";
  return o;
}

bool operator<(const var_tuple& v1,const var_tuple& v2) {
  int diff = v1.vi - v2.vi;
  if(diff == 0) diff = v1.rl - v2.rl;
  if(diff == 0) diff = v1.tl - v2.tl;
  return diff < 0;
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

std::string current_routine;

std::map<std::string,std::map<var_tuple,int>> reads,writes;
std::map<var_tuple,int> valid_k;
std::map<var_tuple,int> old_valid_k;

extern "C" void diagnosticPreValid();
void diagnosticPreValid() {
  old_valid_k = valid_k;
}

ostream& dumpValid(ostream& os, const int vi) {
  os << "\nValid entries:";
  for(auto it : valid_k) {
    if(vi == -1 || it.first.vi == vi) {
      os << " " << it.first << " " << wstr(valid_k[it.first]);
    }
  }
  return os;
}

extern "C" void diagnosticChanged();
void diagnosticChanged() {
  for(auto entry = valid_k.begin(); valid_k.end() != entry; ++entry) {
    auto f = old_valid_k.find(entry->first);
    if(f == old_valid_k.end()) {
#ifdef PRESYNC_DEBUG
      std::cout << " NEW: " << entry->first << " = " << wstr(entry->second) << std::endl;
#endif
    } else if(entry->second != f->second) {
#ifdef PRESYNC_DEBUG
      std::cout << " CHANGED: " << entry->first << " = " << wstr(f->second) << " -> " << wstr(entry->second) << std::endl;
#endif
    }
  }
}

extern "C" void TraverseReads(const char *func_name,void(*trace_func)(int,int,int)) {
    std::string f{func_name};
    auto r = reads[f];
    for(auto i = r.begin();i != r.end(); ++i) {
        trace_func(i->first.vi,i->first.tl,i->second);
    }
}

extern "C" void TraverseWrites(const char *func_name,void(*trace_func)(int,int,int)) {
    std::string f{func_name};
    auto r = writes[f];
    for(auto i = r.begin();i != r.end(); ++i) {
        trace_func(i->first.vi,i->first.tl,i->second);
    }
}

void PostCheckValid(cFunctionData *attribute, cGH *cctkGH, vector<int> const &sync_groups) {
  DECLARE_CCTK_PARAMETERS;

  if (not use_psync)
    return;

  std::string r;
  r += attribute->thorn;
  r += "::";
  r += attribute->routine;
  current_routine = "";
  std::map<var_tuple,int>& writes_m = writes[r];
  for(auto i = writes_m.begin();i != writes_m.end(); ++i) {
    var_tuple vi = i->first;
    var_tuple vt{vi.vi,reflevel,vi.tl};
    if(i->second == WH_INTERIOR)
      valid_k[vt] = WH_INTERIOR;
    else
      valid_k[vt] |= i->second;
  }


  for(auto g = sync_groups.begin();g != sync_groups.end();++g) {
    int gi = *g;
    int i0 = CCTK_FirstVarIndexI(gi);
    int iN = i0+CCTK_NumVarsInGroupI(gi);
    for(int vi=i0;vi<iN;vi++) {
      int& w = writes[r][var_tuple(vi,-1,0)];
      if(w == WH_INTERIOR) {
        var_tuple vt{vi,reflevel,0};
        valid_k[vt] |= WH_GHOSTS;
#ifdef PRESYNC_DEBUG
        std::cout << "SYNC: " << CCTK_FullVarName(vt.vi) << " " << wstr(valid_k[vt]) << std::endl;
#endif
      }
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
        var_tuple vt{vi,reflevel,0};
        int wh = valid_k[vt];
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
    for(int vi=i0;vi<iN;vi++) {
      int type = CCTK_GroupTypeFromVarI(vi);
      if(type == CCTK_GF && CCTK_VarTypeSize(CCTK_VarTypeI(vi)) == sizeof(CCTK_REAL)) {
        var_tuple vt{vi,reflevel,0};
        int wh = valid_k[vt];
        if(on(wh,WH_EXTERIOR)) {
          continue;
        }
        if(!on(wh,WH_INTERIOR)) {// and !silent_psync) {
          std::ostringstream msg;
          msg << "SYNC of variable with invalid interior. Name: "
            << CCTK_FullVarName(vi) << " before: " << current_routine;
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
        valid_k[vt] |= WH_GHOSTS;
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
          if(valid_k[var_tuple(vi,reflevel,0)] != WH_EVERYWHERE) {
            std::ostringstream msg;
            msg << "Required: Valid Everywhere, Observed: Valid " << wstr(valid_k[var_tuple(vi,reflevel,0)]) << " " << CCTK_FullVarName(vi);
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

std::map<var_tuple,int> tmp_read, tmp_write;

bool hasAccess(const std::map<var_tuple,int>& m, const var_tuple& vt) {
  auto i2 = m.find(vt);
  if(i2 == m.end())
    return false;
  return i2->second != WH_NOWHERE;
}

/**
 * Determine whether the currently executed function
 * has either read or write access to the variable
 * in question. If it does not, access is blocked
 * by causing CCTKi_VarDataPtrI() to return null.
 * Note that only REAL grid functions are affected.
 */
extern "C"
int Carpet_hasAccess(const cGH *cctkGH,int var_index) {
  DECLARE_CCTK_PARAMETERS;
  if(!psync_error)
    return true;
  int type = CCTK_GroupTypeFromVarI(var_index);
  if(type == CCTK_GF && CCTK_VarTypeSize(CCTK_VarTypeI(var_index)) == sizeof(CCTK_REAL)) {
    var_tuple vi{var_index,-1,0};
    if(hasAccess(reads[current_routine],vi))
      return true;
    if(hasAccess(writes[current_routine],vi))
      return true;
    if(hasAccess(tmp_read,vi)) {
#ifdef PRESYNC_DEBUG
      std::cout << "tmp acccess for " << CCTK_FullVarName(var_index) << std::endl;
#endif
      return true;
    }
    if(hasAccess(tmp_write,vi)) {
#ifdef PRESYNC_DEBUG
      std::cout << "tmp acccess for " << CCTK_FullVarName(var_index) << std::endl;
#endif
      return true;
    }
    return false;
  } else {
    return true;
  }
}

std::map<int,int> attempted_readwrites;

extern "C" void attempt_readwrite(const char *thorn,const char *var, int spec);
void attempt_readwrite(int gf,int spec);
extern "C" void clear_readwrites();
extern "C" void check_readwrites();

extern "C" void attempt_readwrite(const char *thorn,const char *var, int spec) {
  std::string v;
  v += thorn;
  v += "::";
  v += var;
  int vi = CCTK_VarIndex(v.c_str());
  if(vi < 0) abort();
  attempt_readwrite(vi,spec);
}
void attempt_readwrite(int gf,int spec) {
  if(gf >= 0)
    attempted_readwrites[gf] |= spec;
}
extern "C" void clear_readwrites() {
  attempted_readwrites.clear();
}
extern "C" void check_readwrites() {
  for(auto i=attempted_readwrites.begin(); i != attempted_readwrites.end(); ++i) {
    if((i->second & 0x01)==0x01 && !hasAccess(reads[current_routine],var_tuple(i->first,-1,0))) {
        std::cerr << "Undeclared access: " << current_routine << " read name='" << CCTK_FullVarName(i->first) << "'" << std::endl;
    }
    if((i->second & 0x10)==0x10 && !hasAccess(writes[current_routine],var_tuple(i->first,-1,0))) {
        std::cerr << "Undeclared access: " << current_routine << " write name='" << CCTK_FullVarName(i->first) << "'" << std::endl;
    }
  }
}

/**
 * Request temporary read access to a variable,
 * thus allowing Carpet_hasAccess() to return true.
 * This needs to be called <i>before</i> the 
 * CCTK_DECLARE_VARIABLES macro in order to have
 * the desired effect.
 */
extern "C"
void Carpet_requestAccess(int var_index,int read_spec,int write_spec) {
  assert(var_index >= 0);
  var_tuple vi{var_index,-1,0};
  tmp_read[vi]  |= read_spec;
  tmp_write[vi] |= write_spec;
}

void dump_clauses(std::map<var_tuple,int>& reads_m,std::map<var_tuple,int>& writes_m) {
  std::cout << "ROUTINE: " << current_routine << " {" << std::endl;
  for(auto i = reads_m.begin(); i != reads_m.end(); ++i) {
    std::cout << "  >>READ:  " << CCTK_FullVarName(i->first.vi) << "," << wstr(i->second) << std::endl;
  }
  for(auto i = writes_m.begin(); i != writes_m.end(); ++i) {
    std::cout << "  >>WRITE: " << CCTK_FullVarName(i->first.vi) << "," << wstr(i->second) << std::endl;
  }
  std::cout << "}" << std::endl;
}

/**
 * Computes which groups need to be presync'd.
 */
void PreCheckValid(cFunctionData *attribute,cGH *cctkGH,std::set<int>& pregroups) {
  DECLARE_CCTK_PARAMETERS;
  if(cctkGH == 0) return;
  if(attribute == 0) return;
  if(!use_psync) return;
  tmp_read.erase(tmp_read.begin(),tmp_read.end());
  tmp_write.erase(tmp_write.begin(),tmp_write.end());

  std::string r;
  r += attribute->thorn;
  r += "::";
  r += attribute->routine;
  current_routine = r;

  if(writes.find(r) == writes.end()) {

    std::map<var_tuple,int>& writes_m  = writes[r];
    std::map<var_tuple,int>& reads_m  = reads [r];

    for(int i=0;i < attribute->n_RDWR;i++) {
      RDWR_entry *rdwr = &attribute->RDWR[i];
      var_tuple vt{rdwr->var_id, -1, rdwr->time_level};
      if(rdwr->where_wr != 0) {
        writes_m[vt] = rdwr->where_wr;
      }
      if(rdwr->where_rd != 0) {
        reads_m[vt] = rdwr->where_rd;
      }
    }
  }

  std::map<var_tuple,int>& reads_m  = reads [r];

  for(auto i = reads_m.begin();i != reads_m.end(); ++i) {
    var_tuple vtemp = i->first;
    var_tuple vt{vtemp.vi,reflevel,vtemp.tl}; // clauses to not refer to reflevel but the valid states
                                              // have them so inject them here
    int vi = vt.vi;
    int type = CCTK_GroupTypeFromVarI(vi);
    //std::cout << "-->|Type " << type << " " << CCTK_FullVarName(vi) << "\n";
    if(type != CCTK_GF) {
      //std::cout << "-->|Skip!\n";
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
    if(!on(valid_k[vt],WH_INTERIOR)) // and !silent_psync) 
    {
      // If the read spec is everywhere and we only have
      // interior, that's ok. The system will sync.
      std::ostringstream msg; 
      msg << "Required read for " << vt 
          << " not satisfied. Invalid interior"
          << " at the start of routine " << r;
      dumpValid(msg, vt.vi);
      int level = psync_error ? 0 : 1;
      static bool have_warned = false;
      if(not have_warned) {
        CCTK_WARN(level,msg.str().c_str()); 
        have_warned = true;
      }
    } else if(vt.tl > 0) {
      if(!on(valid_k[vt],i->second)) {
        // If the read spec is everywhere, that's not
        // okay for previous time levels as they won't sync.
        std::ostringstream msg; 
        msg << "Required read for " << vt 
          << " not satisfied. Wanted '" << wstr(i->second)
          << "' found '" << wstr(valid_k[vt])
          << "' at the start of routine " << r;
        dumpValid(msg, vt.vi);
        CCTK_Error(__LINE__,__FILE__,CCTK_THORNSTRING,msg.str().c_str());
      }
    }
    if(i->second == WH_EVERYWHERE) {
      if(on(valid_k[vt],WH_INTERIOR) && valid_k[vt] != WH_EVERYWHERE)
      {
        if(vt.tl != 0) {
          std::ostringstream msg;
          msg << "Attempt to sync previous time level (tl=" << vt.tl << ")"
              << " for " << CCTK_FullVarName(vt.vi);
          CCTK_Error(__LINE__,__FILE__,CCTK_THORNSTRING,msg.str().c_str());
        }

        int g = CCTK_GroupIndexFromVarI(vt.vi);
        pregroups.insert(g);
      } else if(!on(valid_k[vt],WH_INTERIOR))
      {
        std::ostringstream msg; 
        msg << "Cannot sync " << CCTK_FullVarName(vt.vi)
            << " because it is not valid in the interior.";
        dumpValid(msg, vt.vi);
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
 * reverse the order of timelevel information
 */
void flip_rdwr(const cGH *cctkGH, int vi) {
  int const cactus_tl = CCTK_ActiveTimeLevelsVI(cctkGH, vi);
  for(int t = 0; t < cactus_tl-1; t++) {
    var_tuple vt{vi,reflevel,t};
    var_tuple vt_flip{vi,reflevel,cactus_tl-t};
    int tmpdata = valid_k[vt];
    valid_k[vt] = valid_k[vt_flip];
    valid_k[vt_flip] = tmpdata;
  }
}

/**
 * mark a timelvel as invalid
 */
void invalidate_rdwr(const cGH *cctkGH, int vi, int tl) {
  var_tuple vt{vi,reflevel,tl};
  valid_k[vt] = WH_NOWHERE;
}

/**
 * Inverse rotate timelevels:
 * Copy all knowledge about read/write levels up one level,
 * then mark the current level as valid nowhere.
 */
void uncycle_rdwr(const cGH *cctkGH) {
#ifdef PRESYNC_DEBUG
  std::cout << "CYCLE" << std::endl;
#endif
  int num = CCTK_NumVars();
  for(int vi=0;vi<num;vi++) {
    int const cactus_tl = CCTK_ActiveTimeLevelsVI(cctkGH, vi);
    if(cactus_tl > 1) {
      int type = CCTK_GroupTypeFromVarI(vi);
      if(type == CCTK_GF && CCTK_VarTypeSize(CCTK_VarTypeI(vi)) == sizeof(CCTK_REAL)) {
        var_tuple first{vi,reflevel,0};
        int first_valid = valid_k[first]; 
        for(int t = 0; t < cactus_tl-1; t++) {
          var_tuple vold{vi,reflevel,t};
          var_tuple vnew{vi,reflevel,t+1};
          valid_k[vold] = valid_k[vnew];
        }
        var_tuple vt{vi,reflevel,cactus_tl-1};
        valid_k[vt] = first_valid;
      }
    }
  }
#ifdef PRESYNC_DEBUG
  std::cout << "::UNROTATE::" << std::endl;
#endif
}

/**
 * Rotate timelevels:
 * Copy all knowledge about read/write levels down one level,
 * wrapping around at the end.
 */
void cycle_rdwr(const cGH *cctkGH) {
#ifdef PRESYNC_DEBUG
  std::cout << "CYCLE" << std::endl;
#endif
  int num = CCTK_NumVars();
  for(int vi=0;vi<num;vi++) {
    int const cactus_tl = CCTK_ActiveTimeLevelsVI(cctkGH, vi);
    if(cactus_tl > 1) {
      int type = CCTK_GroupTypeFromVarI(vi);
      if(type == CCTK_GF && CCTK_VarTypeSize(CCTK_VarTypeI(vi)) == sizeof(CCTK_REAL)) {
        var_tuple last{vi,reflevel,cactus_tl-1};
        int last_valid = valid_k[last];
        for(int t = cactus_tl - 1; t > 0; t--) {
          var_tuple vold{vi,reflevel,t};
          var_tuple vnew{vi,reflevel,t-1};
          valid_k[vold] = valid_k[vnew];
        }
        var_tuple vt{vi,reflevel,0};
        valid_k[vt] = last_valid;
      }
    }
  }
#ifdef PRESYNC_DEBUG
  std::cout << "::ROTATE::" << std::endl;
#endif
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
  const_cast<cGH*>(cctkGH)->cctk_time = tt->get_time(mglevel, reflevel, timelevel);
  int ierr = SyncProlongateGroups(cctkGH, sync_groups, attribute);
  assert(!ierr);
  const_cast<cGH*>(cctkGH)->cctk_time = old_cctk_time;
  timelevel_offset = old_timelevel_offset;
  timelevel = old_timelevel;
  do_allow_past_timelevels = old_do_allow_past_timelevels;
}

/**
 * Given a variable and a timelevel, set the region
 * of the grid where that variable is valid (i.e. the where_spec).
 */
// TODO: expand to take a reflevel argument?
extern "C" void SetValidRegion(int vi,int tl,int wh) {
  var_tuple vt(vi,reflevel,tl);
  valid_k[vt] = wh;
}

/**
 * Given a variable and a timelevel, return the region
 * of the grid where that variable is valid (i.e. the where_spec).
 */
// TODO: expand to take a reflevel argument?
extern "C" int GetValidRegion(int vi,int tl) {
  var_tuple vt(vi,reflevel,tl);
  return valid_k[vt];
}

/**
 * Apply manual synchronization to variable vi. If this variable
 * is already valid everywhere, this routine does nothing. When
 * the routine finishes, it will be valid everywhere.
 */
extern "C" void ManualSyncGF(CCTK_POINTER_TO_CONST cctkGH_,int tl,int vi) {
  // Do nothing if this is not a GF
  if(CCTK_GroupTypeFromVarI(vi) != CCTK_GF) {
    return;
  }
  const cGH *cctkGH = static_cast<const cGH*>(cctkGH_);
  var_tuple vt{vi,reflevel,tl};
  auto f = valid_k.find(vt);
  if(f == valid_k.end()) {
    dumpValid(std::cerr, vi) << std::endl;
    CCTK_VERROR("Could not find validity information for %s rl=%d tl=%d", CCTK_FullVarName(vi), reflevel, tl);
  }
  assert(f != valid_k.end());
  // Check if anything needs to be done
  if(f->second == WH_EVERYWHERE) {
    return;
  }
  if((f->second & WH_INTERIOR) != WH_INTERIOR) {
    dumpValid(std::cerr, vi) << std::endl;
    CCTK_VERROR("SYNC requires valid data in interior %s rl=%d tl=%d", CCTK_FullVarName(vi), reflevel, tl);
  }
  assert((f->second & WH_INTERIOR) == WH_INTERIOR);

  // Update valid region info
  int gi = CCTK_GroupIndexFromVarI(vi);
  int i0 = CCTK_FirstVarIndexI(gi);
  int iN = i0+CCTK_NumVarsInGroupI(gi);
  for(int vi2=i0;vi2<iN;vi2++) {
    var_tuple vt{vi2,reflevel,tl};
    if(on(valid_k[vt],WH_INTERIOR)) {
      valid_k[vt] = WH_EVERYWHERE;
    }
  }

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
  valid_k[vt] = WH_EVERYWHERE;
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
  const cGH *cctkGH);

typedef CCTK_INT (*sym_iface_boundary_function)(
  CCTK_POINTER_TO_CONST cctkGH);

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

std::map<std::string,Func> boundary_functions;
std::map<std::string,SymFunc> symmetry_functions;
std::map<int,std::vector<Bound>> boundary_conditions;

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

CCTK_INT Carpet_SelectVarForBCI(
    const cGH *cctkGH,
    const CCTK_INT faces,
    const CCTK_INT width,
    CCTK_INT table_handle,
    CCTK_INT var_index,
    const char *bc_name) {
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

extern "C"
CCTK_INT SelectVarForBC(
    const CCTK_POINTER_TO_CONST cctkGH_,
    CCTK_INT faces,
    CCTK_INT width,
    CCTK_INT table_handle,
    const char *var_name,
    const char *bc_name) {
  const cGH *cctkGH = static_cast<const cGH*>(cctkGH_);
  CCTK_INT i = CCTK_VarIndex(var_name);
  return Carpet_SelectVarForBCI(cctkGH,faces,width,table_handle,i,bc_name);
}

extern "C"
CCTK_INT SelectGroupForBC(
    const CCTK_POINTER_TO_CONST cctkGH_,
    CCTK_INT faces,
    CCTK_INT width,
    CCTK_INT table_handle,
    const char *group_name,
    const char *bc_name) {
  const cGH *cctkGH = static_cast<const cGH*>(cctkGH_);
  int group = CCTK_GroupIndex(group_name);
  int vstart = CCTK_FirstVarIndexI(group);
  int vnum   = CCTK_NumVarsInGroupI(group);
  for(int i=vstart;i<vstart+vnum;i++) {
    Carpet_SelectVarForBCI(cctkGH,faces,width,table_handle,i,bc_name);
  }
  return 0;
}

extern "C"
//CCTK_INT Carpet_SelectedGVs(CCTK_POINTER_TO_CONST _GH, CCTK_INT array_size, CCTK_INT *var_indices, CCTK_INT *faces,
//                            CCTK_INT *widths, CCTK_INT *table_handles) {
CCTK_INT Carpet_SelectedGVs() {
  if(bnd_vi == -1) {
#ifdef PRESYNC_DEBUG
    std::cout << "No variable is currently having boundary conditions applied, but a boundary condition is attempting to run." << std::endl;
#endif
    return -1;
  }
  return bnd_vi;
}

/**
 * Apply boundary conditions for a single variable.
 */
extern "C"
void Carpet_ApplyPhysicalBCsForVarI(const cGH *cctkGH, int var_index) {
  DECLARE_CCTK_PARAMETERS;
  if(!use_psync) return;
  auto bc = boundary_conditions;
  if(bc.find(var_index) == bc.end()) {
#ifdef PRESYNC_DEBUG
    std::cout << "ApplyBC: No bc's for " << CCTK_FullVarName(var_index) << std::endl;
#endif
    return;
  }
  bnd_vi = var_index;
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
            ierr = (*fsym.func)(cctkGH);
          }
          var_tuple vt{var_index,reflevel,0};
          valid_k[vt] |= WH_BOUNDARY;
        }
      }
    }
    END_LOCAL_COMPONENT_LOOP;
  }
  END_LOCAL_MAP_LOOP;
  bnd_vi = -1;
}

/**
 * Apply boundary conditions for a group. Called from inside SyncProlongateGroups.
 */
extern "C"
void Carpet_ApplyPhysicalBCsForGroupI(const cGH *cctkGH, int group_index) {
  int vstart = CCTK_FirstVarIndexI(group_index);
  int vnum   = CCTK_NumVarsInGroupI(group_index);
  for(int var_index=vstart;var_index<vstart+vnum;var_index++) {
    Carpet_ApplyPhysicalBCsForVarI(cctkGH,var_index);
  }
}

/**
 * Deprecated: This routine only exists to compare with the old version
 * provided by the Boundary thorn.
 */
extern "C"
void Carpet_ApplyPhysicalBCs(const cGH *cctkGH) {
#ifdef PRESYNC_DEBUG
  std::cout << "Carpet_ApplyPhysicalBCs()" << std::endl;
  for(int n=0;n <= 1;n++) {
    for(auto i=boundary_conditions.begin(); i != boundary_conditions.end(); ++i) {
      std::cout << "n=" << n << " vi=" << i->first << std::endl;
    }
  }
  std::cout << " done" << std::endl;
#endif
}

}

void ShowValid() {
  for(auto i = Carpet::valid_k.begin();i != Carpet::valid_k.end(); ++i) {
    const Carpet::var_tuple& vt = i->first;
    int wh = i->second;
    std::cout << "  valid: " << CCTK_FullVarName(vt.vi) << " tl=" << vt.tl << " wh=" << Carpet::wstr(wh) << std::endl;
  }
}
