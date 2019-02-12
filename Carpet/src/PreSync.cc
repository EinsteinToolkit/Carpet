#include <set>
#include <cstring>
#include <map>
#include <vector>
#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <iostream>
#include <carpet.hh>
#include <math.h>
#include <array>
#include "variables.hh"
#include "modes.hh"
#include "PreSync.h"

extern "C" void ShowValid();

namespace Carpet {

bool msg1 = true, msg2 = true, msg3 = true;

struct var_tuple {
  int vi; // var index
  int rl; // refinement level
  int tl; // time level;
  var_tuple() : vi(-1), rl(-1), tl(-1) {}
  var_tuple(int vi_,int rl_,int tl_) : vi(vi_), rl(rl_), tl(tl_) {}
};

std::ostream& operator<<(std::ostream& o,const var_tuple& vt) {
  o << CCTK_FullVarName(vt.vi);
  for(int i=0;i<vt.tl;i++)
    o << "_p";
  return o;
}

bool operator<(const var_tuple& v1,const var_tuple& v2) {
  if(v1.vi < v2.vi) return true;
  if(v1.vi > v2.vi) return false;
  if(v1.rl < v2.rl) return true;
  if(v1.rl > v2.rl) return false;
  if(v1.tl < v2.tl) return true;
  if(v1.tl > v2.tl) return false;
  return false;
}

inline bool on(int flags,int flag) {
  return (flags & flag) == flag;
}

inline void cctk_assert_(int line,const char *file,const char *thorn,const char *str) {
  std::ostringstream msg; 
  msg << "Assertion Failed: " << str;
  CCTK_Error(line,file,thorn,msg.str().c_str());
}

#define CCTK_ASSERT(X) if(!(X)) cctk_assert_(__LINE__,__FILE__,CCTK_THORNSTRING,#X);

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

/*
 * A read/write clause.
 */
struct RWClause {
  int where; 
  int tl; // time-level
  std::string name;
  RWClause() : where(0), tl(0), name() {}
  ~RWClause() {}
  void parse_tl() {
    int s = name.size();
    while(s > 2 && name[s-1]=='p' && name[s-2] == '_') {
      tl++;
      s -= 2;
      name.resize(s);
    }
  }
};

inline void tolower(std::string& s) {
  for(auto si=s.begin();si != s.end();++si) {
    if(*si >= 'A' && *si <= 'Z') {
      *si = *si + 'a' - 'A';
    }
  }
}

inline void toupper(std::string& s) {
  for(auto si=s.begin();si != s.end();++si) {
    if(*si >= 'a' && *si <= 'z') {
      *si = *si + 'A' - 'a';
    }
  }
}

inline bool iswhite(char c) {
  return c == ' ' || c == '\n' || c == '\r' || c == '\t' || c == '\b' || c == '\f';
}

inline std::ostream& operator<<(std::ostream& o,const RWClause& rwc) {
  return o << "RWC(" << rwc.name << ",tl=" << rwc.tl << ",wh=" << wstr(rwc.where) << ")";
}

/**
 * Parse a single input clause, e.g. READS: foo::bar(Everywhere)
 * and store the result in a vector of RWClause objects. Note
 * that this method supports "All" as a shorthand for "Everywhere"
 * and "In" as a shorthand for "Interior."
 */
void parse_rwclauses(const char *input,std::vector<RWClause>& vec,int default_where) {
  std::string current;
  RWClause rwc;
  bool parsing_where = false;
  const char *input_sav = input;
  bool parsed = false;
  while(true) {
    char c = *input++;
    if(c >= 'A' && c <= 'Z')
      c = c - 'A' + 'a';
    if(c == 0) {
      rwc.name = current;
      if(rwc.name != "") {
        rwc.where = default_where;
        vec.push_back(rwc);
        parsed = true;
      }
      break;
    } else if(iswhite(c)) {
      CCTK_ASSERT(false);
    } else if(c == '(' && !parsing_where) {
      rwc.name = current;
      rwc.parse_tl();
      current = "";
      parsing_where = true;
      CCTK_ASSERT(rwc.name != "");
    } else if(c == ')' && parsing_where) {
      CCTK_ASSERT(rwc.name != "");
      if(current == "everywhere" || current == "all") {
        rwc.where = WH_EVERYWHERE;
      } else if(current == "interior" || current == "in") {
        rwc.where = WH_INTERIOR;
      } else if(current == "boundary") {
        rwc.where = WH_BOUNDARY;
      } else {
        abort();
      }
      vec.push_back(rwc);
      parsed = true;
      break;
      parsing_where = 0;
    } else {
      current += c;
    }
  }
  if(!parsed) std::cout << "PARSE FAIL: (" << input_sav << ") (" << input << ") " << rwc<< std::endl;
  CCTK_ASSERT(parsed);
}

std::string current_routine;

/**
 * This method is convert an array of textual descriptions of read/write clauses to
 * a map that takes a var tuple (var index, time level) and produces a where spec.
 */
void compute_clauses(const int num_strings,const char **strings,std::map<var_tuple,int>& routine_m,int default_where) {
  std::vector<RWClause> rwvec;
  for(int i=0;i< num_strings; ++i) {

    // Parse the where clause (if any)
    parse_rwclauses(strings[i],rwvec,default_where);
  }
  for(auto i=rwvec.begin();i != rwvec.end();++i) {
    RWClause& rwc = *i;
    int vi = CCTK_VarIndex(rwc.name.c_str());
    int where_val = rwc.where;
    if(where_val == 0) {
      int type = CCTK_GroupTypeFromVarI(vi);
      if(type == CCTK_GF && CCTK_VarTypeSize(CCTK_VarTypeI(vi)) == sizeof(CCTK_REAL)) {
        std::ostringstream msg;
        msg << "Missing Where Spec: Var: " << rwc.name << " Routine: " << current_routine;
        CCTK_Error(__LINE__,__FILE__,CCTK_THORNSTRING,msg.str().c_str());
      }
    }

    if(vi >= 0) {
      var_tuple vt{vi, -1, rwc.tl};
      routine_m[vt] = where_val;
    } else {
      // If looking up a specific variable failed, then
      // we lookup everything on the group.
      int gi = CCTK_GroupIndex(rwc.name.c_str());
      if(gi >= 0) {
        if(!CCTK_NumVarsInGroupI(gi)) {
          // Some groups have array length set at runtime, and
          // they can be set to 0 (inactive), in which case
          // the rd/wr should be ignored.
          continue;
        }
        int i0 = CCTK_FirstVarIndexI(gi);
        int iN = i0+CCTK_NumVarsInGroupI(gi);
        for(vi=i0;vi<iN;vi++) {
          var_tuple vt{vi,-1,0};
          routine_m[vt] = where_val;
        }
        assert(i0 < iN);
      } else {
        std::string vname = rwc.name.c_str();
        vname += "[0]";
        vi = CCTK_VarIndex(vname.c_str());
        gi = CCTK_GroupIndexFromVarI(vi);
        if(gi >= 0) {
          int i0 = CCTK_FirstVarIndexI(gi);
          int iN = i0+CCTK_NumVarsInGroupI(gi);
          for(vi=i0;vi<iN;vi++) {
            var_tuple vt{vi,-1,0};
            routine_m[vt] = where_val;
          }
          assert(i0 < iN);
        }
      }
    }
  }
}

std::map<std::string,std::map<var_tuple,int>> reads,writes;
std::map<var_tuple,int> valid_k;
std::map<var_tuple,int> old_valid_k;

extern "C" void diagnosticPreValid();
void diagnosticPreValid() {
  old_valid_k = valid_k;
}

extern "C" void diagnosticChanged();
void diagnosticChanged() {
  for(auto entry = valid_k.begin(); valid_k.end() != entry; ++entry) {
    auto f = old_valid_k.find(entry->first);
    if(f == old_valid_k.end()) {
      std::cout << " NEW: " << entry->first << " = " << wstr(entry->second) << std::endl;
    } else if(entry->second != f->second) {
      std::cout << " CHANGED: " << entry->first << " = " << wstr(f->second) << " -> " << wstr(entry->second) << std::endl;
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

void PostCheckValid(cFunctionData *attribute, vector<int> const &sync_groups) {
  std::string r;
  r += attribute->thorn;
  r += "::";
  r += attribute->routine;
  current_routine = "";
  std::map<var_tuple,int>& writes_m = writes[r];
  for(auto i = writes_m.begin();i != writes_m.end(); ++i) {
    const var_tuple& vi = i->first;
    valid_k[vi] |= i->second;
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
      }
    }
  }
  const char *vname = "TMUNUBASE::eTtt";
  static int vi_ = CCTK_VarIndex(vname);
  static var_tuple vt_{vi_,-1,0};
  std::cout << "   " << vname << " := " << wstr(valid_k[vt_]) << std::endl;
}

/**
 * Do the actual presync of the groups.
 **/
void PreSyncGroups(cFunctionData *attribute,cGH *cctkGH,const std::set<int>& pregroups) {
  DECLARE_CCTK_PARAMETERS;
  std::vector<int> sync_groups;
  bool use_psync = (CCTK_ParameterValInt("use_psync","Cactus") != 0);
  bool psync_error = (CCTK_ParameterValInt("psync_error","Cactus") != 0);

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
          if(msg1) {
            msg1 = false;
            CCTK_WARN(level,msg.str().c_str());
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
    for (size_t sgi=0;sgi<sync_groups.size();sgi++) {
      int i0 = CCTK_FirstVarIndexI(sync_groups[sgi]);
      int iN = i0+CCTK_NumVarsInGroupI(sync_groups[sgi]);
      for (int vi=i0;vi<iN;vi++) {
        std::ostringstream msg;
        msg << "  Presync: syncing for ";
        msg << attribute->thorn << "::" << attribute->routine
          << " in/at " << attribute->where
          << " variable " << CCTK_FullVarName(vi);
        CCTK_INFO(msg.str().c_str());
      }
    }
    if(use_psync) {
      SyncProlongateGroups(cctkGH, sync_groups, attribute);
      for (size_t sgi=0;sgi<sync_groups.size();sgi++) {
        int i0 = CCTK_FirstVarIndexI(sync_groups[sgi]);
        int iN = i0+CCTK_NumVarsInGroupI(sync_groups[sgi]);
        for (int vi=i0;vi<iN;vi++) {
          if(valid_k[var_tuple(vi,-1,0)] != WH_EVERYWHERE) {
            std::ostringstream msg;
            msg << "Not Valid Everywhere:" << wstr(valid_k[var_tuple(vi,-1,0)]) << " " << CCTK_FullVarName(vi);
            msg << " Routine: " << attribute->thorn << "::" << attribute->routine;
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
  bool psync_error = (CCTK_ParameterValInt("psync_error","Cactus") != 0);
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
      std::cout << "tmp acccess for " << CCTK_FullVarName(var_index) << std::endl;
      return true;
    }
    if(hasAccess(tmp_write,vi)) {
      std::cout << "tmp acccess for " << CCTK_FullVarName(var_index) << std::endl;
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
  bool use_psync = (CCTK_ParameterValInt("use_psync","Cactus") != 0);
  bool psync_error = (CCTK_ParameterValInt("psync_error","Cactus") != 0);
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

    compute_clauses(
        attribute->n_WritesClauses,
        attribute->WritesClauses,
        writes_m,WH_INTERIOR);
    compute_clauses(
        attribute->n_ReadsClauses,
        attribute->ReadsClauses,
        reads_m,WH_EVERYWHERE);
  }

  std::map<var_tuple,int>& reads_m  = reads [r];

  for(auto i = reads_m.begin();i != reads_m.end(); ++i) {
    var_tuple vt = i->first;
    vt.rl = reflevel; // clauses to not refer to reflevel but the valid states
                      // have them so inject them here
    int vi = vt.vi;
    int type = CCTK_GroupTypeFromVarI(vi);
    //std::cout << "-->|Type " << type << " " << CCTK_FullVarName(vi) << "\n";
    if(type != CCTK_GF) {
      //std::cout << "-->|Skip!\n";
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
      int level = psync_error ? 0 : 1;
      if(msg2) {
        CCTK_WARN(level,msg.str().c_str()); 
        msg2 = false;
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
      } else if(!on(valid_k[vt],WH_INTERIOR)) // and !silent_psync) 
      {
        std::ostringstream msg; 
        msg << "Cannot sync " << CCTK_FullVarName(vt.vi)
            << " because it is not valid in the interior.";
        int level = psync_error ? 0 : 1;
        if(msg3) {
          CCTK_WARN(level,msg.str().c_str()); 
          msg3 = false;
        }
      }
    }
  }
}

/**
 * Rotate timelevels:
 * Copy all knowledge about read/write levels down one level,
 * then mark the current level as valid nowhere.
 */
void cycle_rdwr(const cGH *cctkGH) {
  int num = CCTK_NumVars();
  for(int vi=0;vi<num;vi++) {
    int const cactus_tl = CCTK_ActiveTimeLevelsVI(cctkGH, vi);
    if(cactus_tl > 1) {
      int type = CCTK_GroupTypeFromVarI(vi);
      if(type == CCTK_GF && CCTK_VarTypeSize(CCTK_VarTypeI(vi)) == sizeof(CCTK_REAL)) {
        for(int t = cactus_tl - 1; t > 0; t--) {
          var_tuple vold{vi,reflevel,t};
          var_tuple vnew{vi,reflevel,t-1};
          valid_k[vold] = valid_k[vnew];
        }
        var_tuple vt{vi,reflevel,0};
        valid_k[vt] = WH_NOWHERE;
      }
    }
  }
  std::cout << "::ROTATE::" << std::endl;
}

/**
 * Called by ManualSyncGF.
 */
void Sync1(const cGH *cctkGH,int gi) {
  std::vector<int> sync_groups;
  sync_groups.push_back(gi);
  cFunctionData *attribute = 0;
  int ierr = SyncProlongateGroups(cctkGH, sync_groups, attribute);
  assert(!ierr);
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
extern "C" void ManualSyncGF(CCTK_POINTER_TO_CONST cctkGH_,int vi) {
  const cGH *cctkGH = static_cast<const cGH*>(cctkGH_);
  var_tuple vt{vi,reflevel,0};
  auto f = valid_k.find(vt);
  std::cout << "ManualSyncGF(" << CCTK_FullVarName(vi) << ")" << std::endl;
  CCTK_ASSERT(f != valid_k.end());
  // Check if anything needs to be done
  if(f->second == WH_EVERYWHERE) {
    return;
  }
  CCTK_ASSERT(f->second == WH_INTERIOR);

  // Update valid region info
  int gi = CCTK_GroupIndexFromVarI(vi);
  int i0 = CCTK_FirstVarIndexI(gi);
  int iN = i0+CCTK_NumVarsInGroupI(gi);
  for(int vi2=i0;vi2<iN;vi2++) {
    var_tuple vt{vi2,reflevel,0};
    if(on(valid_k[vt],WH_INTERIOR)) {
      valid_k[vt] = WH_EVERYWHERE;
    }
  }

  // Take action by mode
  if(is_level_mode()) {
    Sync1(cctkGH,gi);
  } else if(is_global_mode()) {
    BEGIN_REFLEVEL_LOOP(cctkGH) {
      Sync1(cctkGH,gi);
    }
    END_REFLEVEL_LOOP;
  } else if(is_meta_mode()) {
    BEGIN_MGLEVEL_LOOP(cctkGH) {
      BEGIN_REFLEVEL_LOOP(cctkGH) {
        Sync1(cctkGH,gi);
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

struct Bound {
  std::string bc_name;
  int faces;
  int width;
  int table_handle;
};

struct Func {
  boundary_function func;
  int before;
};

struct SymFunc {
  boundary_function func;
  int handle;
  int faces;
  int width;
};

std::map<std::string,Func> boundary_functions;
std::map<std::string,SymFunc> symmetry_functions;
/**
 * The index into the array is the same as the "before"
 * argument defined when registering a BC.
 */
std::array<std::map<int,std::vector<Bound>>,2> boundary_conditions;

extern "C"
CCTK_INT Carpet_RegisterPhysicalBC(
    CCTK_POINTER_TO_CONST /*cctkGH_*/,
    CCTK_POINTER_TO_CONST func_,
    CCTK_INT before,
    CCTK_STRING bc_name) {
  boundary_function func = reinterpret_cast<boundary_function>(func_);
  if(before != 0) before = 1;
  if(NULL==func) {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,  
               "Physical Boundary condition '%s' points to NULL.", bc_name);
  }
  Func& f = boundary_functions[bc_name];
  f.func = func;
  f.before = before;
  return 0;
}

extern "C"
CCTK_INT RegisterSymmetryBC(
    const CCTK_POINTER_TO_CONST /*cctkGH_*/,
    boundary_function func,
    const char *bc_name,
    int handle,
    int faces,
    int width) {
  if(NULL==func) {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,  
               "Symmetry Boundary condition '%s' points to NULL.", bc_name);
  }
  SymFunc& f = symmetry_functions[bc_name];
  f.func = func;
  f.handle = handle;
  f.faces = faces;
  f.width = width;
  return 0;
}

CCTK_INT Carpet_SelectVarForBCI(
    const cGH *cctkGH,
    int faces,
    int width,
    int table_handle,
    int var_index,
    const char *bc_name) {
  if(!boundary_functions.count(bc_name)) {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,  
               "Requested BC '%s' not found.", bc_name);
  }
  Func& f = boundary_functions.at(bc_name);
  CCTK_ASSERT(var_index != 0);
  std::vector<Bound>& bv = boundary_conditions[f.before][var_index];
  Bound b;
  b.faces = faces;
  b.width = width;
  b.table_handle = table_handle;
  b.bc_name = bc_name;
  bv.push_back(b);
  return 0;
}

extern "C"
CCTK_INT SelectVarForBC(
    const CCTK_POINTER_TO_CONST cctkGH_,
    int faces,
    int width,
    int table_handle,
    const char *var_name,
    const char *bc_name) {
  const cGH *cctkGH = static_cast<const cGH*>(cctkGH_);
  int i = CCTK_VarIndex(var_name);
  return Carpet_SelectVarForBCI(cctkGH,faces,width,table_handle,i,bc_name);
}

extern "C"
CCTK_INT SelectGroupForBC(
    const CCTK_POINTER_TO_CONST cctkGH_,
    int faces,
    int width,
    int table_handle,
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

/**
 * Apply boundary conditions for a single variable.
 */
extern "C"
void Carpet_ApplyPhysicalBCsForVarI(const cGH *cctkGH, int var_index,int before) {
  if(before != 0) before = 1;
  auto bc = boundary_conditions[before];
  if(bc.find(var_index) == bc.end()) {
    return;
  }
  std::vector<Bound>& bv = bc[var_index];
  BEGIN_LOCAL_MAP_LOOP(cctkGH, CCTK_GF) {
    BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) {
      for(auto j=bv.begin(); j != bv.end(); ++j) {
        Bound& b = *j;
        Func& f = boundary_functions.at(b.bc_name);

        // Disable faces if they don't apply to this
        // computational region.
        int faces = b.faces;
        #if 0
        for(int i=0;i<6;i++) {
          if(cctkGH->cctk_bbox[i] == 0) {
            faces &= ~(1 << i);
          }
        }
        #endif
        if(faces != 0) {
          int ierr = (*f.func)(cctkGH,1,&var_index,&faces,&b.width,&b.table_handle);
          assert(!ierr);
          for (auto iter = symmetry_functions.begin(); iter != symmetry_functions.end(); iter++) {
            std::string name = iter->first;
            SymFunc& fsym = symmetry_functions.at(name);
//            std::cout << name << " BC applied to " << CCTK_FullVarName(var_index) << std::endl;
            ierr = (*fsym.func)(cctkGH,1,&var_index,&fsym.faces,&fsym.width,&fsym.handle);
          }
          var_tuple vt{var_index,reflevel,0};
          valid_k[vt] |= WH_BOUNDARY;
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
extern "C"
void Carpet_ApplyPhysicalBCsForGroupI(const cGH *cctkGH, int group_index,int before) {
  if(before != 0) before = 1;
  int vstart = CCTK_FirstVarIndexI(group_index);
  int vnum   = CCTK_NumVarsInGroupI(group_index);
  for(int var_index=vstart;var_index<vstart+vnum;var_index++) {
    Carpet_ApplyPhysicalBCsForVarI(cctkGH,var_index,before);
  }
}

/**
 * Deprecated: This routine only exists to compare with the old version
 * provided by the Boundary thorn.
 */
extern "C"
void Carpet_ApplyPhysicalBCs(const cGH *cctkGH) {
  std::cout << "Carpet_ApplyPhysicalBCs()" << std::endl;
  for(int n=0;n <= 1;n++) {
    for(auto i=boundary_conditions[n].begin(); i != boundary_conditions[n].end(); ++i) {
      std::cout << "n=" << n << " vi=" << i->first << std::endl;
    }
  }
  std::cout << " before:" << std::endl;
  for(auto i=boundary_conditions[1].begin(); i != boundary_conditions[1].end(); ++i) {
    Carpet_ApplyPhysicalBCsForVarI(cctkGH,i->first,1); // before
  }
  std::cout << " after:" << std::endl;
  for(auto i=boundary_conditions[0].begin(); i != boundary_conditions[0].end(); ++i) {
    Carpet_ApplyPhysicalBCsForVarI(cctkGH,i->first,0); // after
  }
  std::cout << " done" << std::endl;
}

extern "C"
CCTK_INT Bdry2_Boundary_RegisterSymmetryBC(
    CCTK_POINTER_TO_CONST /*cctkGH_*/,
    iface_boundary_function func_,
    CCTK_INT handle,
    CCTK_INT faces,
    CCTK_INT width,
    CCTK_STRING bc_name) {
  boundary_function func = reinterpret_cast<boundary_function>(func_);
  if(NULL==func) {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,  
               "Symmetry Boundary condition '%s' points to NULL.", bc_name);
  }
  SymFunc& f = symmetry_functions[bc_name];
  f.func = func;
  f.handle = handle;
  f.faces = faces;
//  &f.width = width;
//  std::cout << "Register Width of " << f.width[0] << " for " << bc_name << std::endl;
  return 0;
}

}

void ShowValid() {
  for(auto i = Carpet::valid_k.begin();i != Carpet::valid_k.end(); ++i) {
    const Carpet::var_tuple& vt = i->first;
    int wh = i->second;
    std::cout << "  valid: " << CCTK_FullVarName(vt.vi) << " tl=" << vt.tl << " wh=" << Carpet::wstr(wh) << std::endl;
  }
}
