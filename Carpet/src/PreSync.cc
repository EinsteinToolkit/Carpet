#include <set>
#include <map>
#include <vector>
#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <iostream>
#include <carpet.hh>
#include <math.h>
#include <array>

#define monitor(X) monitor_(X,__FILE__,__LINE__)

namespace Carpet {

// The upper left-hand corner of a 2-D simulation.
// +-----+-----+-----+
// |     |     |     | 
// |  B  |  B  | B+G | 
// |     |     |     | 
// +-----+-----+-----+
// |     |     |     | 
// |  B  |  I  |  G  | 
// |     |     |     | 
// +-----+-----+-----+
// |     |     |     | 
// | B+G |  G  |  G  | 
// |     |     |     | 
// +-----+-----+-----+
//  B  = Boundary 
//  I  = Interior
//  G  = Ghost
// B+G = Boundary and Ghost
//
// Read and Write directives in schedule.ccl apply to:
// Interior, Everywhere, or Invalidate.
//
// Boundary routines are registered with
// to either update B+G cells or not. If it
// updates B+G cells, it runs after synchronization.
// If it doesn't update B+G cells, it runs before.

#define WH_EVERYWHERE          0x7
#define WH_INTERIOR            0x4
#define WH_BOUNDARY            0x2 /* Describes B+W cells */
#define WH_GHOSTS              0x1
#define WH_INVALIDATE_EXTERIOR 0x8

inline bool on(int flags,int flag) {
  return (flags & flag) == flag;
}

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

struct RWClause {
  int where;
  std::string name;
  RWClause() : where(0), name() {}
  ~RWClause() {}
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

void parse_rwclauses(const char *input,std::vector<RWClause>& vec) {
  std::string current;
  RWClause rwc;
  bool parsing_where = false;
  while(true) {
    char c = *input++;
    if(c >= 'A' && c <= 'Z')
      c = c - 'A' + 'a';
    if(iswhite(c)||c==0) {
      if(current != "") {
        rwc.name = current;
        vec.push_back(rwc);
      }
      rwc.name = "";
      rwc.where = 0;
      current = "";
    } else if(c == '(' && !parsing_where) {
      rwc.name = current;
      current = "";
      parsing_where = true;
      assert(rwc.name != "");
    } else if(c == ')' && parsing_where) {
      assert(rwc.name != "");
      if(current == "everywhere") {
        rwc.where = WH_EVERYWHERE;
      } else if(current == "interior-exterior" || current == "in-out") {
        rwc.where = WH_INTERIOR|WH_INVALIDATE_EXTERIOR;
        abort();
      } else if(current == "interior" || current == "in") {
        rwc.where = WH_INTERIOR;
      } else if(current == "boundary") {
        rwc.where = WH_BOUNDARY;
        abort();
      } else if(current == "interior+ghosts" || current == "interiorwithghosts") {
        rwc.where = WH_INTERIOR|WH_GHOSTS;
        abort();
      } else if(current == "boundary+ghosts" || current == "boundarywithghosts") {
        rwc.where = WH_BOUNDARY|WH_GHOSTS;
        abort();
      } else {
        abort();
        // Is it a parameters?
        int pos = current.find("::");
        if(pos >= 0 && rwc.name == "parameter") {
          std::cerr << "Looks like a parameter: (" << current << ")" << std::endl;
          int type;
          std::string thorn = current.substr(0,pos);
          toupper(thorn);
          const void *val = CCTK_ParameterGet(
            current.substr(pos+2).c_str(),
            current.substr(0,pos).c_str(),
            &type);
          if(val != nullptr && type == PARAMETER_STRING) {
            const char* const val_str =
              * ( (const char*const *)val);
            std::cerr << "Parse string: " << val_str << std::endl;
            return;
          }
          std::cerr << "val=" << val << " " << (type == PARAMETER_STRING) << " type=" << type << std::endl;
        }
        std::cerr << "Bad where spec[" << current << "]";
        abort();
      }
      vec.push_back(rwc);
      rwc.name = "";
      rwc.where = 0;
      current = "";
      parsing_where = 0;
    } else {
      current += c;
    }
    if(c == 0)
      break;
  }
}

std::string current_routine;

void compute_clauses(const int num_strings,const char **strings,std::map<int,int>& routine_m) {
  std::vector<RWClause> rwvec;
  for(int i=0;i< num_strings; ++i) {

    // Parse the where clause (if any)
    parse_rwclauses(strings[i],rwvec);
  }
  for(auto i=rwvec.begin();i != rwvec.end();++i) {
    RWClause& rwc = *i;
    int vi = CCTK_VarIndex(rwc.name.c_str());
    int where_val = rwc.where;
    if(where_val == 0) {
      int type = CCTK_GroupTypeFromVarI(vi);
      if(type == CCTK_GF && CCTK_VarTypeSize(CCTK_VarTypeI(vi)) == sizeof(CCTK_REAL)) {
        std::cerr << "Missing Where Spec: Var: " << rwc.name << " Routine: " << current_routine << std::endl;
        assert(0=="missing where spec");
      }
    }
    if(vi >= 0) {
      routine_m[vi] = where_val;
    } else {
      // If looking up a specific variable failed, then
      // we lookup everything on the group.
      int gi = CCTK_GroupIndex(rwc.name.c_str());
      if(gi >= 0) {
        int i0 = CCTK_FirstVarIndexI(gi);
        int iN = i0+CCTK_NumVarsInGroupI(gi);
        for(vi=i0;vi<iN;vi++) {
          routine_m[vi] = where_val;
        }
      } else {
        std::cerr << "error: Could not find (" << rwc.name << ")" << std::endl;
        abort();
      }
    }
  }
}

std::map<std::string,std::map<int,int>> reads,writes;
std::map<int,int> valid_k;
std::map<int,std::string> where_k;

void monitor_(int vi,const char *file,int line) {
  //std::cout << file << ":" << line << ": " << CCTK_FullName(vi) << "=" << wstr(valid_k[vi]) << std::endl;
}

std::map<std::string,std::set<int>> access_allowed;

extern "C" int is_access_allowed(int vi) {
  if(current_routine == "")
    return 1;
  auto& s = access_allowed[current_routine];
  if(s.find(-1) != s.end())
    return 1;
  if(s.find(vi) != s.end())
    return 1;
  return 0;
}

extern "C" void assert_safe() {
  auto& s = access_allowed[current_routine];
  if(s.find(-1) == s.end()) {
    std::cerr << "ASSERT SAFE: (" << current_routine << ")" << std::endl;
    s.insert(-1);
  }
}

void PostCheckValid(cFunctionData *attribute, vector<int> const &sync_groups) {
  std::string r;
  r += attribute->thorn;
  r += "::";
  r += attribute->routine;
  current_routine = "";
  std::map<int,int>& writes_m = writes[r];
  for(auto i = writes_m.begin();i != writes_m.end(); ++i) {
    int vi = i->first;
    int wh = valid_k[vi];
    monitor(vi);
    int wh_before = wh;
    if(on(wh_before,WH_INVALIDATE_EXTERIOR)) {
      wh = WH_INTERIOR;
    } else {
      wh = i->second;
    }
    valid_k[vi] = wh;
    monitor(vi);
    where_k[vi] = r;
  }
  for(auto g = sync_groups.begin();g != sync_groups.end();++g) {
    int gi = *g;
    int i0 = CCTK_FirstVarIndexI(gi);
    int iN = i0+CCTK_NumVarsInGroupI(gi);
    for(int vi=i0;vi<iN;vi++) {
      int& w = writes[r][vi];
      if(w == WH_INTERIOR)
        w = WH_EVERYWHERE;
        valid_k[vi] = w;
        monitor(w);
    }
  }
}

void PreSyncGroups(cFunctionData *attribute,cGH *cctkGH,const std::set<int>& pregroups) {
  std::vector<int> sync_groups;
  for(auto i=pregroups.begin();i != pregroups.end();++i) {
    int gi = *i;
    int i0 = CCTK_FirstVarIndexI(gi);
    int iN = i0+CCTK_NumVarsInGroupI(gi);
    bool push = true;
    for(int vi=i0;vi<iN;vi++) {
      int wh = valid_k[vi];
      monitor(vi);
      if(on(wh,WH_GHOSTS)) {
        continue;
      }
      if(!on(wh,WH_INTERIOR)) {
        std::cerr << "SYNC of variable with invalid interior. Name: " << CCTK_FullName(vi) << " after: " << current_routine << std::endl;
        monitor(vi);
        assert(false);
      }
      if(push) {
        sync_groups.push_back(gi);
        std::cout << "SYNC GROUP: " << CCTK_GroupName(gi) << std::endl;
        push = false;
      }
      valid_k[vi] |= WH_GHOSTS;
      monitor(vi);
    }
  }
  if(sync_groups.size()>0)
    SyncProlongateGroups(cctkGH, sync_groups, attribute);
}

void PreCheckValid(cFunctionData *attribute,cGH *cctkGH,std::set<int>& pregroups) {
  if(cctkGH == 0) return;
  if(attribute == 0) return;

  std::string r;
  r += attribute->thorn;
  r += "::";
  r += attribute->routine;
  current_routine = r;

  if(writes.find(r) == writes.end()) {

    std::map<int,int>& writes_m  = writes[r];
    std::map<int,int>& reads_m  = reads [r];

    compute_clauses(
        attribute->n_WritesClauses,
        attribute->WritesClauses,
        writes_m);
    compute_clauses(
        attribute->n_ReadsClauses,
        attribute->ReadsClauses,
        reads_m);

    for(auto i = writes_m.begin();i != writes_m.end(); ++i) {
      int vi = i->first;
      access_allowed[r].insert(vi);
    }
    for(auto i = reads_m.begin();i != reads_m.end(); ++i) {
      int vi = i->first;
      access_allowed[r].insert(vi);
    }
  }

  std::map<int,int>& reads_m  = reads [r];

  for(auto i = reads_m.begin();i != reads_m.end(); ++i) {
    int vi = i->first;
    monitor(vi);
    if(i->second == WH_EVERYWHERE) {
      if(on(valid_k[vi],WH_INTERIOR) && valid_k[vi] != WH_EVERYWHERE) {
        int g = CCTK_GroupIndexFromVarI(vi);
        pregroups.insert(g);
        //std::cout << "SYNC:: " << r << " " << CCTK_VarName(vi) << " from=" << where_k[vi] << std::endl;
      } else if(!on(valid_k[vi],WH_INTERIOR)) {
        std::cerr << "Cannot sync " << CCTK_FullName(vi) << " because it is not valid in the interior." << std::endl;
        abort();
      }
      monitor(vi);
    }
  }
}

/**
 * TODO: Currently, we only track the validity of the
 * GF on time level zero. We could track based on arbitrary
 * timelevel.
 */
void cycle_rdwr(const cGH *cctkGH) {
  unsigned int num = CCTK_NumVars();
  for(unsigned int vi=0;vi<num;vi++) {
    int const cactus_tl = CCTK_ActiveTimeLevelsVI(cctkGH, vi);
    if(cactus_tl > 1) {
      int type = CCTK_GroupTypeFromVarI(vi);
      if(type == CCTK_GF && CCTK_VarTypeSize(CCTK_VarTypeI(vi)) == sizeof(CCTK_REAL)) {
        valid_k[vi] = 0; // Nowhere valid
        monitor(vi);
      }
    }
  }
}

void Sync1(const cGH *cctkGH,int gi) {
  std::vector<int> sync_groups;
  sync_groups.push_back(gi);
  cFunctionData *attribute = 0;
  SyncProlongateGroups(cctkGH, sync_groups, attribute);
}

extern "C" void ManualSyncGF(const cGH *cctkGH,int vi) {

  auto f = valid_k.find(vi);
  assert(f != valid_k.end());
  monitor(vi);
  // Check if anything needs to be done
  if(f->second == WH_EVERYWHERE) {
    return;
  }
  assert(f->second == WH_INTERIOR);

  // Update valid region info
  int gi = CCTK_GroupIndexFromVarI(vi);
  int i0 = CCTK_FirstVarIndexI(gi);
  int iN = i0+CCTK_NumVarsInGroupI(gi);
  for(int vi2=i0;vi2<iN;vi2++) {
    if(on(valid_k[vi2],WH_INTERIOR)) {
      valid_k[vi2] = WH_EVERYWHERE;
      monitor(vi2);
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
}

typedef void (*boundary_function)(
  const cGH *cctkGH,
  int num_vars,
  int *var_indices,
  int *faces,
  int *widths,
  int *table_handles);

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

std::map<std::string,Func> boundary_functions;
/**
 * The index into the array is the same as the "before"
 * argument defined when registering a BC.
 */
std::array<std::map<int,std::vector<Bound>>,2> boundary_conditions;

extern "C"
void Carpet_RegisterPhysicalBC(
    const cGH *cctkGH,
    boundary_function func,
    const char *bc_name,
    int before) {
  if(before != 0) before = 1;
  Func& f = boundary_functions[bc_name];
  f.func = func;
  f.before = before;
}

extern "C"
void Carpet_SelectVarForBCI(
    const cGH *cctkGH,
    int faces,
    int width,
    int table_handle,
    int var_index,
    const char *bc_name) {
  Func& f = boundary_functions[bc_name];
  assert(var_index != 0);
  std::vector<Bound>& bv = boundary_conditions[f.before][var_index];
  Bound b;
  b.faces = faces;
  b.width = width;
  b.table_handle = table_handle;
  b.bc_name = bc_name;
  bv.push_back(b);
}

extern "C"
void Carpet_SelectGroupForBC(
    const cGH *cctkGH,
    int faces,
    int width,
    int table_handle,
    const char *group_name,
    const char *bc_name) {
  int group = CCTK_GroupIndex(group_name);
  int vstart = CCTK_FirstVarIndexI(group);
  int vnum   = CCTK_NumVarsInGroupI(group);
  for(int i=vstart;i<vstart+vnum;i++) {
    Carpet_SelectVarForBCI(cctkGH,faces,width,table_handle,i,bc_name);
  }
}

extern "C"
void Carpet_ClearBCForVarI(
    const cGH *cctkGH,
    int var_index) {
  assert(var_index != 0);
  boundary_conditions[0][var_index].resize(0);
  boundary_conditions[1][var_index].resize(0);
}

extern "C"
void Carpet_ApplyPhysicalBCsForVarI(const cGH *cctkGH, int var_index,int before) {
  if(before != 0) before = 1;
  auto bc = boundary_conditions[before];
  if(bc.find(var_index) == bc.end()) {
    return;
  } else {
    std::cout << "Carpet_ApplyPhysicalBCsForVarI(" << CCTK_FullName(var_index) << ") -> **FOUND**" << std::endl;
  }
  std::vector<Bound>& bv = bc[var_index];
  BEGIN_LOCAL_MAP_LOOP(cctkGH, CCTK_GF) {
    BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) {
      for(auto j=bv.begin(); j != bv.end(); ++j) {
        Bound& b = *j;
        Func& f = boundary_functions[b.bc_name];

        // Disable faces if they don't apply to this
        // computational region.
        int faces = b.faces;
        for(int i=0;i<6;i++) {
          if(cctkGH->cctk_bbox[i] == 0) {
            faces &= ~(1 << i);
          }
        }

        if(faces != 0)
          (*f.func)(cctkGH,1,&var_index,&faces,&b.width,&b.table_handle);
        valid_k[var_index] |= WH_BOUNDARY;
      }
    }
    END_LOCAL_COMPONENT_LOOP;
  }
  END_LOCAL_MAP_LOOP;
}

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

}
