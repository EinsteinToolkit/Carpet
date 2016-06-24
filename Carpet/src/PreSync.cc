#include <set>
#include <map>
#include <vector>
#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <iostream>
#include <carpet.hh>
#include <math.h>

namespace Carpet {

#define WH_EVERYWHERE          0x7
#define WH_INTERIOR            0x4
#define WH_BOUNDARY            0x2
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
      } else if(current == "interior" || current == "in") {
        rwc.where = WH_INTERIOR;
      } else if(current == "boundary") {
        rwc.where = WH_BOUNDARY;
      } else if(current == "boundary+ghosts" || current == "boundarywithghosts") {
        rwc.where = WH_BOUNDARY|WH_GHOSTS;
      } else {
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
    std::cerr << "RWC: " << rwc.name << " " << wstr(rwc.where) << std::endl;
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

std::map<std::string,std::set<int>> access_allowed;

extern "C" int is_access_allowed(int vi) {
  if(current_routine == "")
    return 1;
  auto& s = access_allowed[current_routine];
  if(s.find(-1) != s.end())
    return 1;
  if(s.find(vi) != s.end())
    return 1;
  //std::cerr << "BLOCKED: " << CCTK_FullName(vi) << "/" << CCTK_GroupNameFromVarI(vi) << " in " << current_routine << std::endl;
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
    int wh_before = wh;
    if(on(wh_before,WH_INVALIDATE_EXTERIOR)) {
      wh = WH_INTERIOR;
    } else {
      wh |= i->second;
    }
    #if 1
    //if(wh_before != 0 && wh != wh_before) 
    if(i->second == WH_INTERIOR && !on(wh_before,WH_INTERIOR) && (on(wh_before,WH_BOUNDARY)||on(wh_before,WH_GHOSTS)))
    {
      std::cerr << "UPGRADE: " << CCTK_FullName(vi) << " from " << wstr(wh_before) << " to " << wstr(wh) << " using " << wstr(i->second) << " in " << r << std::endl;
    }
    #endif
    valid_k[vi] = wh;
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
    }
  }
}

void PreSyncGroups(cFunctionData *attribute,cGH *cctkGH,const std::set<int>& pregroups) {
  std::vector<int> sync_groups;
  for(auto i=pregroups.begin();i != pregroups.end();++i) {
    int gi = *i;
    sync_groups.push_back(gi);
    int i0 = CCTK_FirstVarIndexI(gi);
    int iN = i0+CCTK_NumVarsInGroupI(gi);
    for(int vi=i0;vi<iN;vi++) {
      if(on(valid_k[vi],WH_INTERIOR)) {
        std::cerr << "SYNC of variable with invalid interior. Name: " << CCTK_FullName(vi) << " after: " << current_routine << std::endl;
        assert(false);
      }
      valid_k[vi] = WH_EVERYWHERE;
    }
  }
  if(not sync_groups.empty())
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
    if(i->second == WH_EVERYWHERE) {
      if(valid_k[vi] == WH_INTERIOR) {
        int g = CCTK_GroupIndexFromVarI(vi);
        pregroups.insert(g);
        std::cerr << "SYNC:: " << r << " " << CCTK_VarName(vi) << " from=" << where_k[vi] << std::endl;
      }
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
      }
    }
  }
}

}
