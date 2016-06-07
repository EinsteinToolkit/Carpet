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

  #define WH_EVERYWHERE 0x11
  #define WH_INTERIOR   0x01
  #define WH_EXTERIOR   0x10

inline void tolower(std::string& s) {
  for(auto si=s.begin();si != s.end();++si) {
    if(*si >= 'A' && *si <= 'Z') {
      *si = *si + 'a' - 'A';
    }
  }
}

void compute_clauses(const int num_strings,const char **strings,std::map<int,int>& routine_m) {
  for(int i=0;i< num_strings; ++i) {

    // Parse the where clause (if any)
    const char *clause = strings[i];
    const char *openp  = strchr(clause,'(');
    const char *closep = strchr(clause,')');
    const char *end_impl = strchr(clause,':');
    if(end_impl == 0) {
      std::cerr << "bad str=" << strings[i] << std::endl;
    }
    assert(end_impl != 0);
    std::string where, str;
    int where_val=0;
    if(openp != 0) {
      str.assign(clause,openp-clause);
      where.assign(openp+1,closep-openp-1);
      tolower(where);
      std::string imp(clause,end_impl-clause);
      tolower(imp);

      if(where == "everywhere") {
        where_val = WH_EVERYWHERE;
      } else if(where == "interior") {
        where_val = WH_INTERIOR;
      } else if(where == "boundary") {
        where_val = WH_EXTERIOR;
      } else {
        std::cout << "error in where clause for " << str << "=" << where << std::endl;
        assert(false);
      }
    } else {
      str = clause;
      where_val = WH_EVERYWHERE;
    }

    int vi = CCTK_VarIndex(str.c_str());
    if(vi >= 0) {
      routine_m[vi] = where_val;
    } else {
      // If looking up a specific variable failed, then
      // we lookup everything on the group.
      int gi = CCTK_GroupIndex(str.c_str());
      if(gi >= 0) {
        int i0 = CCTK_FirstVarIndexI(gi);
        int iN = i0+CCTK_NumVarsInGroupI(gi);
        for(vi=i0;vi<iN;vi++) {
          routine_m[vi] = where_val;
        }
      } else {
        std::cerr << "error: Could not find (" << str << ")" << std::endl;
      }
    }
  }
}

std::map<std::string,std::map<int,int>> reads,writes;
std::map<std::string,std::map<int,int>> valid;
std::map<int,std::string> where_k;

void PostSyncGroups(cFunctionData *attribute, vector<int> const &sync_groups) {
  std::string r;
  r += attribute->thorn;
  r += "::";
  r += attribute->routine;
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

void PreSyncGroups(cFunctionData *attribute,cGH *cctkGH,std::set<int>& pregroups) {
  if(cctkGH == 0) return;
  if(attribute == 0) return;

  if(pregroups.size() == 0) {
    std::string r;
    r += attribute->thorn;
    r += "::";
    r += attribute->routine;

    std::map<int,int>& writes_m = writes[r];
    std::map<int,int>& reads_m  = reads [r];

    if(writes_m.size()==0 && reads_m.size()==0) {
      compute_clauses(
          attribute->n_WritesClauses,
          attribute->WritesClauses,
          writes_m);
      compute_clauses(
          attribute->n_ReadsClauses,
          attribute->ReadsClauses,
          reads_m);
    }

    std::map<int,int>& valid_k = valid[r];
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
    for(auto i = writes_m.begin();i != writes_m.end(); ++i) {
      int vi = i->first;
      int wh = i->second;
      valid_k[vi] = wh;
      where_k[vi] = r;
    }
  }
  if(pregroups.size()>0) {
    std::vector<int> sync_groups;
    for(auto i=pregroups.begin();i != pregroups.end();++i) {
      sync_groups.push_back(*i);
    }
    SyncProlongateGroups(cctkGH, sync_groups, attribute);
  }
}

}
