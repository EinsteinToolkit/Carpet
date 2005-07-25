#include <cassert>
#include <cmath>
#include <cstdlib>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "ggf.hh"
#include "gh.hh"

#include "carpet.hh"



namespace Carpet {
  
  using namespace std;
  
  // restricts a set of groups which all have the same vartype
  static void RestrictGroups (const cGH* cgh, group_set& groups);
  
  
  void Restrict (const cGH* cgh)
  {
    DECLARE_CCTK_PARAMETERS;
    
    assert (is_level_mode());
    
    if (suppress_restriction) {
      Checkpoint ("Restriction suppressed");
      return;
    }
    
    Checkpoint ("Restrict");
    
    // all but the finest level are restricted
    if (reflevel == reflevels-1) {
      return;
    }

    // sort all grid functions into sets of the same vartype
    vector<group_set> groups;

    for (int group = 0; group < CCTK_NumGroups(); ++group) {
      if (CCTK_GroupTypeI(group) == CCTK_GF
          && CCTK_NumVarsInGroupI(group) > 0
          && CCTK_QueryGroupStorageI(cgh, group)) {

        group_set newset;
        const int firstvar = CCTK_FirstVarIndexI (group);
        newset.vartype = CCTK_VarTypeI (firstvar);
        assert (newset.vartype >= 0);
        int c;
        for (c = 0; c < groups.size(); c++) {
          if (newset.vartype == groups[c].vartype) {
            break;
          }
        }
        if (c == groups.size()) {
          groups.push_back (newset);
        }
        groups[c].members.push_back (group);
      }
    }

    // Restrict
    for (int c = 0; c < groups.size(); c++) {
      RestrictGroups (cgh, groups[c]);
    }

    // Synchronise
    for (int c = 0; c < groups.size(); c++) {
      SyncGroups (cgh, groups[c]);
    }
  }
  

  // restricts a set of groups which all have the same vartype
  static void RestrictGroups (const cGH* cgh, group_set& groups) {
    DECLARE_CCTK_PARAMETERS;

    const int tl = 0;

    for (comm_state state(groups.vartype); ! state.done(); state.step()) {
      for (int c = 0; c < groups.members.size(); ++c) {
        const int group = groups.members[c];
        for (int m=0; m<(int)arrdata.at(group).size(); ++m) {

          // use background time here (which may not be modified
          // by the user)
          const CCTK_REAL time = vtt.at(m)->time (tl, reflevel, mglevel);

          const CCTK_REAL time1 = vtt.at(m)->time (0, reflevel, mglevel);
          const CCTK_REAL time2
            = (cgh->cctk_time - cctk_initial_time) / fabs(delta_time);
          assert (fabs(time1 - time2) / (fabs(time1) + fabs(time2) + fabs(cgh->cctk_delta_time)) < 1e-12);

          for (int v = 0; v < arrdata.at(group).at(m).data.size(); ++v) {
            ggf *const gv = arrdata.at(group).at(m).data.at(v);
            for (int c = 0; c < vhh.at(m)->components(reflevel); ++c) {
              gv->ref_restrict (state, tl, reflevel, c, mglevel, time);
            }
          }
        }
      } // loop over groups
    } // for state
  }

} // namespace Carpet

