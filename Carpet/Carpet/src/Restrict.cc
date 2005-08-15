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
  
  // restricts a set of groups
  static void RestrictGroups (const cGH* cgh, const vector<int>& groups);
  
  
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

    // remove all groups with are non-GFs, empty, or have no storage assigned
    vector<int> groups;

    for (int group = 0; group < CCTK_NumGroups(); ++group) {
      if (CCTK_GroupTypeI(group) == CCTK_GF
          && CCTK_NumVarsInGroupI(group) > 0
          && CCTK_QueryGroupStorageI(cgh, group)) {
        groups.push_back (group);
      }
    }

    // Restrict
    RestrictGroups (cgh, groups);

    // Synchronise
    SyncGroups (cgh, groups);
  }
  

  // restricts a set of groups which all have the same vartype
  static void RestrictGroups (const cGH* cgh, const vector<int>& groups) {
    DECLARE_CCTK_PARAMETERS;

    const int tl = 0;

    for (comm_state state; ! state.done(); state.step()) {
      for (int c = 0; c < groups.size(); ++c) {
        const int group = groups[c];
        for (int m=0; m<(int)arrdata.at(group).size(); ++m) {

          // use background time here (which may not be modified
          // by the user)
          const CCTK_REAL time = vtt.at(m)->time (tl, reflevel, mglevel);

          const CCTK_REAL time1 = vtt.at(m)->time (0, reflevel, mglevel);
          const CCTK_REAL time2
            = (cgh->cctk_time - cctk_initial_time) / delta_time;
          assert (abs(time1 - time2) / (abs(time1) + abs(time2) + abs(cgh->cctk_delta_time)) < 1e-12);

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

