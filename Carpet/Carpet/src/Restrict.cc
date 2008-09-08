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
  static void RestrictGroups (const cGH* cctkGH, const vector<int>& groups);
  
  
  void Restrict (const cGH* cctkGH)
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
    groups.reserve (CCTK_NumGroups());

    for (int group = 0; group < CCTK_NumGroups(); ++group) {
      if (CCTK_GroupTypeI(group) == CCTK_GF
          and CCTK_NumVarsInGroupI(group) > 0
          and CCTK_QueryGroupStorageI(cctkGH, group)) {
        groups.push_back (group);
      }
    }

    // Restrict
    RestrictGroups (cctkGH, groups);

    // Synchronise
    SyncGroups (cctkGH, groups);
  }
  

  // restricts a set of groups which all have the same vartype
  static void RestrictGroups (const cGH* cctkGH, const vector<int>& groups) {
    DECLARE_CCTK_PARAMETERS;

    const int tl = 0;

    for (comm_state state; not state.done(); state.step()) {
      for (int group = 0; group < (int)groups.size(); ++group) {
        const int g = groups[group];
        for (int m=0; m<(int)arrdata.at(g).size(); ++m) {

          // use background time here (which may not be modified
          // by the user)
          const CCTK_REAL time = vtt.at(m)->time (tl, reflevel, mglevel);

          const CCTK_REAL time1 = vtt.at(m)->time (0, reflevel, mglevel);
          const CCTK_REAL time2 =
            (cctkGH->cctk_time - cctk_initial_time) / delta_time;
          const CCTK_REAL time0 =
            abs(time1) + abs(time2) + abs(cctkGH->cctk_delta_time);
          const CCTK_REAL eps = 1.0e-12;
          assert (abs(time1 - time2) <= eps * time0);

          for (int v = 0; v < (int)arrdata.at(g).at(m).data.size(); ++v) {
            ggf *const gv = arrdata.at(g).at(m).data.at(v);
            gv->ref_restrict_all (state, tl, reflevel, mglevel, time);
          }
        }
      } // loop over groups
    } // for state
  }

} // namespace Carpet

