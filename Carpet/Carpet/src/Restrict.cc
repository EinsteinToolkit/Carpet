#include <cassert>
#include <cmath>
#include <cstdlib>

#include <cctk.h>
#include <cctk_Parameters.h>

#include <Requirements.hh>

#include <ggf.hh>
#include <gh.hh>

#include <carpet.hh>
#include <Timers.hh>



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
    
#ifdef REQUIREMENTS_HH
    Requirements::Restrict(groups, cctkGH->cctk_iteration, reflevel);
#endif
    
    // Restrict
    {
      static Timer timer ("Restrict");
      timer.start();
      RestrictGroups (cctkGH, groups);
      timer.stop();
    }

    // Synchronise
    {
      static Timer timer ("RestrictSync");
      timer.start();
      SyncGroups (cctkGH, groups);
      timer.stop();
    }
  }
  

  // restrict a set of groups
  static void RestrictGroups (const cGH* cctkGH, const vector<int>& groups) {
    DECLARE_CCTK_PARAMETERS;

    static vector<Timer*> timers;
    if (timers.empty()) {
      timers.push_back(new Timer("comm_state[0].create"));
      for (astate state = static_cast<astate>(0);
           state != state_done;
           state = static_cast<astate>(static_cast<int>(state)+1))
      {
        ostringstream name1;
        name1 << "comm_state[" << timers.size() << "]"
              << "." << tostring(state) << ".user";
        timers.push_back(new Timer(name1.str()));
        ostringstream name2;
        name2 << "comm_state[" << timers.size() << "]"
              << "." << tostring(state) << ".step";
        timers.push_back(new Timer(name2.str()));
      }
    }

    vector<Timer*>::iterator ti = timers.begin();
    (*ti)->start();
    for (comm_state state; not state.done(); state.step()) {
      (*ti)->stop(); ++ti; (*ti)->start();
      for (int group = 0; group < (int)groups.size(); ++group) {
        const int g = groups.AT(group);
        const int active_tl = CCTK_ActiveTimeLevelsGI (cctkGH, g);
        assert (active_tl>=0);
        const int tl = active_tl > 1 ? timelevel : 0;
        for (int m=0; m<(int)arrdata.AT(g).size(); ++m) {
          for (int v = 0; v < (int)arrdata.AT(g).AT(m).data.size(); ++v) {
            ggf *const gv = arrdata.AT(g).AT(m).data.AT(v);
            gv->ref_restrict_all (state, tl, reflevel, mglevel);
          }
        }
      } // loop over groups
      (*ti)->stop(); ++ti; (*ti)->start();
    } // for state
    (*ti)->stop();
    ++ti; assert(ti == timers.end());
  }

} // namespace Carpet

