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
  
  
  
  void Restrict (const cGH* cgh)
  {
    DECLARE_CCTK_PARAMETERS;
    
    assert (is_level_mode());
    
    if (suppress_restriction) {
      Checkpoint ("Restriction suppressed");
      return;
    }
    
    Checkpoint ("Restrict");
    
    // Restrict
    if (reflevel < reflevels-1) {

      // make the comm_state loop the innermost
      // in order to minimise the number of outstanding communications
      if (minimise_outstanding_communications) {
        for (int group=0; group<CCTK_NumGroups(); ++group) {
          if (CCTK_GroupTypeI(group) == CCTK_GF) {
            if (CCTK_QueryGroupStorageI(cgh, group)) {
              
              const int tl = 0;
              
              for (int m=0; m<(int)arrdata.at(group).size(); ++m) {
                
                // use background time here (which may not be modified
                // by the user)
                const CCTK_REAL time = vtt.at(m)->time (tl, reflevel, mglevel);
                
                const CCTK_REAL time1 = vtt.at(m)->time (0, reflevel, mglevel);
                const CCTK_REAL time2
                  = (cgh->cctk_time - cctk_initial_time) / delta_time;
                assert (fabs(time1 - time2) / (fabs(time1) + fabs(time2) + fabs(cgh->cctk_delta_time)) < 1e-12);
                
                for (int var=0; var<(int)arrdata.at(group).at(m).data.size(); ++var) {
                  for (int c=0; c<vhh.at(m)->components(reflevel); ++c) {
                    for (comm_state state; !state.done(); state.step()) {
                      arrdata.at(group).at(m).data.at(var)->ref_restrict
                        (state, tl, reflevel, c, mglevel, time);
                    }
                  }
                }
              }
              
            } // if group has storage
          } // if grouptype == CCTK_GF
        } // loop over groups
      } else {
        for (comm_state state; !state.done(); state.step()) {
          for (int group=0; group<CCTK_NumGroups(); ++group) {
            if (CCTK_GroupTypeI(group) == CCTK_GF) {
              if (CCTK_QueryGroupStorageI(cgh, group)) {
              
                const int tl = 0;
              
                for (int m=0; m<(int)arrdata.at(group).size(); ++m) {
                
                  // use background time here (which may not be modified
                  // by the user)
                  const CCTK_REAL time = vtt.at(m)->time (tl, reflevel, mglevel);
                
                  const CCTK_REAL time1 = vtt.at(m)->time (0, reflevel, mglevel);
                  const CCTK_REAL time2
                    = (cgh->cctk_time - cctk_initial_time) / delta_time;
                  assert (fabs(time1 - time2) / (fabs(time1) + fabs(time2) + fabs(cgh->cctk_delta_time)) < 1e-12);
                
                  for (int var=0; var<(int)arrdata.at(group).at(m).data.size(); ++var) {
                    for (int c=0; c<vhh.at(m)->components(reflevel); ++c) {
                      arrdata.at(group).at(m).data.at(var)->ref_restrict
                        (state, tl, reflevel, c, mglevel, time);
                    }
                  }
                }
              
              } // if group has storage
            } // if grouptype == CCTK_GF
          } // loop over groups
        } // for state
      } // if minimise_outstanding_communications
    } // if not finest refinement level
    
    
    
    // Sync
    if (reflevel < reflevels-1) {

      // make the comm_state loop the innermost
      // in order to minimise the number of outstanding communications
      if (minimise_outstanding_communications) {
        for (int group=0; group<CCTK_NumGroups(); ++group) {
          if (CCTK_GroupTypeI(group) == CCTK_GF) {
            if (CCTK_QueryGroupStorageI(cgh, group)) {
              
              const int tl = 0;
              
              for (int m=0; m<(int)arrdata.at(group).size(); ++m) {
                for (int var=0; var<(int)arrdata.at(group).at(m).data.size(); ++var) {
                  for (int c=0; c<vhh.at(m)->components(reflevel); ++c) {
                    for (comm_state state; !state.done(); state.step()) {
                      arrdata.at(group).at(m).data.at(var)->sync
                        (state, tl, reflevel, c, mglevel);
                    }
                  }
                }
              }
              
            } // if group has storage
          } // if grouptype == CCTK_GF
        } // loop over groups
      } else {
        for (comm_state state; !state.done(); state.step()) {
          for (int group=0; group<CCTK_NumGroups(); ++group) {
            if (CCTK_GroupTypeI(group) == CCTK_GF) {
              if (CCTK_QueryGroupStorageI(cgh, group)) {
              
                const int tl = 0;
              
                for (int m=0; m<(int)arrdata.at(group).size(); ++m) {
                  for (int var=0; var<(int)arrdata.at(group).at(m).data.size(); ++var) {
                    for (int c=0; c<vhh.at(m)->components(reflevel); ++c) {
                      arrdata.at(group).at(m).data.at(var)->sync
                        (state, tl, reflevel, c, mglevel);
                    }
                  }
                }
              
              } // if group has storage
            } // if grouptype == CCTK_GF
          } // loop over groups
        } // for state
      } // if minimise_outstanding_communications
    } // if not finest refinement level
    
  }
  
} // namespace Carpet

