#include <assert.h>
#include <stdlib.h>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "ggf.hh"
#include "gh.hh"

#include "carpet.hh"

extern "C" {
  static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/Comm.cc,v 1.25 2004/02/03 14:34:34 schnetter Exp $";
  CCTK_FILEVERSION(Carpet_Carpet_Comm_cc);
}



namespace Carpet {
  
  using namespace std;
  
  
  
  int SyncGroup (const cGH* cgh, const char* groupname)
  {
    DECLARE_CCTK_PARAMETERS;
    
    const int group = CCTK_GroupIndex(groupname);
    assert (group>=0 && group<CCTK_NumGroups());
    assert (group<(int)arrdata.size());
    
    Checkpoint ("SyncGroup %s time=%g", groupname, (double)cgh->cctk_time);
    
    const int grouptype = CCTK_GroupTypeI(group);
    
    if (grouptype == CCTK_GF) {
      if (reflevel == -1) {
        CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Cannot synchronise in global mode "
                    "(Tried to synchronise group \"%s\")",
                    groupname);
      }
      if (map != -1) {
        if (maps == 1) {
          CCTK_VWarn (2, __LINE__, __FILE__, CCTK_THORNSTRING,
                      "Synchronising group \"%s\" in singlemap mode",
                      groupname);
        } else {
          CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                      "Cannot synchronise in singlemap mode "
                      "(Tried to synchronise group \"%s\")",
                      groupname);
        }
      }
      if (component != -1) {
        if (maps == 1 && vhh.at(map)->local_components(reflevel) == 1) {
          CCTK_VWarn (2, __LINE__, __FILE__, CCTK_THORNSTRING,
                      "Synchronising group \"%s\" in local mode",
                      groupname);
        } else {
          CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                      "Cannot synchronise in local mode "
                      "(Tried to synchronise group \"%s\")",
                      groupname);
        }
      }
    }
    
    if (! CCTK_QueryGroupStorageI(cgh, group)) {
      CCTK_VWarn (2, __LINE__, __FILE__, CCTK_THORNSTRING,
		  "Cannot synchronise group \"%s\" because it has no storage",
		  groupname);
      return -1;
    }
    
    if (CCTK_NumVarsInGroupI(group) == 0) return 0;
    
    const int n0 = CCTK_FirstVarIndexI(group);
    assert (n0>=0);
    const int num_tl = CCTK_NumTimeLevelsFromVarI(n0);
    assert (num_tl>0);
    const int tl = 0;
    
    // Prolongate the boundaries
    if (do_prolongate) {
      switch (grouptype) {
        
      case CCTK_GF:
        assert (reflevel>=0 && reflevel<reflevels);
        if (reflevel > 0) {
          
          // use the current time here (which may be modified by the
          // user)
          const CCTK_REAL time = ((cgh->cctk_time - cctk_initial_time)
                                  / (delta_time * mglevelfact));
          
          for (comm_state<dim> state; !state.done(); state.step()) {
            for (int m=0; m<maps; ++m) {
              for (int var=0; var<CCTK_NumVarsInGroupI(group); ++var) {
                for (int c=0; c<vhh.at(m)->components(reflevel); ++c) {
                  arrdata.at(group).at(m).data.at(var)->ref_bnd_prolongate
                    (state, tl, reflevel, c, mglevel, time);
                }
              }
            }
          } // for state
        } // if reflevel>0
        break;
        
      case CCTK_SCALAR:
      case CCTK_ARRAY:
        // do nothing
        break;
        
      default:
        assert (0);
      } // switch grouptype
    } // if do_prolongate
    
    // Sync
    for (comm_state<dim> state; !state.done(); state.step()) {
      switch (CCTK_GroupTypeI(group)) {
        
      case CCTK_GF:
        for (int m=0; m<maps; ++m) {
          for (int var=0; var<CCTK_NumVarsInGroupI(group); ++var) {
            for (int c=0; c<vhh.at(m)->components(reflevel); ++c) {
              arrdata.at(group).at(m).data.at(var)->sync
                (state, tl, reflevel, c, mglevel);
            }
          }
        }
        break;
        
      case CCTK_SCALAR:
      case CCTK_ARRAY:
        for (int var=0; var<(int)arrdata.at(group).at(0).data.size(); ++var) {
          arrdata.at(group).at(0).data.at(var)->sync (state, 0, 0, 0, 0);
        }
        break;
        
      default:
        assert (0);
      } // switch grouptype
    } // for state
    
    return 0;
  }
  
  
  
  int EnableGroupComm (const cGH* cgh, const char* groupname)
  {
    // Communication is always enabled
    return 0;
  }
  
  int DisableGroupComm (const cGH* cgh, const char* groupname)
  {
    // Communication is always enabled
    return -1;
  }
  
} // namespace Carpet
