#include <assert.h>
#include <stdlib.h>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "ggf.hh"
#include "gh.hh"

#include "carpet.hh"

extern "C" {
  static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/Comm.cc,v 1.23 2003/11/05 16:18:37 schnetter Exp $";
  CCTK_FILEVERSION(Carpet_Carpet_Comm_cc);
}



namespace Carpet {
  
  using namespace std;
  
  
  
  int SyncGroup (const cGH* cgh, const char* groupname)
  {
    DECLARE_CCTK_PARAMETERS;
    
    const int group = CCTK_GroupIndex(groupname);
    assert (group>=0 && group<CCTK_NumGroups());
    
    Checkpoint ("%*sSyncGroup %s time=%g", 2*reflevel, "",
                groupname, cgh->cctk_time);
    
    const int grouptype = CCTK_GroupTypeI(group);
    
    if (grouptype == CCTK_GF) {
      if (reflevel == -1) {
        CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Cannot synchronise grid functions in global mode "
                    "(Tried to synchronise group \"%s\")",
                    groupname);
      }
      if (hh->local_components(reflevel) != 1 && component != -1) {
        CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Cannot synchronise grid functions in local mode "
                    "(Tried to synchronise group \"%s\")",
                    groupname);
      }
      if (hh->local_components(reflevel) != 1) assert (component == -1);
    }
    
    if (! CCTK_QueryGroupStorageI(cgh, group)) {
      CCTK_VWarn (2, __LINE__, __FILE__, CCTK_THORNSTRING,
		  "Cannot synchronise group \"%s\" because it has no storage",
		  groupname);
      return -1;
    }
    
    if (CCTK_NumVarsInGroupI(group)==0) return 0;
    
    const int n0 = CCTK_FirstVarIndexI(group);
    assert (n0>=0);
    const int num_tl = CCTK_NumTimeLevelsFromVarI(n0);
    assert (num_tl>0);
    const int tl = 0;
    
    // Prolongate the boundaries
    if (reflevel>0) {
      if (grouptype == CCTK_GF) {
        if (do_prolongate) {
          if (arrdata[group].do_transfer) {
            for (comm_state<dim> state; !state.done(); state.step()) {
              assert (group<(int)arrdata.size());
              for (int var=0; var<(int)arrdata[group].data.size(); ++var) {
                // use the current time here (which may be modified by
                // the user)
                const CCTK_REAL time = (cgh->cctk_time - cctk_initial_time) / delta_time;
#if 0
                const CCTK_REAL time1 = tt->time (tl, reflevel, mglevel);
                const CCTK_REAL time2 = (cgh->cctk_time - cctk_initial_time) / delta_time;
                assert (fabs((time1 - time2) / (fabs(time1) + fabs(time2) + fabs(cgh->cctk_delta_time))) < 1e-12);
#endif
                for (int c=0; c<arrdata[group].hh->components(reflevel); ++c) {
                  arrdata[group].data[var]->ref_bnd_prolongate
                    (state, tl, reflevel, c, mglevel, time);
                }
              } // for var
            } // for state
          } else {
            Checkpoint ("%*s(no prolongating for group %s)",
                        2*reflevel, "", groupname);
          } // if ! do_transfer
        } // if do_prolongate
      } // if grouptype == CCTK_GF
    } // if reflevel>0
    
    // Sync
    for (comm_state<dim> state; !state.done(); state.step()) {
      assert (group<(int)arrdata.size());
      for (int var=0; var<(int)arrdata[group].data.size(); ++var) {
        if (grouptype == CCTK_GF) {
          for (int c=0; c<arrdata[group].hh->components(reflevel); ++c) {
            arrdata[group].data[var]->sync (state, tl, reflevel, c, mglevel);
          }
        } else {                // grouptype != CCTK_GF
          arrdata[group].data[var]->sync (state, 0, 0, 0, 0);
        } // grouptype != CCTK_GF
      } // for var
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
