#include <cassert>
#include <cstdlib>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "ggf.hh"
#include "gh.hh"

#include "carpet.hh"



namespace Carpet {
  
  using namespace std;
  
  
  static int CheckSyncGroupConsistency ( const cGH* cgh,
                                          const char *groupname );
  static void ProlongateGroupBoundaries ( const cGH* cgh,
                                          CCTK_REAL initial_time, int group );
  
  int SyncGroup (const cGH* cgh, const char* groupname)
  {
    DECLARE_CCTK_PARAMETERS;
    int retval = 0;
    const int group = CCTK_GroupIndex(groupname);
    assert (group>=0 and group<CCTK_NumGroups());
    assert (group<(int)arrdata.size());
    const int grouptype = CCTK_GroupTypeI(group);
    assert(grouptype == CCTK_GF ||
           grouptype == CCTK_SCALAR || grouptype == CCTK_ARRAY);
    assert (CCTK_NumVarsInGroupI(group) != 0);

    Checkpoint ("SyncGroup \"%s\" time=%g", groupname, (double)cgh->cctk_time);
    
    retval = CheckSyncGroupConsistency (cgh, groupname);
    
    if (retval == 0) {
    
      // Prolongate the boundaries
      if (do_prolongate && grouptype == CCTK_GF) {
        assert (reflevel>=0 and reflevel<reflevels);
        if (reflevel > 0) {
          ProlongateGroupBoundaries ( cgh, cctk_initial_time, group );
        }
      }
    
      // Sync
      SyncGVGroup ( cgh, group );
    }
    return retval;
  }
  
  void ProlongateGroupBoundaries ( const cGH* cgh, CCTK_REAL initial_time,
                                  int group )
  {
    DECLARE_CCTK_PARAMETERS;

    // use the current time here (which may be modified by the user)
    const CCTK_REAL time = (cgh->cctk_time - initial_time) / delta_time;
    const int tl = 0;
    
    // make the comm_state loop the innermost
    // in order to minimise the number of outstanding communications
    if (minimise_outstanding_communications) {
      for (int m=0; m<(int)arrdata.at(group).size(); ++m) {
        for (int var=0; var<CCTK_NumVarsInGroupI(group); ++var) {
          for (int c=0; c<vhh.at(m)->components(reflevel); ++c) {
            for (comm_state state; !state.done(); state.step()) {
              arrdata.at(group).at(m).data.at(var)->ref_bnd_prolongate
                (state, tl, reflevel, c, mglevel, time);
            }
          }
        }
      }
    } else {
      for (comm_state state; !state.done(); state.step()) {
        for (int m=0; m<(int)arrdata.at(group).size(); ++m) {
          for (int var=0; var<CCTK_NumVarsInGroupI(group); ++var) {
            for (int c=0; c<vhh.at(m)->components(reflevel); ++c) {
              arrdata.at(group).at(m).data.at(var)->ref_bnd_prolongate
                (state, tl, reflevel, c, mglevel, time);
            }
          }
        }
      }
    }
  }

  void SyncGVGroup ( const cGH* cgh, int group )
  {
    DECLARE_CCTK_PARAMETERS;
    const int tl = 0;

    // make the comm_state loop the innermost
    // in order to minimise the number of outstanding communications
    if (minimise_outstanding_communications) {
      for (int m=0; m<(int)arrdata.at(group).size(); ++m) {
                        CCTK_GroupName (group));
        for (int var=0; var<(int)arrdata.at(group).at(m).data.size(); ++var) {
          for (int c=0; c<vhh.at(m)->components(reflevel); ++c) {
            for (comm_state state; ! state.done(); state.step()) {
              arrdata.at(group).at(m).data.at(var)->sync
                (state, tl, reflevel, c, mglevel);
            }
          }
        }
                        CCTK_GroupName (group));
      }
    } else {
      for (comm_state state; ! state.done(); state.step()) {
        for (int m=0; m<(int)arrdata.at(group).size(); ++m) {
          for (int var=0; var<(int)arrdata.at(group).at(m).data.size(); ++var) {
            for (int c=0; c<vhh.at(m)->components(reflevel); ++c) {
              arrdata.at(group).at(m).data.at(var)->sync
                (state, tl, reflevel, c, mglevel);
            }
          }
        }
      }
    }
  }

  int CheckSyncGroupConsistency ( const cGH* cgh,const char *groupname )
  {
    int retval = 0;
    const int group = CCTK_GroupIndex(groupname);
    const int grouptype = CCTK_GroupTypeI(group);
    
    if (grouptype == CCTK_GF) {
      if (reflevel == -1) {
        CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Cannot synchronise in global mode "
                    "(Tried to synchronise group \"%s\")",
                    groupname);
      }
      if (map != -1 and component == -1) {
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
        if (maps == 1 and vhh.at(map)->local_components(reflevel) == 1) {
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
      retval = -1;
    }
    return retval;
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
