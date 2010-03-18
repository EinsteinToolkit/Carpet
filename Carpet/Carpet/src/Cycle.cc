#include <cassert>
#include <cstdlib>

#include <cctk.h>
#include <cctk_Parameters.h>

#include <ggf.hh>
#include <gh.hh>

#include <carpet.hh>



namespace Carpet {
  
  using namespace std;
  
  
  
  void CycleTimeLevels (cGH* const cctkGH)
  {
    DECLARE_CCTK_PARAMETERS;
    
    Checkpoint ("CycleTimeLevels");
    assert (is_level_mode());
    
    assert (timelevel == 0);
    tt->advance_time (mglevel, reflevel);
    cctkGH->cctk_time = tt->get_time (mglevel, reflevel, timelevel);
    
    for (int group=0; group<CCTK_NumGroups(); ++group) {
      if (CCTK_QueryGroupStorageI(cctkGH, group)) {
        
        int const activetimelevels = CCTK_ActiveTimeLevelsGI (cctkGH, group);
        if (activetimelevels > 1) {
          if (activetimelevels < prolongation_order_time+1) {
            char * const groupname = CCTK_GroupName (group);
            CCTK_VWarn (CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                        "Group \"%s\" has %d only active time levels.  Groups with more than one active time level need at least %d active time levels for prolongation_order_time=%d",
                        groupname,
                        activetimelevels, int(prolongation_order_time+1),
                        int(prolongation_order_time));
            free (groupname);
          }
        }
        
        switch (CCTK_GroupTypeI(group)) {
          
        case CCTK_GF:
          assert (reflevel>=0 and reflevel<reflevels);
	  for (int m=0; m<(int)arrdata.AT(group).size(); ++m) {
	    for (int var=0; var<CCTK_NumVarsInGroupI(group); ++var) {
              arrdata.AT(group).AT(m).data.AT(var)->
                cycle_all (reflevel, mglevel);
            }
          }
          break;
          
        case CCTK_SCALAR:
        case CCTK_ARRAY:
          if (do_global_mode) {
	    int const numtimelevels = CCTK_NumTimeLevelsI (group);
	    int const firstvarindex = CCTK_FirstVarIndexI (group);
            for (int var=0; var<CCTK_NumVarsInGroupI(group); ++var) {
              arrdata.AT(group).AT(0).data.AT(var)->cycle_all (0, mglevel);
	      {
                int const varindex = firstvarindex + var;
		for (int tl=0; tl<numtimelevels; ++tl) {
		  cctkGH->data[varindex][tl]
		    = (tl < groupdata.AT(group).info.activetimelevels
		       ? ((*arrdata.AT(group).AT(0).data.AT(var))
			  (tl, 0, 0, 0)->storage())
		       : NULL);
		}
	      }
            }
          }
          break;
          
        default:
          assert (0);
        } // switch grouptype
      } // if storage
    } // for group
  }
  
  
  
  void UncycleTimeLevels (cGH* const cctkGH)
  {
    DECLARE_CCTK_PARAMETERS;
    
    Checkpoint ("UncycleTimeLevels");
    assert (is_level_mode());
    
    assert (timelevel == 0);
    tt->retreat_time (mglevel, reflevel);
    cctkGH->cctk_time = tt->get_time (mglevel, reflevel, timelevel);
    
    for (int group=0; group<CCTK_NumGroups(); ++group) {
      if (CCTK_QueryGroupStorageI(cctkGH, group)) {
        
        int const activetimelevels = CCTK_ActiveTimeLevelsGI (cctkGH, group);
        if (activetimelevels > 1) {
          if (activetimelevels < prolongation_order_time+1) {
            char * const groupname = CCTK_GroupName (group);
            CCTK_VWarn (CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                        "Group \"%s\" has %d only active time levels.  Groups with more than one active time level need at least %d active time levels for prolongation_order_time=%d",
                        groupname,
                        activetimelevels, int(prolongation_order_time+1),
                        int(prolongation_order_time));
            free (groupname);
          }
        }
        
        switch (CCTK_GroupTypeI(group)) {
          
        case CCTK_GF:
          assert (reflevel>=0 and reflevel<reflevels);
	  for (int m=0; m<(int)arrdata.AT(group).size(); ++m) {
	    for (int var=0; var<CCTK_NumVarsInGroupI(group); ++var) {
              arrdata.AT(group).AT(m).data.AT(var)->
                uncycle_all (reflevel, mglevel);
            }
          }
          break;
          
        case CCTK_SCALAR:
        case CCTK_ARRAY:
          if (do_global_mode) {
	    int const numtimelevels = CCTK_NumTimeLevelsI (group);
	    int const firstvarindex = CCTK_FirstVarIndexI (group);
            for (int var=0; var<CCTK_NumVarsInGroupI(group); ++var) {
              arrdata.AT(group).AT(0).data.AT(var)->uncycle_all (0, mglevel);
	      {
                int const varindex = firstvarindex + var;
		for (int tl=0; tl<numtimelevels; ++tl) {
		  cctkGH->data[varindex][tl]
		    = (tl < groupdata.AT(group).info.activetimelevels
		       ? ((*arrdata.AT(group).AT(0).data.AT(var))
			  (tl, 0, 0, 0)->storage())
		       : NULL);
		}
	      }
            }
          }
          break;
          
        default:
          assert (0);
        } // switch grouptype
      } // if storage
    } // for group
  }
  
  
  
  void FlipTimeLevels (cGH* const cctkGH)
  {
    DECLARE_CCTK_PARAMETERS;
    
    Checkpoint ("FlipTimeLevels");
    assert (is_level_mode());
    
    assert (timelevel == 0);
    tt->flip_timelevels (mglevel, reflevel);
    cctkGH->cctk_time = tt->get_time (mglevel, reflevel, timelevel);
    cctkGH->cctk_delta_time *= -1;
    
    for (int group=0; group<CCTK_NumGroups(); ++group) {
      if (CCTK_QueryGroupStorageI(cctkGH, group)) {
        
        int const activetimelevels = CCTK_ActiveTimeLevelsGI (cctkGH, group);
        if (activetimelevels > 1) {
          if (activetimelevels < prolongation_order_time+1) {
            char * const groupname = CCTK_GroupName (group);
            CCTK_VWarn (CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                        "Group \"%s\" has %d only active time levels.  Groups with more than one active time level need at least %d active time levels for prolongation_order_time=%d",
                        groupname,
                        activetimelevels, int(prolongation_order_time+1),
                        int(prolongation_order_time));
            free (groupname);
          }
        }
        
        switch (CCTK_GroupTypeI(group)) {
          
        case CCTK_GF:
          assert (reflevel>=0 and reflevel<reflevels);
          for (int m=0; m<(int)arrdata.AT(group).size(); ++m) {
            for (int var=0; var<CCTK_NumVarsInGroupI(group); ++var) {
              arrdata.AT(group).AT(m).data.AT(var)->
                flip_all (reflevel, mglevel);
            }
          }
          break;
          
        case CCTK_SCALAR:
        case CCTK_ARRAY:
          if (do_global_mode) {
	    int const numtimelevels = CCTK_NumTimeLevelsI (group);
	    int const firstvarindex = CCTK_FirstVarIndexI (group);
            for (int var=0; var<CCTK_NumVarsInGroupI(group); ++var) {
              arrdata.AT(group).AT(0).data.AT(var)->flip_all (0, mglevel);
	      {
                int const varindex = firstvarindex + var;
		for (int tl=0; tl<numtimelevels; ++tl) {
		  cctkGH->data[varindex][tl]
		    = (tl < groupdata.AT(group).info.activetimelevels
		       ? ((*arrdata.AT(group).AT(0).data.AT(var))
			  (tl, 0, 0, 0)->storage())
		       : NULL);
		}
	      }
            }
          }
          break;
          
        default:
          assert (0);
        } // switch grouptype
      } // if storage
    } // for group
  }
  
  
  
  void FillTimeLevels (const cGH* const cctkGH)
  {
    Checkpoint ("FillTimeLevels");
    assert (is_level_mode());
    
    for (int group=0; group<CCTK_NumGroups(); ++group) {
      if (CCTK_QueryGroupStorageI(cctkGH, group)) {
        switch (CCTK_GroupTypeI(group)) {
          
        case CCTK_GF:
          assert (reflevel>=0 and reflevel<reflevels);
	  for (int m=0; m<(int)arrdata.AT(group).size(); ++m) {
	    for (int var=0; var<CCTK_NumVarsInGroupI(group); ++var) {
              arrdata.AT(group).AT(m).data.AT(var)->
                fill_all (reflevel, mglevel);
            }
          }
          break;
          
        case CCTK_SCALAR:
        case CCTK_ARRAY:
          if (do_global_mode) {
            for (int var=0; var<CCTK_NumVarsInGroupI(group); ++var) {
              arrdata.AT(group).AT(0).data.AT(var)->fill_all (0, mglevel);
            }
          }
          break;
          
        default:
          assert (0);
        } // switch grouptype
      } // if storage
    } // for group
  }
  
} // namespace Carpet
