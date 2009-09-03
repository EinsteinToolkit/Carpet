#include <cassert>
#include <cstdlib>

#include <cctk.h>

#include <ggf.hh>
#include <gh.hh>

#include <carpet.hh>



namespace Carpet {
  
  using namespace std;
  
  
  
  void CycleTimeLevels (const cGH* cgh)
  {
    Checkpoint ("CycleTimeLevels");
    assert (is_level_mode());
    
    for (int group=0; group<CCTK_NumGroups(); ++group) {
      if (CCTK_QueryGroupStorageI(cgh, group)) {
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
		  cgh->data[varindex][tl]
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
  
  
  
  void FlipTimeLevels (const cGH* cgh)
  {
    Checkpoint ("FlipTimeLevels");
    assert (is_level_mode());
    
    for (int group=0; group<CCTK_NumGroups(); ++group) {
      if (CCTK_QueryGroupStorageI(cgh, group)) {
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
		  cgh->data[varindex][tl]
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
  
  
  
  void FillTimeLevels (const cGH* cgh)
  {
    Checkpoint ("FillTimeLevels");
    assert (is_level_mode());
    
    for (int group=0; group<CCTK_NumGroups(); ++group) {
      if (CCTK_QueryGroupStorageI(cgh, group)) {
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
