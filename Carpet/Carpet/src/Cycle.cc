#include <cassert>
#include <cstdlib>

#include "cctk.h"

#include "ggf.hh"
#include "gh.hh"

#include "carpet.hh"



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
	  for (int m=0; m<(int)arrdata.at(group).size(); ++m) {
	    for (int var=0; var<CCTK_NumVarsInGroupI(group); ++var) {
              arrdata.at(group).at(m).data.at(var)->
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
              arrdata.at(group).at(0).data.at(var)->cycle_all (0, mglevel);
	      {
                int const varindex = firstvarindex + var;
		int const c = CCTK_MyProc(cgh);
		for (int tl=0; tl<numtimelevels; ++tl) {
		  cgh->data[varindex][tl]
		    = (tl < groupdata.at(group).info.activetimelevels
		       ? ((*arrdata.at(group).at(0).data.at(var))
			  (tl, 0, c, 0)->storage())
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
          for (int m=0; m<(int)arrdata.at(group).size(); ++m) {
            for (int var=0; var<CCTK_NumVarsInGroupI(group); ++var) {
              arrdata.at(group).at(m).data.at(var)->
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
              arrdata.at(group).at(0).data.at(var)->flip_all (0, mglevel);
	      {
                int const varindex = firstvarindex + var;
		int const c = CCTK_MyProc(cgh);
		for (int tl=0; tl<numtimelevels; ++tl) {
		  cgh->data[varindex][tl]
		    = (tl < groupdata.at(group).info.activetimelevels
		       ? ((*arrdata.at(group).at(0).data.at(var))
			  (tl, 0, c, 0)->storage())
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
	  for (int m=0; m<(int)arrdata.at(group).size(); ++m) {
	    for (int var=0; var<CCTK_NumVarsInGroupI(group); ++var) {
              arrdata.at(group).at(m).data.at(var)->
                fill_all (reflevel, mglevel);
            }
          }
          break;
          
        case CCTK_SCALAR:
        case CCTK_ARRAY:
          if (do_global_mode) {
            for (int var=0; var<CCTK_NumVarsInGroupI(group); ++var) {
              arrdata.at(group).at(0).data.at(var)->fill_all (0, mglevel);
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
