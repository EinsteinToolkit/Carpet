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
              for (int c=0; c<arrdata.at(group).at(m).hh->components(reflevel); ++c) {
                arrdata.at(group).at(m).data.at(var)->cycle (reflevel, c, mglevel);
              }
            }
          }
          break;
          
        case CCTK_SCALAR:
        case CCTK_ARRAY:
          if (do_global_mode) {
	    int const numtimelevels = CCTK_NumTimeLevelsI (group);
	    int const firstvarindex = CCTK_FirstVarIndexI (group);
            for (int var=0; var<CCTK_NumVarsInGroupI(group); ++var) {
              for (int c=0; c<arrdata.at(group).at(0).hh->components(0); ++c) {
                arrdata.at(group).at(0).data.at(var)->cycle (0, c, mglevel);
              }
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
        const int num_vars = CCTK_NumVarsInGroupI(group);
        if (num_vars>0) {
          const int var0 = CCTK_FirstVarIndexI(group);
          assert (var0>=0);
          
          switch (CCTK_GroupTypeI(group)) {
            
          case CCTK_GF:
            for (int m=0; m<(int)arrdata.at(group).size(); ++m) {
              for (int var=0; var<CCTK_NumVarsInGroupI(group); ++var) {
                for (int c=0; c<arrdata.at(group).at(m).hh->components(reflevel); ++c) {
                  arrdata.at(group).at(m).data.at(var)->flip (reflevel, c, mglevel);
                }
              }
            }
            break;
            
          case CCTK_SCALAR:
          case CCTK_ARRAY:
            if (do_global_mode) {
              for (int var=0; var<CCTK_NumVarsInGroupI(group); ++var) {
                for (int c=0; c<arrdata.at(group).at(0).hh->components(0); ++c) {
                  arrdata.at(group).at(0).data.at(var)->flip (0, c, mglevel);
                }
              }
            }
            break;
            
          default:
            assert (0);
          } // switch grouptype
          
        } // if num_vars>0
      } // if storage
    } // for group
  }
  
  
  
  void FillTimeLevels (const cGH* cgh)
  {
    Checkpoint ("CopyTimeLevels");
    assert (is_level_mode());
    
    for (int group=0; group<CCTK_NumGroups(); ++group) {
      if (CCTK_QueryGroupStorageI(cgh, group)) {
        switch (CCTK_GroupTypeI(group)) {
          
        case CCTK_GF:
          assert (reflevel>=0 and reflevel<reflevels);
	  for (int m=0; m<(int)arrdata.at(group).size(); ++m) {
	    for (int var=0; var<CCTK_NumVarsInGroupI(group); ++var) {
              for (int c=0; c<arrdata.at(group).at(m).hh->components(reflevel); ++c) {
                arrdata.at(group).at(m).data.at(var)->fill (reflevel, c, mglevel);
              }
            }
          }
          break;
          
        case CCTK_SCALAR:
        case CCTK_ARRAY:
          if (do_global_mode) {
            for (int var=0; var<CCTK_NumVarsInGroupI(group); ++var) {
              for (int c=0; c<arrdata.at(group).at(0).hh->components(0); ++c) {
                arrdata.at(group).at(0).data.at(var)->fill (0, c, mglevel);
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
  
} // namespace Carpet
