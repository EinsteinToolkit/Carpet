#include <assert.h>
#include <stdlib.h>

#include "cctk.h"

#include "Carpet/CarpetLib/src/ggf.hh"
#include "Carpet/CarpetLib/src/gh.hh"

#include "carpet.hh"

extern "C" {
  static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/Cycle.cc,v 1.12 2002/11/16 19:10:50 schnetter Exp $";
  CCTK_FILEVERSION(Carpet_Carpet_Cycle_cc);
}



namespace Carpet {
  
  using namespace std;
  
  
  
  void CycleTimeLevels (const cGH* cgh)
  {
    Checkpoint ("%*sCycleTimeLevels", 2*reflevel, "");
    
    for (int group=0; group<CCTK_NumGroups(); ++group) {
      if (reflevel<arrdata[group].hh->reflevels()
	  && CCTK_QueryGroupStorageI(cgh, group)) {
	for (int var=0; var<CCTK_NumVarsInGroupI(group); ++var) {
	  
	  assert (group<(int)arrdata.size());
	  assert (var<(int)arrdata[group].data.size());
	  for (int c=0; c<arrdata[group].hh->components(reflevel); ++c) {
	    arrdata[group].data[var]->cycle (reflevel, c, mglevel);
	  }
	  
	}
      }
    }
  }
  
  
  
  void FlipTimeLevels (const cGH* cgh)
  {
    Checkpoint ("FlipTimeLevels");
    
    for (int group=0; group<CCTK_NumGroups(); ++group) {
      if (CCTK_QueryGroupStorageI(cgh, group)) {
	
	const int var0 = CCTK_FirstVarIndexI(group);
	const int num_tl = CCTK_NumTimeLevelsFromVarI(var0);
	switch (num_tl) {
	case 1:
	  // Do nothing
	  break;
	case 3:
	  // Flip
	  for (int var=0; var<CCTK_NumVarsInGroupI(group); ++var) {
	    
	    assert (group<(int)arrdata.size());
	    assert (var<(int)arrdata[group].data.size());
	    for (int rl=0; rl<hh->reflevels(); ++rl) {
	      for (int c=0; c<arrdata[group].hh->components(rl); ++c) {
		arrdata[group].data[var]->flip (rl, c, mglevel);
	      }
	    }
	    
	  }
	  break;
	default:
	  // Error
	  assert (0);
	} // switch
	
      }
    }
  }
  
} // namespace Carpet
