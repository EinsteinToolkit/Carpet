#include <assert.h>
#include <stdlib.h>

#include "cctk.h"

#include "Carpet/CarpetLib/src/ggf.hh"
#include "Carpet/CarpetLib/src/gh.hh"

#include "carpet.hh"

static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/Cycle.cc,v 1.3 2001/11/02 10:58:58 schnetter Exp $";



namespace Carpet {
  
  using namespace std;
  
  
  
  void CycleTimeLevels (const cGH* cgh)
  {
    Checkpoint ("%*sCycleTimeLevels", 2*reflevel, "");
    
    for (int group=0; group<CCTK_NumGroups(); ++group) {
      if (CCTK_QueryGroupStorageI((cGH*)cgh, group)) {
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
  
} // namespace Carpet
