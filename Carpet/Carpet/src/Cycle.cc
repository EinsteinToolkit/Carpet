#include <assert.h>
#include <stdlib.h>

#include "cctk.h"

#include "Carpet/CarpetLib/src/ggf.hh"
#include "Carpet/CarpetLib/src/gh.hh"

#include "carpet.hh"

static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/Cycle.cc,v 1.1 2001/07/04 12:29:46 schnetter Exp $";



namespace Carpet {
  
  using namespace std;
  
  
  
  void CycleTimeLevels (cGH* cgh)
  {
    Checkpoint ("%*sCycleTimeLevels", 2*reflevel, "");
    
    for (int group=0; group<CCTK_NumGroups(); ++group) {
      if (CCTK_QueryGroupStorageI(cgh, group)) {
	for (int var=0; var<CCTK_NumVarsInGroupI(group); ++var) {
	  const int n = CCTK_FirstVarIndexI(group) + var;
	  switch (CCTK_GroupTypeFromVarI(n)) {
	  case CCTK_SCALAR: {
	    assert (group<(int)scdata.size());
	    assert (var<(int)scdata[group].data.size());
	    const int num_tl = CCTK_NumTimeLevelsFromVarI(n);
	    assert (num_tl>0);
	    void* tmpdata = scdata[group].data[var][reflevel][0];
	    for (int tl=1; tl<num_tl; ++tl) {
	      scdata[group].data[var][reflevel][tl]
		= scdata[group].data[var][reflevel][tl-1];
	    }
	    scdata[group].data[var][reflevel][0] = tmpdata;
	    tmpdata = 0;
	    break;
	  }
	  case CCTK_ARRAY: {
	    assert (group<(int)arrdata.size());
	    assert (var<(int)arrdata[group].data.size());
	    for (int c=0; c<arrdata[group].hh->components(reflevel); ++c) {
	      arrdata[group].data[var]->cycle (reflevel, c, mglevel);
	    }
	    break;
	  }
	  case CCTK_GF: {
	    assert (group<(int)gfdata.size());
	    assert (var<(int)gfdata[group].data.size());
	    for (int c=0; c<hh->components(reflevel); ++c) {
	      gfdata[group].data[var]->cycle (reflevel, c, mglevel);
	    }
	    break;
	  }
	  default:
	    abort();
	  }
	}
      }
    }
  }
  
} // namespace Carpet
