#include <assert.h>
#include <stdlib.h>

#include "cctk.h"

#include "Carpet/CarpetLib/src/ggf.hh"
#include "Carpet/CarpetLib/src/gh.hh"

#include "carpet.hh"

static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/Restrict.cc,v 1.5 2001/12/09 16:41:53 schnetter Exp $";



namespace Carpet {
  
  using namespace std;
  
  
  
  void Restrict (const cGH* cgh)
  {
    assert (component == -1);
    
    Checkpoint ("%*sRestrict", 2*reflevel, "");
    
    if (reflevel == hh->reflevels()-1) return;
    
    for (int group=0; group<CCTK_NumGroups(); ++group) {
      
      // Restrict only groups with storage
      if (CCTK_QueryGroupStorageI(cgh, group)) {
	
	const int tl = 0;
	
	for (int var=0; var<(int)arrdata[group].data.size(); ++var) {
	  for (int c=0; c<hh->components(reflevel); ++c) {
	    arrdata[group].data[var]->ref_restrict (tl, reflevel, c, mglevel);
	  }
	  for (int c=0; c<arrdata[group].hh->components(reflevel); ++c) {
	    arrdata[group].data[var]->sync (tl, reflevel, c, mglevel);
	  }
	  
	  for (int c=0; c<arrdata[group].hh->components(reflevel+1); ++c) {
	    arrdata[group].data[var]->ref_bnd_prolongate
	      (tl, reflevel+1, c, mglevel);
	  }
	  // TODO: is this necessary?
	  for (int c=0; c<arrdata[group].hh->components(reflevel+1); ++c) {
	    arrdata[group].data[var]->sync (tl, reflevel+1, c, mglevel);
	  }
	}
	
      }	// if group has storage
      
    } // loop over groups
  }
  
} // namespace Carpet
