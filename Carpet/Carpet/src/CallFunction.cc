#include <assert.h>
#include <stdlib.h>

#include <algorithm>

#include "cctk.h"
#include "cctki_GHExtensions.h"

#include "gh.hh"

#include "carpet.hh"
  
extern "C" {
  static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/CallFunction.cc,v 1.11 2003/07/22 20:09:05 schnetter Exp $";
  CCTK_FILEVERSION(Carpet_Carpet_CallFunction_cc);
}



namespace Carpet {
  
  using namespace std;
  
  /** Traverse one function on all components of one refinement level
      of one multigrid level.  */
  int CallFunction (void* function, cFunctionData* attribute, void* data)
  {
//     Checkpoint ("%*sStarting CallFunction...", 2*reflevel, "");
    
    cGH* cgh = (cGH*)data;
    
    if (attribute->global || reflevel==-1) {
      // Global operation: call once
      
      if (do_global_mode) {
        assert (component == -1);
        const int saved_mglevel = mglevel;
        if (mglevel!=-1) set_mglevel (cgh, -1);
        const int saved_reflevel = reflevel;
        if (reflevel!=-1) set_reflevel (cgh, -1);
        Waypoint ("Global mode call at %s to %s::%s",
                  attribute->where, attribute->thorn, attribute->routine);
        const int res = CCTK_CallFunction (function, attribute, data);
        assert (res==0);
        if (reflevel!=saved_reflevel) set_reflevel (cgh, saved_reflevel);
        if (mglevel!=saved_mglevel) set_mglevel (cgh, saved_mglevel);
      }
      
    } else if (attribute->level) {
      // Level operation: call once per refinement level
      
      Waypoint ("%*sLevel mode call at %s to %s::%s", 2*reflevel, "",
		attribute->where, attribute->thorn, attribute->routine);
      const int res = CCTK_CallFunction (function, attribute, data);
      assert (res==0);
      
    } else {
      // Local operation: call once per component
      
      BEGIN_LOCAL_COMPONENT_LOOP(cgh, CCTK_GF) {
	
	Waypoint ("%*sLocal mode call on component %d at %s to %s::%s",
		  2*reflevel, "", component,
		  attribute->where, attribute->thorn, attribute->routine);
	const int res = CCTK_CallFunction (function, attribute, data);
	assert (res==0);
	
      }	END_LOCAL_COMPONENT_LOOP;
      
    }
    
//     Checkpoint ("%*sdone with CallFunction.", 2*reflevel, "");
    
    // The return value indicates whether the grid functions have been
    // synchronised.
    // return 0: let the flesh do the synchronisation, if necessary
    return 0;
  }
  
} // namespace Carpet
