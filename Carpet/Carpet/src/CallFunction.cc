#include <assert.h>
#include <stdlib.h>

#include "cctk.h"
#include "cctki_GHExtensions.h"

#include "Carpet/CarpetLib/src/gh.hh"

#include "carpet.hh"
  
extern "C" {
  static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/CallFunction.cc,v 1.7 2003/05/07 10:03:21 schnetter Exp $";
  CCTK_FILEVERSION(Carpet_Carpet_CallFunction_cc);
}



namespace Carpet {
  
  using namespace std;
  
  int CallFunction (void* function, cFunctionData* attribute, void* data)
  {
    // Traverse one function on all components of one refinement level
    // of one multigrid level
    
    assert (mglevel>=0);
    assert (reflevel>=0);
    
//     Checkpoint ("%*sStarting CallFunction...", 2*reflevel, "");
    
    cGH* cgh = (cGH*)data;
    
    if (attribute->global) {
      // Global operation: call once
      
      if (reflevel==0) {
        Waypoint ("%*sGlobal mode call at %s to %s::%s", 2*reflevel, "",
                  attribute->where, attribute->thorn, attribute->routine);
        const int res = CCTK_CallFunction (function, attribute, data);
        assert (res==0);
      }
      
    } else if (attribute->level) {
      // Level operation: call once per refinement level
      
      Waypoint ("%*sLevel mode call at %s to %s::%s", 2*reflevel, "",
		attribute->where, attribute->thorn, attribute->routine);
      const int res = CCTK_CallFunction (function, attribute, data);
      assert (res==0);
      
    } else {
      // Local operation: call once per component
      
      BEGIN_LOCAL_COMPONENT_LOOP(cgh) {
	
	Waypoint ("%*sLocal mode call on component %d at %s to %s::%s",
		  2*reflevel, "", component,
		  attribute->where, attribute->thorn, attribute->routine);
	const int res = CCTK_CallFunction (function, attribute, data);
	assert (res==0);
	
      }	END_LOCAL_COMPONENT_LOOP(cgh);
      
    }
    
//     Checkpoint ("%*sdone with CallFunction.", 2*reflevel, "");
    
    // The return value indicates whether the grid functions have been
    // synchronised.
    // return 0: let the flesh do the synchronisation, if necessary
    return 0;
  }
  
} // namespace Carpet
