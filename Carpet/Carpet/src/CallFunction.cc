#include <assert.h>
#include <stdlib.h>

#include "cctk.h"
#include "cctki_GHExtensions.h"

#include "Carpet/CarpetLib/src/gh.hh"

#include "carpet.hh"
  
extern "C" {
  static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/CallFunction.cc,v 1.5 2002/10/24 10:39:37 schnetter Exp $";
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
      // Global operation: call once per refinement level
      
      const int res = CCTK_CallFunction (function, attribute, data);
      assert (res==0);
      
    } else {
      // Local operation: call once per component
      
      BEGIN_LOCAL_COMPONENT_LOOP(cgh) {
	
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
