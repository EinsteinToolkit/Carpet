#include <assert.h>
#include <stdlib.h>

#include <algorithm>

#include "cctk.h"
#include "cctki_GHExtensions.h"

#include "gh.hh"

#include "carpet.hh"
  
extern "C" {
  static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/CallFunction.cc,v 1.15 2004/03/11 10:16:50 cott Exp $";
  CCTK_FILEVERSION(Carpet_Carpet_CallFunction_cc);
}



namespace Carpet {
  
  using namespace std;
  
  /// Traverse one function on all components of one refinement level
  /// of one multigrid level.
  int CallFunction (void* function, ///< the function to call
                    cFunctionData* attribute, ///< attributes of the function
                    void* data) ///< ???
  {
//     Checkpoint ("Starting CallFunction...");
    
    cGH* cgh = (cGH*)data;
    
// TODO: disable temporarily
//    if (attribute->meta || is_meta_mode()) {
    if (is_meta_mode()) {
      // Convtest operation
      
      if (do_meta_mode) {
        assert (is_meta_mode());
        Checkpoint ("Meta mode call at %s to %s::%s",
                    attribute->where, attribute->thorn, attribute->routine);
	BEGIN_META_MODE(cgh) {
        const int res = CCTK_CallFunction (function, attribute, data);
        assert (res==0);
	} END_META_MODE;
      }
      
    } else if (attribute->global || is_global_mode()) {
      // Global operation: call once
      
      if (do_global_mode) {
        Checkpoint ("Global mode call at %s to %s::%s",
                    attribute->where, attribute->thorn, attribute->routine);
        BEGIN_GLOBAL_MODE(cgh) {
          const int res = CCTK_CallFunction (function, attribute, data);
          assert (res==0);
        } END_GLOBAL_MODE;
      }
      
    } else if (attribute->level) {
      // Level operation: call once per refinement level
      
      Checkpoint ("Level mode call at %s to %s::%s",
                  attribute->where, attribute->thorn, attribute->routine);
      const int res = CCTK_CallFunction (function, attribute, data);
      assert (res==0);
      
// TODO: disable temporarily
//    } else if (attribute->singlemap) {
    } else if (false) {
      // Single map operation: call once per refinement level and map
      
      BEGIN_MAP_LOOP(cgh, CCTK_GF) {
        
        Checkpoint ("Singlemap mode call at %s to %s::%s",
                    attribute->where, attribute->thorn, attribute->routine);
        const int res = CCTK_CallFunction (function, attribute, data);
        assert (res==0);
        
      }	END_MAP_LOOP;
      
    } else {
      // Local operation: call once per component
      
      BEGIN_MAP_LOOP(cgh, CCTK_GF) {
        BEGIN_LOCAL_COMPONENT_LOOP(cgh, CCTK_GF) {
          
          Checkpoint ("Local mode call at %s to %s::%s",
                      attribute->where, attribute->thorn, attribute->routine);
          const int res = CCTK_CallFunction (function, attribute, data);
          assert (res==0);
          
        } END_LOCAL_COMPONENT_LOOP;
      }	END_MAP_LOOP;
      
    }
    
//     Checkpoint ("done with CallFunction.");
    
    // The return value indicates whether the grid functions have been
    // synchronised.
    // return 0: let the flesh do the synchronisation, if necessary
    return 0;
  }
  
} // namespace Carpet
