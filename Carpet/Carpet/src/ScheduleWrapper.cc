#include <cassert>
#include <list>

#include "cctk.h"
#include "cctk_Arguments.h"

#include "carpet.hh"

namespace Carpet {
  
  using namespace std;
  
  
  
  typedef CCTK_INT (* func) (CCTK_POINTER_TO_CONST cctkGH);
  typedef list <func> flist;
  
  static flist func_befores, func_afters;
  
  
  
  extern "C"
  CCTK_INT
  Carpet_RegisterScheduleWrapper (func const func_before,
                                  func const func_after)
  {
    // Add functions
    if (func_before) func_befores.push_back  (func_before);
    if (func_after ) func_afters .push_front (func_after );
    return 0;
  }
  
  extern "C"
  CCTK_INT
  Carpet_UnRegisterScheduleWrapper (func const func_before,
                                    func const func_after)
  {
    // Remove functions
    if (func_before) func_befores.remove (func_before);
    if (func_after ) func_afters .remove (func_after );
    return 0;
  }
  
  
  
  void
  CallBeforeRoutines (cGH const * restrict const cctkGH)
  {
    for (flist::const_iterator
           fli = func_befores.begin(); fli != func_befores.end(); ++ fli)
    {
      (* fli) (cctkGH);
    }
  }
  
  void
  CallAfterRoutines (cGH const * restrict const cctkGH)
  {
    for (flist::const_iterator
           fli = func_afters.begin(); fli != func_afters.end(); ++ fli)
    {
      (* fli) (cctkGH);
    }
  }
  
} // namespace Carpet
