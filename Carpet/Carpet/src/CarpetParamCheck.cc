#include <assert.h>
#include <stdlib.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "carpet.hh"

extern "C" {
  static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/CarpetParamCheck.cc,v 1.7 2003/07/22 20:09:05 schnetter Exp $";
  CCTK_FILEVERSION(Carpet_Carpet_CarpetParamCheck_cc);
}



namespace Carpet {
  
  using namespace std;
  
  /** Ensure that the parameters have legal values.
   *
   * Note that this checking happens only after most of Carpet has
   * already been set up.
   */
  int CarpetParamCheck (CCTK_ARGUMENTS)
  {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;
    
    if (CCTK_ParameterQueryTimesSet ("periodic", "Carpet")
	|| CCTK_ParameterQueryTimesSet ("periodic_x", "Carpet")
	|| CCTK_ParameterQueryTimesSet ("periodic_y", "Carpet")
	|| CCTK_ParameterQueryTimesSet ("periodic_z", "Carpet")) {
      CCTK_PARAMWARN ("Some of the parameters \"Carpet::periodic*\" have been set.  These parameters are there for compatibility reasons only and must not be used.");
    }
    
    return 0;
  }
  
} // namespace Carpet
