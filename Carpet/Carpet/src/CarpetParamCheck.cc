#include <cassert>
#include <cstdlib>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "carpet.hh"



namespace Carpet {
  
  using namespace std;
  
  /** Ensure that the parameters have legal values.
   *
   * Note that this checking happens only after most of Carpet has
   * already been set up.
   */
  void CarpetParamCheck (CCTK_ARGUMENTS)
  {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;
    
    if (CCTK_ParameterQueryTimesSet ("periodic", "Carpet")
	|| CCTK_ParameterQueryTimesSet ("periodic_x", "Carpet")
	|| CCTK_ParameterQueryTimesSet ("periodic_y", "Carpet")
	|| CCTK_ParameterQueryTimesSet ("periodic_z", "Carpet")) {
      CCTK_PARAMWARN ("Some of the parameters \"Carpet::periodic*\" have been set.  These parameters are there for compatibility reasons only and must not be used.");
    }
  }
  
} // namespace Carpet
