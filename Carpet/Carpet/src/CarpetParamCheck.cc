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
	or CCTK_ParameterQueryTimesSet ("periodic_x", "Carpet")
	or CCTK_ParameterQueryTimesSet ("periodic_y", "Carpet")
	or CCTK_ParameterQueryTimesSet ("periodic_z", "Carpet")) {
      CCTK_PARAMWARN ("Some of the parameters \"Carpet::periodic*\" have been set.  These parameters are there for compatibility reasons only and must not be used.");
    }
    
    if (adaptive_stepsize and max_refinement_levels > 1) {
      CCTK_PARAMWARN ("Adaptive time step sizes do not work with mesh refinement yet.  Please use only a single level, and set max_refinement_levels=1.");
    }
  }
  
} // namespace Carpet
