#include <assert.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "carpet.hh"
#include "regrid.hh"



namespace CarpetRegrid {
  
  using namespace std;
  using namespace Carpet;
  
  void CarpetRegridParamcheck (CCTK_ARGUMENTS)
  {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;
    
    if (refinement_levels > maxreflevels) {
      CCTK_PARAMWARN ("The parameter CarpetRegrid::refinement_levels is larger than Carpet::max_refinement_levels");
    }

    if (smart_outer_boundaries) {
      int type;
      const CCTK_INT * const domain_from_coordbase
        = (const CCTK_INT *) CCTK_ParameterGet ("domain_from_coordbase", "Carpet", &type);
      assert (domain_from_coordbase);
      assert (type == PARAMETER_BOOLEAN);
      if (! *domain_from_coordbase) {
        CCTK_PARAMWARN ("The parameter CarpetRegrid::smart_outer_boundaries can only be used when Carpet::domain_from_coordbase=yes");
      }
      if (CCTK_Equals(refined_regions, "manual-coordinate-list")) {
        // do nothing
      } else {
        CCTK_PARAMWARN ("The parameter CarpetRegrid::smart_outer_boundaries can currently only be used when CarpetRegrid::refined_regions is set to \"manual-coordinate-list\"");
      }
    }
  }
  
}
