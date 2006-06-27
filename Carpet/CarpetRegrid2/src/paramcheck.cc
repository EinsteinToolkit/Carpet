#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "carpet.hh"



namespace CarpetRegrid2 {
  
  using namespace Carpet;
  
  
  
  extern "C" {
    void CarpetRegrid2_ParamCheck (CCTK_ARGUMENTS);
  }
  
  
  
  void CarpetRegrid2_ParamCheck (CCTK_ARGUMENTS)
  {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;
    
    if (num_centres >= 1 and num_levels_1 > maxreflevels or
        num_centres >= 2 and num_levels_2 > maxreflevels or
        num_centres >= 3 and num_levels_3 > maxreflevels)
    {
      CCTK_PARAMWARN ("The number of requested refinement levels is larger than the maximum number of levels specified by Carpet::max_refinement_levels");
    }
  }
  
} // namespace CarpetRegrid2
