#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <carpet.hh>



namespace CarpetRegrid2 {
  
  using namespace Carpet;
  
  
  
  extern "C" {
    void
    CarpetRegrid2_ParamCheck (CCTK_ARGUMENTS);
  }
  
  
  
  void
  CarpetRegrid2_ParamCheck (CCTK_ARGUMENTS)
  {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;
    
    if ((num_centres >= 1 and num_levels_1 > maxreflevels) or
        (num_centres >= 2 and num_levels_2 > maxreflevels) or
        (num_centres >= 3 and num_levels_3 > maxreflevels) or
        (num_centres >= 4 and num_levels_4 > maxreflevels) or
        (num_centres >= 5 and num_levels_5 > maxreflevels) or
        (num_centres >= 6 and num_levels_6 > maxreflevels) or
        (num_centres >= 7 and num_levels_7 > maxreflevels) or
        (num_centres >= 8 and num_levels_8 > maxreflevels) or
        (num_centres >= 9 and num_levels_9 > maxreflevels) or
        (num_centres >= 10 and num_levels_10 > maxreflevels))
    {
      CCTK_PARAMWARN ("The number of requested refinement levels is larger than the maximum number of levels specified by Carpet::max_refinement_levels");
    }
    
    if (symmetry_rotating90) {
      CCTK_PARAMWARN ("symmetry_rotating90 is not yet implemented");
    }
  }
  
} // namespace CarpetRegrid2
