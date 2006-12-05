#include <cassert>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include "indexing.hh"



namespace CarpetRegrid2 {
  
  extern "C" {
    void
    CarpetRegrid2_Initialise (CCTK_ARGUMENTS);
  }
  
  
  
  void
  CarpetRegrid2_Initialise (CCTK_ARGUMENTS)
  {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;
    
    // Initialise meta-information
    * last_iteration = -1;
    * last_map = -1;
    
    // Initialise refinement information
    num_levels[0] = 0;
    num_levels[1] = 0;
    num_levels[2] = 0;
    
    int lsh[2];
    getvectorindex2 (cctkGH, "CarpetRegrid2::radii", lsh);
    
    if (num_centres >= 1) {
      num_levels[0] = num_levels_1;
      position_x[0] = position_x_1;
      position_y[0] = position_y_1;
      position_z[0] = position_z_1;
      for (int rl = 0; rl < num_levels[0]; ++ rl) {
        radius[index2 (lsh, rl, 0)] = radius_1[rl];
      }
    }
    
    if (num_centres >= 2) {
      num_levels[1] = num_levels_2;
      position_x[1] = position_x_2;
      position_y[1] = position_y_2;
      position_z[1] = position_z_2;
      for (int rl = 0; rl < num_levels[1]; ++ rl) {
        radius[index2 (lsh, rl, 1)] = radius_2[rl];
      }
    }
    
    if (num_centres >= 3) {
      num_levels[2] = num_levels_3;
      position_x[2] = position_x_3;
      position_y[2] = position_y_3;
      position_z[2] = position_z_3;
      for (int rl = 0; rl < num_levels[2]; ++ rl) {
        radius[index2 (lsh, rl, 2)] = radius_3[rl];
      }
    }
  }
  
} // namespace CarpetRegrid2
