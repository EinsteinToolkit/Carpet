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
    
    if (num_centres >= 4) {
      num_levels[3] = num_levels_4;
      position_x[3] = position_x_4;
      position_y[3] = position_y_4;
      position_z[3] = position_z_4;
      for (int rl = 0; rl < num_levels[3]; ++ rl) {
        radius[index2 (lsh, rl, 3)] = radius_4[rl];
      }
    }
    
    if (num_centres >= 5) {
      num_levels[4] = num_levels_5;
      position_x[4] = position_x_5;
      position_y[4] = position_y_5;
      position_z[4] = position_z_5;
      for (int rl = 0; rl < num_levels[4]; ++ rl) {
        radius[index2 (lsh, rl, 4)] = radius_5[rl];
      }
    }
    
    if (num_centres >= 6) {
      num_levels[5] = num_levels_6;
      position_x[5] = position_x_6;
      position_y[5] = position_y_6;
      position_z[5] = position_z_6;
      for (int rl = 0; rl < num_levels[5]; ++ rl) {
        radius[index2 (lsh, rl, 5)] = radius_6[rl];
      }
    }
    
    if (num_centres >= 7) {
      num_levels[6] = num_levels_7;
      position_x[6] = position_x_7;
      position_y[6] = position_y_7;
      position_z[6] = position_z_7;
      for (int rl = 0; rl < num_levels[6]; ++ rl) {
        radius[index2 (lsh, rl, 6)] = radius_7[rl];
      }
    }
    
    if (num_centres >= 8) {
      num_levels[7] = num_levels_8;
      position_x[7] = position_x_8;
      position_y[7] = position_y_8;
      position_z[7] = position_z_8;
      for (int rl = 0; rl < num_levels[7]; ++ rl) {
        radius[index2 (lsh, rl, 7)] = radius_8[rl];
      }
    }
    
    if (num_centres >= 9) {
      num_levels[8] = num_levels_9;
      position_x[8] = position_x_9;
      position_y[8] = position_y_9;
      position_z[8] = position_z_9;
      for (int rl = 0; rl < num_levels[8]; ++ rl) {
        radius[index2 (lsh, rl, 8)] = radius_9[rl];
      }
    }
    
    if (num_centres >= 10) {
      num_levels[9] = num_levels_10;
      position_x[9] = position_x_10;
      position_y[9] = position_y_10;
      position_z[9] = position_z_10;
      for (int rl = 0; rl < num_levels[9]; ++ rl) {
        radius[index2 (lsh, rl, 9)] = radius_10[rl];
      }
    }
    
    for (int n = 0; n < num_centres; ++ n) {
      old_position_x[n] = position_x[n];
      old_position_y[n] = position_y[n];
      old_position_z[n] = position_z[n];
    }
  }
  
} // namespace CarpetRegrid2
