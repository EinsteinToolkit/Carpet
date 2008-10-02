#include <cmath>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include "mask_excluded.hh"



namespace CarpetMask {
  
  using namespace std;
  
  
  
  /**
   * Set the weight in the excluded regions to zero.
   */
  
  void
  CarpetExcludedSetup (CCTK_ARGUMENTS)
  {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;
    
    CCTK_REAL * const weight =
      static_cast <CCTK_REAL *>
      (CCTK_VarDataPtr (cctkGH, 0, "CarpetReduce::weight"));
    
    if (not weight) {
      CCTK_WARN (CCTK_WARN_ABORT,
                 "CarpetReduce is not active, or CarpetReduce::mask does not have storage");
    }
    
    for (int n = 0; n < 10; ++ n) {
      
      CCTK_REAL const r0 = excluded_radius[n];
      if (r0 >= 0.0) {
        
        CCTK_REAL const x0 = excluded_centre_x[n];
        CCTK_REAL const y0 = excluded_centre_y[n];
        CCTK_REAL const z0 = excluded_centre_z[n];
        
        CCTK_REAL const r2 = pow (r0, 2);
        
        bool const exterior = exclude_exterior[n];
        
        for (int k = 0; k < cctk_lsh[2]; ++ k) {
          for (int j = 0; j < cctk_lsh[1]; ++ j) {
            for (int i = 0; i < cctk_lsh[0]; ++ i) {
              int const ind = CCTK_GFINDEX3D (cctkGH, i, j, k);
              
              CCTK_REAL const dx2 = pow (x[ind] - x0, 2);
              CCTK_REAL const dy2 = pow (y[ind] - y0, 2);
              CCTK_REAL const dz2 = pow (z[ind] - z0, 2);
              
              if (exterior ?
                  dx2 + dy2 + dz2 >= r2 :
                  dx2 + dy2 + dz2 <= r2)
              {
                weight[ind] = 0.0;
              }
              
            }
          }
        }
        
      } // if r>=0
    }   // for n
    
  }
  
} // namespace CarpetMask
