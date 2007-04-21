#include <cassert>
#include <cmath>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include "mask_surface.hh"



namespace CarpetMask {
  
  using namespace std;
  
  
  
  /**
   * Set the weight in the excluded regions to zero.
   */
  
  void
  CarpetSurfaceSetup (CCTK_ARGUMENTS)
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
      
      int const sn = excluded_surface[n];
      if (sn >= 0) {
        
        assert (sn < nsurfaces);
        
        CCTK_REAL const shrink_factor = excluded_surface_factor[n];
        
        if (sf_valid[sn] > 0) {
          
          CCTK_REAL const x0 = sf_origin_x[sn];
          CCTK_REAL const y0 = sf_origin_y[sn];
          CCTK_REAL const z0 = sf_origin_z[sn];
          
          CCTK_REAL const theta0 = sf_origin_theta[sn];
          CCTK_REAL const phi0   = sf_origin_phi  [sn];
          CCTK_REAL const dtheta = sf_delta_theta[sn];
          CCTK_REAL const dphi   = sf_delta_phi  [sn];
          
          int const ntheta = sf_ntheta[sn];
          int const nphi   = sf_nphi  [sn];
          
          for (int k = 0; k < cctk_lsh[2]; ++ k) {
            for (int j = 0; j < cctk_lsh[1]; ++ j) {
              for (int i = 0; i < cctk_lsh[0]; ++ i) {
                int const ind = CCTK_GFINDEX3D (cctkGH, i, j, k);
                
                CCTK_REAL const dx = x[ind] - x0;
                CCTK_REAL const dy = y[ind] - y0;
                CCTK_REAL const dz = z[ind] - z0;
                
                CCTK_REAL const rho =
                  sqrt (pow (dx, 2) + pow (dy, 2) + pow (dz, 2));
                if (rho < 1.0e-12) {
                  // Always excise the surface origin
                  weight[ind] = 0.0;
                } else {
                  CCTK_REAL const theta = acos (dz / rho);
                  CCTK_REAL const phi   = atan2 (dy, dx);
                  int const a = floor ((theta - theta0) / dtheta + 0.5);
                  int const b = floor ((phi   - phi0  ) / dphi   + 0.5);
                  
                  assert (a >= 0 and a < ntheta);
                  assert (b >= 0 and b < nphi  );
                  
                  CCTK_REAL const dr = sf_radius[a + maxntheta * b];
                  if (rho <= dr * shrink_factor) {
                    weight[ind] = 0.0;
                  }
                }
                
              }
            }
          }
          
        } // if surface is valid
        
      } // if r>=0
    }   // for n
    
  }
  
} // namespace CarpetMask
