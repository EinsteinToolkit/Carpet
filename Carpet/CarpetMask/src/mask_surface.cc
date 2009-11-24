#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include "mask_surface.hh"
#include "loopcontrol.h"



namespace CarpetMask {
  
  using namespace std;
  
  
  
  // Number of excluded regions which are defined in param.ccl
  int const num_excluded = 10;
  
  
  
  /**
   * Set the weight in the excluded regions to zero.
   */
  
  void
  CarpetSurfaceSetup (CCTK_ARGUMENTS)
  {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;
    
    // Some state verbose output
    static vector <int> last_output;
    if (last_output.empty()) {
      last_output.resize (num_excluded, -1);
    }
    
    CCTK_REAL * const weight =
      static_cast <CCTK_REAL *>
      (CCTK_VarDataPtr (cctkGH, 0, "CarpetReduce::weight"));
    
    if (not weight) {
      CCTK_WARN (CCTK_WARN_ABORT,
                 "CarpetReduce is not active, or CarpetReduce::mask does not have storage");
    }
    
    for (int n = 0; n < num_excluded; ++ n) {
      
      int const sn = excluded_surface[n];
      if (sn >= 0) {
        
        assert (sn < nsurfaces);
        
        if (sf_active[sn]) {
          
          CCTK_REAL const shrink_factor = excluded_surface_factor[n];
          
          CCTK_REAL const x0 = sf_origin_x[sn];
          CCTK_REAL const y0 = sf_origin_y[sn];
          CCTK_REAL const z0 = sf_origin_z[sn];
          
          CCTK_REAL const theta0 = sf_origin_theta[sn];
          CCTK_REAL const phi0   = sf_origin_phi  [sn];
          CCTK_REAL const dtheta = sf_delta_theta[sn];
          CCTK_REAL const dphi   = sf_delta_phi  [sn];
          
          int const ntheta = sf_ntheta[sn];
          int const nphi   = sf_nphi  [sn];
          
          if (verbose) {
            if (cctk_iteration > last_output.at(n)) {
              last_output.at(n) = cctk_iteration;
              CCTK_VInfo (CCTK_THORNSTRING,
                          "Setting mask from surface %d at [%g,%g,%g], average radius %g",
                          sn,
                          double(x0), double(y0), double(z0),
                          double(sf_mean_radius[sn]));
            }
          }
          
              #pragma omp parallel
              LC_LOOP3 (CarpetSurfaceSetup,
                      i,j,k, 0,0,0, cctk_lsh[0], cctk_lsh[1], cctk_lsh[2],
                      cctk_lssh[CCTK_LSSH_IDX(0,0)],cctk_lssh[CCTK_LSSH_IDX(0,1)],cctk_lssh[CCTK_LSSH_IDX(0,2)])
              {
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
                  CCTK_REAL theta =
                    acos (min (CCTK_REAL (+1.0),
                               max (CCTK_REAL (-1.0), dz / rho)));
                  if (symmetric_z[sn]) {
                    if (theta > M_PI/2.0) {
                      theta = M_PI - theta;
                    }
                  }
                  
                  assert (not isnan (theta));
                  assert (theta >= 0);
                  assert (theta <= M_PI);
                  CCTK_REAL phi =
                    fmod (atan2 (dy, dx) + CCTK_REAL (2 * M_PI),
                          CCTK_REAL (2 * M_PI));
                  if (symmetric_x[sn] or symmetric_y[sn]) {
                    if (symmetric_x[sn] and symmetric_y[sn]) {
                      if (phi > M_PI / 2.0) {
                        phi = M_PI - phi;
                      }
                    } else {
                      if (phi > M_PI) {
                        phi = 2 * M_PI - phi;
                      }
                    }
                  }
                  assert (not isnan (phi));
                  assert (phi >= 0);
                  assert (phi < 2 * M_PI);
                  int const a = floor ((theta - theta0) / dtheta + 0.5);
                  assert (a >= 0);
                  assert (a < ntheta);
                  int const b = floor ((phi   - phi0  ) / dphi   + 0.5);
                  assert (b >= 0);
                  assert (b < nphi);
                  
                  assert (a >= 0 and a < ntheta);
                  assert (b >= 0 and b < nphi  );
                  
                  CCTK_REAL const dr =
                    sf_radius[a + maxntheta * (b + maxnphi * sn)];
                  if (rho <= dr * shrink_factor) {
                    weight[ind] = 0.0;
                  }
                }
                
              }
              LC_ENDLOOP3 (CarpetSurfaceSetup);
          
        } else {
          
          if (verbose) {
            if (cctk_iteration > last_output.at(n)) {
              last_output.at(n) = cctk_iteration;
              CCTK_VInfo (CCTK_THORNSTRING, "Surface %d is not active", sn);
            }
          }
          
        } // if surface is not active
        
      } // if sn >= 0
    }   // for n
    
  }
  
} // namespace CarpetMask
