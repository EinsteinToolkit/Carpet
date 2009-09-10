#include <cassert>
#include <sstream>
#include <string>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"



namespace CarpetTracker {

using namespace std;
  
  
  
  // Maximum number of tracked surfaces
  int const num_surfaces = 10;
  
  
  
  extern "C" {
    void
    CarpetTracker_SetPositions (CCTK_ARGUMENTS);
  }
  
  
  
  void
  CarpetTracker_SetPositions (CCTK_ARGUMENTS)
  {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;
    
    static CCTK_INT cctk_iteration_done = -1;

    if (cctk_iteration == cctk_iteration_done) return;
    cctk_iteration_done = cctk_iteration;

    for (int n = 0; n < num_surfaces; ++ n) {
      int const sn = surface[n];
      if (sn >= 0) {
        assert (sn >= 0 and sn < nsurfaces);
        
        if (sf_active[sn]) {
          if (sf_valid[sn] > 0) {
            
            if (verbose) {
              CCTK_VInfo (CCTK_THORNSTRING,
                          "Setting position of refined region #%d from surface #%d from (%g,%g,%g) to (%g,%g,%g)",
                          n + 1, sn,
                          static_cast <double> (position_x[n]),
                          static_cast <double> (position_y[n]),
                          static_cast <double> (position_z[n]),
                          static_cast <double> (sf_centroid_x[sn]),
                          static_cast <double> (sf_centroid_y[sn]),
                          static_cast <double> (sf_centroid_z[sn]));
            }
            
            // Activate region
            active[n] = 1;
            
            // Set position in CarpetRegrid2
            position_x[n] = sf_centroid_x[sn];
            position_y[n] = sf_centroid_y[sn];
            position_z[n] = sf_centroid_z[sn];
            
          } else {
            
            if (verbose) {
              CCTK_VInfo (CCTK_THORNSTRING,
                          "No position information available for refined region #%d from surface #%d",
                          n + 1, sn);
            }
            
          } // if not valid
          
        } else {
          
          if (verbose) {
            CCTK_VInfo (CCTK_THORNSTRING,
                        "Refined region #%d (depending on surface #%d) is inactive",
                        n + 1, sn);
          }
          
          // Deactivate region
          active[n] = 0;
          
        } // if not active
        
      } // if sn > 0
    } // for n
  }
  
} // namespace CarpetTracker
