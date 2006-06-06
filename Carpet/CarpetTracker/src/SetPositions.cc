#include <cassert>
#include <sstream>
#include <string>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"



namespace CarpetTracker {

using namespace std;
  
  
  
  extern "C" {
    void
    CarpetTracker_SetPositions (CCTK_ARGUMENTS);
  }
  
  
  
  void
  SetParameter (char const * const dir, int const n, CCTK_REAL const pos)
  {
    ostringstream name_buf;
    name_buf << "position_" << dir << "_" << n + 1;
    string const name (name_buf.str());
    ostringstream pos_buf;
    pos_buf << pos;
    string const pos_str (pos_buf.str());
    int const ierr =
      CCTK_ParameterSet (name.c_str(), "CarpetRegrid2", pos_str.c_str());
    assert (not ierr);
  }
  
  
  
  void
  CarpetTracker_SetPositions (CCTK_ARGUMENTS)
  {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;
    
    for (int n = 0; n < 3; ++ n) {
      int const sn = surface[n];
      if (sn >= 0) {
        assert (sn >= 0 and sn < nsurfaces);
        
        if (sf_valid[sn] > 0) {
          
          if (verbose) {
            CCTK_VInfo (CCTK_THORNSTRING,
                        "Setting position of refined region #%d from surface #%d to (%g,%g,%g)",
                        n + 1, sn,
                        static_cast <double> (sf_origin_x[sn]),
                        static_cast <double> (sf_origin_y[sn]),
                        static_cast <double> (sf_origin_z[sn]));
          }
          SetParameter ("x", n, sf_origin_x[sn]);
          SetParameter ("y", n, sf_origin_y[sn]);
          SetParameter ("z", n, sf_origin_z[sn]);
          
        } else {
          
          if (verbose) {
            CCTK_VInfo (CCTK_THORNSTRING,
                        "No position information available for refined region #%d from surface #%d",
                        n + 1, sn);
          }
          
        }
        
      } // if
    } // for
  }
  
} // namespace CarpetTracker
