#include <assert.h>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "gh.hh"

#include "carpet.hh"
#include "regrid.hh"

extern "C" {
  static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetRegrid/src/moving.cc,v 1.1 2004/04/14 22:19:44 schnetter Exp $";
  CCTK_FILEVERSION(Carpet_CarpetRegrid_moving_cc);
}



namespace CarpetRegrid {
  
  using namespace std;
  using namespace Carpet;
  
  
  
  int Moving (cGH const * const cctkGH,
              gh<dim> const & hh,
              int const reflevel,
              int const map,
              int const size,
              jjvect const & nboundaryzones,
              jjvect const & is_internal,
              jjvect const & is_staggered,
              jjvect const & shiftout,
              gh<dim>::rexts  & bbsss,
              gh<dim>::rbnds  & obss,
              gh<dim>::rprocs & pss)
  {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;
    
    assert (reflevel>=0 && reflevel<maxreflevels);
    assert (map>=0 && map<maps);
    
    assert (refinement_levels >= 1);
    
    // do nothing if the levels already exist
    if (reflevel == refinement_levels) return 0;
    
    assert (bbsss.size() >= 1);
    
    bbsss.resize (refinement_levels);
    obss.resize (refinement_levels);
    pss.resize (refinement_levels);
    
    ivect rstr = hh.baseextent.stride();
    ivect rlb  = hh.baseextent.lower();
    ivect rub  = hh.baseextent.upper();
    
    for (size_t rl=1; rl<bbsss.size(); ++rl) {
      
      // save old values
      ivect const oldrlb = rlb;
      ivect const oldrub = rub;
      
      // refined boxes have smaller stride
      assert (all(rstr%hh.reffact == 0));
      rstr /= hh.reffact;
      
      // calculate new extent
      rvect pos;
      pos[0] = moving_centre_x + moving_circle_radius * cos(2*M_PI * moving_circle_frequency * cctk_time);
      pos[1] = moving_centre_y + moving_circle_radius * sin(2*M_PI * moving_circle_frequency * cctk_time);
      pos[2] = moving_centre_z;
      ivect const centre = floor(rvect(rub - rlb) * pos / rstr + 0.5) * rstr;
      ivect const radius = floor(rvect(rub - rlb) * moving_region_radius / rstr + 0.5) * rstr;
      rlb = oldrlb + centre - radius;
      rub = oldrlb + centre + radius;
      assert (all(rlb >= oldrlb && rub <= oldrub));
      
      ibbox const bb (rlb, rub, rstr);
      vector<ibbox> bbs (1);
      bbs.at(0) = bb;
      
      bbvect const ob (false);
      gh<dim>::cbnds obs (1);
      obs.at(0) = ob;
      
      // make multiprocessor aware
      gh<dim>::cprocs ps;
      SplitRegions (cctkGH, bbs, obs, ps);
      
      // make multigrid aware
      vector<vector<ibbox> > bbss;
      MakeMultigridBoxes
        (cctkGH,
         size, nboundaryzones, is_internal, is_staggered, shiftout,
         bbs, obs, bbss);
      
      bbsss.at(rl) = bbss;
      obss.at(rl) = obs;
      pss.at(rl) = ps;
      
    } // for rl
    
    return 1;
  }
  
} // namespace CarpetRegrid
