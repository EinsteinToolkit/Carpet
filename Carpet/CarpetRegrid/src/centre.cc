#include <assert.h>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "gh.hh"

#include "carpet.hh"
#include "regrid.hh"

extern "C" {
  static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetRegrid/src/centre.cc,v 1.1 2004/01/25 14:57:30 schnetter Exp $";
  CCTK_FILEVERSION(Carpet_CarpetRegrid_centre_cc);
}



namespace CarpetRegrid {
  
  using namespace std;
  using namespace Carpet;
  
  
  
  int Centre (cGH const * const cctkGH,
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
    DECLARE_CCTK_PARAMETERS;
    
    assert (reflevel>=0 && reflevel<maxreflevels);
    assert (map>=0 && map<maps);
    
    assert (refinement_levels >= 1);
    
    // do nothing if the levels already exist
    if (bbsss.size() == refinement_levels) return 0;
    
    assert (bbsss.size() >= 1);
    
    bbsss.resize (refinement_levels);
    obss.resize (refinement_levels);
    pss.resize (refinement_levels);
    
    bvect const symmetric (symmetry_x, symmetry_y, symmetry_z);
    ivect const zero(0), one(1), two(2);
    
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
      ivect const quarter = (rub - rlb) / 4 / rstr * rstr;
      ivect const half    = (rub - rlb) / 2 / rstr * rstr;
      rlb = oldrlb + symmetric.ifthen(zero, quarter);
      rub = oldrub - symmetric.ifthen(half, quarter);
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
