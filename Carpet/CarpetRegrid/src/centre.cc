#include <assert.h>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "gh.hh"

#include "carpet.hh"
#include "regrid.hh"



namespace CarpetRegrid {
  
  using namespace std;
  using namespace Carpet;
  
  
  
  int Centre (cGH const * const cctkGH,
              gh const & hh,
              gh::mexts  & bbsss,
              gh::rbnds  & obss,
              gh::rprocs & pss)
  {
    DECLARE_CCTK_PARAMETERS;
    
    assert (refinement_levels >= 1);
    
    // do nothing if the levels already exist
    if (reflevel == refinement_levels) return 0;
    
    assert (bbsss.size() >= 1);
    vector<vector<ibbox> > bbss = bbsss.at(0);
    
    bbss.resize (refinement_levels);
    obss.resize (refinement_levels);
    pss.resize (refinement_levels);
    
    bvect const symmetric (symmetry_x, symmetry_y, symmetry_z);
    ivect const zero(0), one(1), two(2);
    
    ivect rstr = hh.baseextent.stride();
    ivect rlb  = hh.baseextent.lower();
    ivect rub  = hh.baseextent.upper();
    
    assert (! smart_outer_boundaries);
    
    for (size_t rl=1; rl<bbss.size(); ++rl) {
      
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
      assert (all(rlb >= oldrlb and rub <= oldrub));
      
      ibbox const bb (rlb, rub, rstr);
      vector<ibbox> bbs (1);
      bbs.at(0) = bb;
      
      bbvect const ob (false);
      gh::cbnds obs (1);
      obs.at(0) = ob;
      
      // make multiprocessor aware
      gh::cprocs ps;
      SplitRegions (cctkGH, bbs, obs, ps);
      
      bbss.at(rl) = bbs;
      obss.at(rl) = obs;
      pss.at(rl) = ps;
      
      // make multigrid aware
      MakeMultigridBoxes (cctkGH, bbss, obss, bbsss);
      
    } // for rl
    
    return 1;
  }
  
} // namespace CarpetRegrid
