#include <assert.h>

#include <sstream>
#include <vector>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "gh.hh"

#include "carpet.hh"
#include "regrid.hh"

extern "C" {
  static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetRegrid/src/manualgridpoints.cc,v 1.5 2004/07/02 10:14:51 tradke Exp $";
  CCTK_FILEVERSION(Carpet_CarpetRegrid_manualgridpoints_cc);
}



namespace CarpetRegrid {
  
  using namespace std;
  using namespace Carpet;
  
  
  
  int ManualGridpoints (cGH const * const cctkGH,
                        gh<dim> const & hh,
                        gh<dim>::rexts  & bbsss,
                        gh<dim>::rbnds  & obss,
                        gh<dim>::rprocs & pss)
  {
    DECLARE_CCTK_PARAMETERS;
    
    if (refinement_levels > 4) {
      CCTK_WARN (0, "Cannot currently specify manual refinement regions for more than 4 refinement levels");
    }
    assert (refinement_levels >= 1 && refinement_levels <= 4);
    
    // do nothing if the levels already exist
    if (reflevel == refinement_levels) return 0;
    
    assert (bbsss.size() >= 1);
    
    bbsss.resize (refinement_levels);
    obss.resize (refinement_levels);
    pss.resize (refinement_levels);
    
    vector<ivect> ilower(3), iupper(3);
    ilower.at(0) = ivect (l1ixmin, l1iymin, l1izmin);
    iupper.at(0) = ivect (l1ixmax, l1iymax, l1izmax);
    ilower.at(1) = ivect (l2ixmin, l2iymin, l2izmin);
    iupper.at(1) = ivect (l2ixmax, l2iymax, l2izmax);
    ilower.at(2) = ivect (l3ixmin, l3iymin, l3izmin);
    iupper.at(2) = ivect (l3ixmax, l3iymax, l3izmax);
    
    assert (! smart_outer_boundaries);
    
    for (size_t rl=1; rl<bbsss.size(); ++rl) {
      
      bbvect const ob (false);
      
      vector<ibbox> bbs;
      gh<dim>::cbnds obs;
      
      ManualGridpoints_OneLevel
        (cctkGH, hh, rl,refinement_levels,
         ilower.at(rl-1), iupper.at(rl-1), ob, bbs, obs);
      
      // make multiprocessor aware
      gh<dim>::cprocs ps;
      SplitRegions (cctkGH, bbs, obs, ps);
      
      // make multigrid aware
      vector<vector<ibbox> > bbss;
      MakeMultigridBoxes (cctkGH, bbs, obs, bbss);
      
      bbsss.at(rl) = bbss;
      obss.at(rl) = obs;
      pss.at(rl) = ps;
      
    } // for rl
    
    return 1;
  }
  
  
  
  void ManualGridpoints_OneLevel (const cGH * const cctkGH,
                                  const gh<dim> & hh,
                                  const int rl,
                                  const int numrl,
                                  const ivect ilower,
                                  const ivect iupper,
                                  const bbvect obound,
                                  vector<ibbox> & bbs,
                                  vector<bbvect> & obs)
  {
    const ivect rstr = hh.baseextent.stride();
    const ivect rlb  = hh.baseextent.lower();
    const ivect rub  = hh.baseextent.upper();
    
    const int levfac = ipow(hh.reffact, rl);
    assert (all (rstr % levfac == 0));
    const ivect str (rstr / levfac);
    const ivect lb  (ilower);
    const ivect ub  (iupper);
    if (! all(lb>=rlb && ub<=rub)) {
      ostringstream buf;
      buf << "The refinement region boundaries for refinement level #" << rl << " are not within the main grid.  Allowed are the grid point boundaries " << rlb << " - " << rub << "; specified were " << lb << " - " << ub << ends;
      CCTK_WARN (0, buf.str().c_str());
    }
    if (! all(lb<=ub)) {
      ostringstream buf;
      buf << "The refinement region boundaries for refinement level #" << rl << " have the upper boundary (" << ub << ") less than the lower boundary (" << lb << ")" << ends;
      CCTK_WARN (0, buf.str().c_str());
    }
    if (! all(lb%str==0 && ub%str==0)) {
      CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
		  "The refinement region boundaries for refinement level #%d are not a multiple of the stride for that level", rl);
    }
    assert (all(lb>=rlb && ub<=rub));
    assert (all(lb<=ub));
    assert (all(lb%str==0 && ub%str==0));
    
    bbs.push_back (ibbox(lb, ub, str));
    obs.push_back (obound);
  }
  
} // namespace CarpetRegrid
