#include <assert.h>

#include <vector>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "gh.hh"

#include "carpet.hh"
#include "regrid.hh"

extern "C" {
  static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetRegrid/src/manualcoordinates.cc,v 1.2 2004/02/27 16:24:06 schnetter Exp $";
  CCTK_FILEVERSION(Carpet_CarpetRegrid_manualcoordinates_cc);
}



namespace CarpetRegrid {
  
  using namespace std;
  using namespace Carpet;
  
  
  
  int ManualCoordinates (cGH const * const cctkGH,
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
    
    if (refinement_levels > 4) {
      CCTK_WARN (0, "Cannot currently specify manual refinement regions for more than 4 refinement levels");
    }
    assert (refinement_levels >= 1 && refinement_levels <= 4);
    
    // do nothing if the levels already exist
    if (bbsss.size() == refinement_levels) return 0;
    
    assert (bbsss.size() >= 1);
    
    bbsss.resize (refinement_levels);
    obss.resize (refinement_levels);
    pss.resize (refinement_levels);
    
    vector<rvect> lower(3), upper(3);
    lower.at(0) = rvect (l1xmin, l1ymin, l1zmin);
    upper.at(0) = rvect (l1xmax, l1ymax, l1zmax);
    lower.at(1) = rvect (l2xmin, l2ymin, l2zmin);
    upper.at(1) = rvect (l2xmax, l2ymax, l2zmax);
    lower.at(2) = rvect (l3xmin, l3ymin, l3zmin);
    upper.at(2) = rvect (l3xmax, l3ymax, l3zmax);
    
    for (size_t rl=1; rl<bbsss.size(); ++rl) {
      
      bbvect const ob (false);
      
      vector<ibbox> bbs;
      gh<dim>::cbnds obs;
      
      ManualCoordinates_OneLevel
        (cctkGH, hh, rl, refinement_levels,
         lower.at(rl-1), upper.at(rl-1), ob, bbs, obs);
      
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
  
  
  
  void ManualCoordinates_OneLevel (const cGH * const cctkGH,
                                   const gh<dim> & hh,
                                   const int reflevel,
                                   const int reflevels,
                                   const rvect lower,
                                   const rvect upper,
                                   const bbvect obound,
                                   vector<ibbox> & bbs,
                                   vector<bbvect> & obs)
  {
    if (reflevel >= reflevels) return;
    
    jvect const ilower = pos2int (cctkGH, hh, lower, reflevel);
    jvect const iupper = pos2int (cctkGH, hh, upper, reflevel);
    
    ManualGridpoints_OneLevel
      (cctkGH, hh, reflevel, reflevels, ilower, iupper, obound, bbs, obs);
  }
  
  
  
  ivect delta2int (const cGH * const cctkGH, const gh<dim>& hh,
                   const rvect & rpos, const int reflevel)
  {
    rvect global_lower, global_upper;
    for (int d=0; d<dim; ++d) {
      const int ierr = CCTK_CoordRange
	(cctkGH, &global_lower[d], &global_upper[d], d+1, 0, "cart3d");
      if (ierr<0) {
	global_lower[d] = 0;
	global_upper[d] = 1;
      }
    }
    const ivect global_extent (hh.baseextent.upper() - hh.baseextent.lower());
    
    const rvect scale  = rvect(global_extent) / (global_upper - global_lower);
    const int levfac = ipow(hh.reffact, reflevel);
    assert (all (hh.baseextent.stride() % levfac == 0));
    const ivect istride = hh.baseextent.stride() / levfac;
    
    const ivect ipos
      = ivect(floor(rpos * scale / rvect(istride) + 0.5)) * istride;
    
    const rvect apos = rpos * scale;
    assert (all(abs(apos - rvect(ipos)) < rvect(istride)*0.01));
    
    return ipos;
  }
  
  
  
  ivect pos2int (const cGH* const cctkGH, const gh<dim>& hh,
                 const rvect & rpos, const int reflevel)
  {
    rvect global_lower, global_upper;
    for (int d=0; d<dim; ++d) {
      const int ierr = CCTK_CoordRange
	(cctkGH, &global_lower[d], &global_upper[d], d+1, 0, "cart3d");
      if (ierr<0) {
	global_lower[d] = 0;
	global_upper[d] = 1;
      }
    }
    const ivect global_extent (hh.baseextent.upper() - hh.baseextent.lower());
    
    const rvect scale  = rvect(global_extent) / (global_upper - global_lower);
    const int levfac = ipow(hh.reffact, reflevel);
    assert (all (hh.baseextent.stride() % levfac == 0));
    const ivect istride = hh.baseextent.stride() / levfac;
    
    const ivect ipos
      = (ivect(floor((rpos - global_lower) * scale / rvect(istride) + 0.5))
         * istride);
    
    return ipos;
  }
  
} // namespace CarpetRegrid
