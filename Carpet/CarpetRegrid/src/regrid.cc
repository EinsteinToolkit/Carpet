#include <assert.h>
#include <math.h>
#include <stdarg.h>
#include <stdlib.h>

#include <list>
#include <vector>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "Carpet/CarpetLib/src/bbox.hh"
#include "Carpet/CarpetLib/src/bboxset.hh"
#include "Carpet/CarpetLib/src/gh.hh"
#include "Carpet/CarpetLib/src/gf.hh"
#include "Carpet/CarpetLib/src/vect.hh"

#include "carpet.hh"
#include "regrid.hh"

static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetRegrid/src/regrid.cc,v 1.1 2001/12/14 16:34:39 schnetter Exp $";



// That's a hack
namespace Carpet {
  void Waypoint (const char* fmt, ...);
}

namespace CarpetRegrid {
  
  using namespace std;
  using namespace Carpet;
  
  
  
  int CarpetRegridRegrid (CCTK_ARGUMENTS)
  {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;
    
    // Return if no regridding is desired
    if (regrid_every == -1) return 0;
    
    // Return if this is the finest possible level
    if (reflevel == maxreflevels-1) return 0;
    
    // Return if we want to regrid during initial data only, and this
    // is not the time for initial data
    if (regrid_every == 0 && cctk_iteration != 0) return 0;
    
    // Return if we want to regrid regularly, but not at this time
    if (regrid_every > 0 && cctk_iteration % regrid_every != 0) return 0;
    
    Waypoint ("%*sRegrid", 2*reflevel, "");
    
    list<bbox<int,dim> > bbl;
    
    if (CCTK_EQUALS(refined_regions, "none")) {
      
      MakeRegions_BaseLevel (cctkGH, bbl);
      
    } else if (CCTK_EQUALS(refined_regions, "centre")) {
      
      MakeRegions_RefineCentre (cctkGH, refinement_levels, bbl);
      
    } else if (CCTK_EQUALS(refined_regions, "manual")) {
      
      if (refinement_levels > 4) {
	CCTK_WARN (0, "Cannot currently specify manual refinement regions for more than 4 refinement levels");
      }
      assert (refinement_levels<4);
      vector<vect<int,dim> > lower(3), upper(3);
      lower[0] = vect<int,dim> (l1xmin, l1ymin, l1zmin);
      upper[0] = vect<int,dim> (l1xmax, l1ymax, l1zmax);
      lower[1] = vect<int,dim> (l2xmin, l2ymin, l2zmin);
      upper[1] = vect<int,dim> (l2xmax, l2ymax, l2zmax);
      lower[2] = vect<int,dim> (l3xmin, l3ymin, l3zmin);
      upper[2] = vect<int,dim> (l3xmax, l3ymax, l3zmax);
      MakeRegions_AsSpecified (cctkGH, refinement_levels, lower, upper, bbl);
      
    } else if (CCTK_EQUALS(refined_regions, "automatic")) {
      
      const int vi = CCTK_VarIndex (errorvar);
      assert (vi>=0 && vi<CCTK_NumVars());
      const int gi = CCTK_GroupIndexFromVarI (vi);
      assert (gi>=0 && gi<CCTK_NumGroups());
      
      assert (CCTK_VarTypeI(vi) == CCTK_VARIABLE_REAL);
      
      const gf<CCTK_REAL,dim>& error
	= *(const gf<CCTK_REAL,dim>*)arrdata[gi].data[vi];
      
      MakeRegions_Adaptively (cctkGH, minwidth, minfraction, maxerror, error,
			      bbl);
      
    } else {
      abort();
    }
    
    vector<bbox<int,dim> > bbs;
    bbs.reserve (bbl.size());
    for (list<bbox<int,dim> >::const_iterator it = bbl.begin(); it != bbl.end(); ++it) {
      bbs.push_back (*it);
    }
    
    gh<dim>::cexts bbss;
    const int mglevels = 1;	// arbitrary value
    bbss = hh->make_reflevel_multigrid_boxes(bbs, mglevels);
    
    // TODO: ensure nesting properties
    
    SplitRegions (cctkGH, bbss);
    
    gh<dim>::rexts bbsss = hh->extents;
    if (bbss.size() == 0) {
      // TODO: this might not work
      bbsss.resize(reflevel+1);
    } else {
      assert (reflevel < (int)bbsss.size());
      if (reflevel+1 == (int)bbsss.size()) {
	bbsss.push_back (bbss);
      } else {
	bbsss[reflevel+1] = bbss;
      }
    }
    
    gh<dim>::rprocs pss;
    MakeProcessors (cctkGH, bbsss, pss);
    
    RegisterRecomposeRegions (bbsss, pss);
    
    return 0;
  }
  
  
  
  void MakeRegions_BaseLevel (const cGH* cctkGH, list<bbox<int,dim> >& bbl)
  {
    assert (bbl.empty());
  }
  
  
  
  // This is a helpful helper routine. The user can use it to define
  // how the hierarchy should be refined. But the result of this
  // routine is rather arbitrary.
  void MakeRegions_RefineCentre (const cGH* cctkGH, const int reflevels,
				 list<bbox<int,dim> >& bbl)
  {
    assert (bbl.empty());
    
    if (reflevel+1 >= reflevels) return;
      
    // note: what this routine calls "ub" is "ub+str" elsewhere
    vect<int,dim> rstr = hh->baseextent.stride();
    vect<int,dim> rlb  = hh->baseextent.lower();
    vect<int,dim> rub  = hh->baseextent.upper() + rstr;
    
    for (int rl=0; rl<reflevel+1; ++rl) {
      // save old values
      const vect<int,dim> oldrlb = rlb;
      const vect<int,dim> oldrub = rub;
      // calculate extent and centre
      const vect<int,dim> rextent = rub - rlb;
      const vect<int,dim> rcentre = rlb + (rextent / 2 / rstr) * rstr;
      // calculate new extent
      assert (all(rextent % hh->reffact == 0));
      const vect<int,dim> newrextent = rextent / hh->reffact;
      // refined boxes have smaller stride
      assert (all(rstr%hh->reffact == 0));
      rstr /= hh->reffact;
      // refine (arbitrarily) around the center only
      rlb = rcentre - (newrextent/2 / rstr) * rstr;
      // // refine (arbitrarily) around the lower boundary only
      // rlb = rlb;
      rub = rlb + newrextent;
      // require rub<oldrub because we really want rub-rstr<=oldrub-oldstr
      assert (all(rlb >= oldrlb && rub < oldrub));
    }
    
    bbl.push_back (bbox<int,dim>(rlb, rub-rstr, rstr));
  }
  
  
  
  void MakeRegions_AsSpecified (const cGH* cctkGH, const int reflevels,
				const vector<vect<int,dim> > lower,
				const vector<vect<int,dim> > upper,
				list<bbox<int,dim> >& bbl)
  {
    assert (bbl.empty());
    
    if (reflevel+1 >= reflevels) return;
    
    const vect<int,dim> rstr = hh->baseextent.stride();
    const vect<int,dim> rlb  = hh->baseextent.lower();
    const vect<int,dim> rub  = hh->baseextent.upper();
    
    const int rl = reflevel+1;
    
    const int levfac = floor(pow((double)hh->reffact, rl) + 0.5);
    assert (all (rstr % levfac == 0));
    const vect<int,dim> str (rstr / levfac);
    const vect<int,dim> lb  (lower[rl-1]);
    const vect<int,dim> ub  (upper[rl-1]);
    if (! all(lb>=rlb && ub<=rub)) {
      CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
		  "The refinement region boundaries for refinement level #%d are not within the main grid", rl);
    }
    if (! all(lb<=ub)) {
      CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
		  "The refinement region boundaries for refinement level #%d have the upper boundary less than the lower boundary", rl);
    }
    if (! all(lb%str==0 && ub%str==0)) {
      CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
		  "The refinement region boundaries for refinement level #%d are not a multiple of the stride for that level", rl);
    }
    assert (all(lb>=rlb && ub<=rub));
    assert (all(lb<=ub));
    assert (all(lb%str==0 && ub%str==0));
    
    bbl.push_back (bbox<int,dim>(lb, ub, str));
  }
  
  
  
  static void
  MakeRegions_Adaptively (const cGH* const cctkGH,
			  const int minwidth, const double minfraction,
			  const CCTK_REAL maxerror,
			  const data<CCTK_REAL,dim>* const error,
			  list<bbox<int,dim> >& bbl,
			  const bbox<int,dim> region)
  {
    // Just to be sure
    assert (! region.empty());
    
    // Count grid points that need to be refined
    int cnt = 0;
    for (bbox<int,dim>::iterator it=region.begin(); it!=region.end(); ++it) {
      if ((*error)[*it] > maxerror) ++cnt;
    }
    const CCTK_REAL fraction = (CCTK_REAL)cnt / region.num_points();
    const int width = minval(region.shape());
    
    if (cnt == 0) {
      // Don't refine
    } else if (width < 2*minwidth || fraction >= minfraction) {
      // Refine the whole region
      bbl.push_back (region);
    } else {
      // Split the region and check recursively
      const int d = maxloc(region.shape());
      const vect<int,dim> lo(region.lower());
      const vect<int,dim> up(region.upper());
      const vect<int,dim> str(region.stride());
      vect<int,dim> lo1(lo), lo2(lo);
      vect<int,dim> up1(up), up2(up);
      up1[d] = (((up[d]-lo[d])/str[d]+1)/2-1)*str[d];
      lo2[d] = up1[d]+str[d];
      const bbox<int,dim> region1(lo1,up1,str);
      const bbox<int,dim> region2(lo2,up2,str);
      assert (region1 <= region);
      assert (region2 <= region);
      assert ((region1 & region2).empty());
      assert (region1 + region2 == region);
      MakeRegions_Adaptively (cctkGH, minwidth, minfraction, maxerror,
			      error, bbl, region1);
      MakeRegions_Adaptively (cctkGH, minwidth, minfraction, maxerror,
			      error, bbl, region2);
    }
    
  }
  
  void
  MakeRegions_Adaptively (const cGH* const cctkGH,
			  const int minwidth, const double minfraction,
			  const CCTK_REAL maxerror,
			  const gf<CCTK_REAL,dim>& error,
			  list<bbox<int,dim> >& bbl)
  {
    assert (bbl.empty());
    
    if (reflevel+1 >= maxreflevels) return;
    
    assert (component == -1);
    
    const int rl = reflevel;
    
    // Arbitrary
    const int tl = 0;
    const int ml = 0;
    
    for (int c=0; c<hh->components(rl); ++c) {
      const bbox<int,dim> region = hh->extents[rl][c][ml];
      assert (! region.empty());
      
      const data<CCTK_REAL,dim>* errdata = error(tl,rl,c,ml);
      
      MakeRegions_Adaptively (cctkGH, minwidth, minfraction, maxerror,
			      errdata, bbl, region);
    }
  }
  
} // namespace CarpetRegrid
