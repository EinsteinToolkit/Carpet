#include <assert.h>
#include <math.h>
#include <stdarg.h>
#include <stdlib.h>

#include <list>
#include <vector>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "Carpet/CarpetLib/src/bbox.hh"
#include "Carpet/CarpetLib/src/bboxset.hh"
#include "Carpet/CarpetLib/src/defs.hh"
#include "Carpet/CarpetLib/src/gh.hh"
#include "Carpet/CarpetLib/src/gf.hh"
#include "Carpet/CarpetLib/src/vect.hh"

#include "carpet.hh"
#include "regrid.hh"

static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetRegrid/src/regrid.cc,v 1.8 2002/03/06 17:47:58 schnetter Exp $";



namespace CarpetRegrid {
  
  using namespace std;
  using namespace Carpet;
  
  
  
  int CarpetRegridStartup ()
  {
    RegisterRegridRoutine (CarpetRegridRegrid);
    return 0;
  }
  
  
  
  int CarpetRegridRegrid (const cGH * const cctkGH,
			  gh<dim>::rexts& bbsss,
			  gh<dim>::rprocs& pss)
  {
    DECLARE_CCTK_PARAMETERS;
    
    assert (regrid_every == -1 || regrid_every == 0
	    || regrid_every % maxmglevelfact == 0);
    
    // Return if no regridding is desired
    if (regrid_every == -1) return 0;
    
#if 0
    // Return if this is not the main hierarchy
    if (mglevel != 0) return 0;
#endif
    assert (mglevel == -1);
    
    // Return if this is the finest possible level
    if (reflevel == maxreflevels-1) return 0;
    
    // Return if we want to regrid during initial data only, and this
    // is not the time for initial data
    if (regrid_every == 0 && cctkGH->cctk_iteration != 0) return 0;
    
    // Return if we want to regrid regularly, but not at this time
    if (regrid_every > 0 && cctkGH->cctk_iteration != 0
	&& cctkGH->cctk_iteration % regrid_every != 0) {
      return 0;
    }
    
    list<bbox<int,dim> > bbl;
    
    if (CCTK_EQUALS(refined_regions, "none")) {
      
      MakeRegions_BaseLevel (cctkGH, bbl);
      
    } else if (CCTK_EQUALS(refined_regions, "centre")) {
      
      MakeRegions_RefineCentre (cctkGH, refinement_levels, bbl);
      
    } else if (CCTK_EQUALS(refined_regions, "manual-gridpoints")) {
      
      if (refinement_levels > 4) {
	CCTK_WARN (0, "Cannot currently specify manual refinement regions for more than 4 refinement levels");
      }
      assert (refinement_levels<=4);
      vector<vect<int,dim> > lower(3), upper(3);
      lower[0] = vect<int,dim> (l1ixmin, l1iymin, l1izmin);
      upper[0] = vect<int,dim> (l1ixmax, l1iymax, l1izmax);
      lower[1] = vect<int,dim> (l2ixmin, l2iymin, l2izmin);
      upper[1] = vect<int,dim> (l2ixmax, l2iymax, l2izmax);
      lower[2] = vect<int,dim> (l3ixmin, l3iymin, l3izmin);
      upper[2] = vect<int,dim> (l3ixmax, l3iymax, l3izmax);
      MakeRegions_AsSpecified (cctkGH, refinement_levels, lower, upper, bbl);
      
    } else if (CCTK_EQUALS(refined_regions, "manual-coordinates")) {
      
      if (refinement_levels > 4) {
	CCTK_WARN (0, "Cannot currently specify manual refinement regions for more than 4 refinement levels");
      }
      assert (refinement_levels<=4);
      vector<vect<CCTK_REAL,dim> > lower(3), upper(3);
      lower[0] = vect<CCTK_REAL,dim> (l1xmin, l1ymin, l1zmin);
      upper[0] = vect<CCTK_REAL,dim> (l1xmax, l1ymax, l1zmax);
      lower[1] = vect<CCTK_REAL,dim> (l2xmin, l2ymin, l2zmin);
      upper[1] = vect<CCTK_REAL,dim> (l2xmax, l2ymax, l2zmax);
      lower[2] = vect<CCTK_REAL,dim> (l3xmin, l3ymin, l3zmin);
      upper[2] = vect<CCTK_REAL,dim> (l3xmax, l3ymax, l3zmax);
      MakeRegions_AsSpecified (cctkGH, refinement_levels, lower, upper, bbl);
      
    } else if (CCTK_EQUALS(refined_regions, "automatic")) {
      
      const int vi = CCTK_VarIndex (errorvar);
      assert (vi>=0 && vi<CCTK_NumVars());
      const int gi = CCTK_GroupIndexFromVarI (vi);
      assert (gi>=0 && gi<CCTK_NumGroups());
      const int v1 = CCTK_FirstVarIndexI(gi);
      assert (v1<=vi && v1<CCTK_NumVars());
      
      assert (CCTK_GroupTypeI(gi) == CCTK_GF);
      assert (CCTK_VarTypeI(vi) == CCTK_VARIABLE_REAL);
      assert (CCTK_GroupDimI(gi) == dim);
      
      assert (gi < (int)arrdata.size());
      assert (vi-v1 < (int)arrdata[gi].data.size());
      assert (arrdata[gi].data[vi-v1]);
      const gf<CCTK_REAL,dim>& error
	= *dynamic_cast<const gf<CCTK_REAL,dim>*>(arrdata[gi].data[vi-v1]);
      
      MakeRegions_Adaptively (cctkGH, minwidth, minfraction, maxerror, error,
			      bbl);
      
    } else {
      abort();
    }
    
    vector<bbox<int,dim> > bbs;
    bbs.reserve (bbl.size());
    for (list<bbox<int,dim> >::const_iterator it = bbl.begin();
	 it != bbl.end();
	 ++it) {
      bbs.push_back (*it);
    }
    
    // TODO: ensure nesting properties
    
    SplitRegions (cctkGH, bbs);
    
    gh<dim>::cexts bbss;
    bbss = hh->make_reflevel_multigrid_boxes(bbs, mglevels);
    
    bbsss = hh->extents;
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
    
    MakeProcessors (cctkGH, bbsss, pss);
    
    return 1;			// do recompose
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
#if 0
      // refine (arbitrarily) around the lower boundary only
      rlb = rlb;
#endif
      // honour multigrid factors
      const int mgstr = ipow(hh->mgfact, mglevels);
      rlb = (rlb / mgstr) * mgstr;
      rub = rlb + newrextent;
      // require rub<oldrub because we really want rub-rstr<=oldrub-oldstr
      assert (all(rlb >= oldrlb && rub < oldrub));
    }
    
    bbl.push_back (bbox<int,dim>(rlb, rub-rstr, rstr));
  }
  
  
  
  static void
  MakeRegions_AsSpecified_OneLevel (const cGH* cctkGH, const int reflevels,
				    const vect<int,dim> ilower,
				    const vect<int,dim> iupper,
				    list<bbox<int,dim> >& bbl)
  {
    assert (bbl.empty());
    
    if (reflevel+1 >= reflevels) return;
    
    const vect<int,dim> rstr = hh->baseextent.stride();
    const vect<int,dim> rlb  = hh->baseextent.lower();
    const vect<int,dim> rub  = hh->baseextent.upper();
    
    const int rl = reflevel+1;
    
    const int levfac = ipow(hh->reffact, rl);
    assert (all (rstr % levfac == 0));
    const vect<int,dim> str (rstr / levfac);
    const vect<int,dim> lb  (ilower);
    const vect<int,dim> ub  (iupper);
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
  
  
  
  void MakeRegions_AsSpecified (const cGH* cctkGH, const int reflevels,
				const vector<vect<int,dim> > lower,
				const vector<vect<int,dim> > upper,
				list<bbox<int,dim> >& bbl)
  {
    if (reflevel+1 >= reflevels) return;
    
    const int rl = reflevel+1;
    
    const vect<int,dim> ilower = lower[rl-1];
    const vect<int,dim> iupper = upper[rl-1];
    
    MakeRegions_AsSpecified_OneLevel (cctkGH, reflevels, ilower, iupper, bbl);
  }
  
  
  
  void MakeRegions_AsSpecified (const cGH* cctkGH, const int reflevels,
				const vector<vect<CCTK_REAL,dim> > lower,
				const vector<vect<CCTK_REAL,dim> > upper,
				list<bbox<int,dim> >& bbl)
  {
    if (reflevel+1 >= reflevels) return;
    
    vect<CCTK_REAL,dim> global_lower, global_upper;
    for (int d=0; d<dim; ++d) {
      const int ierr = CCTK_CoordRange
	(cctkGH, &global_lower[d], &global_upper[d], d+1, 0, "cart3d");
      if (ierr<0) {
	global_lower[d] = 0;
	global_upper[d] = 1;
      }
    }
    const vect<int,dim> global_extent (hh->baseextent.upper() - hh->baseextent.lower() + hh->baseextent.stride() * (dd->lghosts + dd->ughosts));
    
    const int rl = reflevel+1;
    
    const vect<int,dim> ilower = vect<int,dim>(map((CCTK_REAL(*)(CCTK_REAL))(floor), (lower[rl-1] - global_lower) / (global_upper - global_lower) * vect<CCTK_REAL,dim>(global_extent) + 0.5));
    const vect<int,dim> iupper = vect<int,dim>(map((CCTK_REAL(*)(CCTK_REAL))(floor), (upper[rl-1] - global_lower) / (global_upper - global_lower) * vect<CCTK_REAL,dim>(global_extent) + 0.5));
    
    MakeRegions_AsSpecified_OneLevel (cctkGH, reflevels, ilower, iupper, bbl);
  }
  
  
  
  static void
  MakeRegions_Adaptively_Recombine (list<bbox<int,dim> >& bbl1,
				    list<bbox<int,dim> >& bbl2,
				    list<bbox<int,dim> >& bbl,
				    const bbox<int,dim>& iface,
				    const int dir)
  {
    assert (!iface.empty());
    assert (iface.lower()[dir] == iface.upper()[dir]);
    
    const int oldnumboxes = bbl.size() + bbl1.size() + bbl2.size();
    int numcombinedboxes = 0;
    
    int oldnumpoints = 0;
    for (list<bbox<int,dim> >::const_iterator ibb = bbl.begin(); ibb != bbl.end(); ++ibb) {
      oldnumpoints += ibb->num_points();
    }
    for (list<bbox<int,dim> >::const_iterator ibb1 = bbl1.begin(); ibb1 != bbl1.end(); ++ibb1) {
      oldnumpoints += ibb1->num_points();
    }
    for (list<bbox<int,dim> >::const_iterator ibb2 = bbl2.begin(); ibb2 != bbl2.end(); ++ibb2) {
      oldnumpoints += ibb2->num_points();
    }
    
#if 0
    // remember old bounding boxes
    bboxset<int,dim> oldboxes;
    for (list<bbox<int,dim> >::const_iterator ibb1 = bbl1.begin(); ibb1 != bbl1.end(); ++ibb1) {
      oldboxes += *ibb1;
    }
    for (list<bbox<int,dim> >::const_iterator ibb2 = bbl2.begin(); ibb2 != bbl2.end(); ++ibb2) {
      oldboxes += *ibb2;
    }
#endif
#if 0
    cout << endl;
    cout << "MakeRegions_Adaptively_Recombine: initial list:" << endl;
    cout << bbl << endl;
    cout << "MakeRegions_Adaptively_Recombine: initial list 1:" << endl;
    cout << bbl1 << endl;
    cout << "MakeRegions_Adaptively_Recombine: initial list 2:" << endl;
    cout << bbl2 << endl;
#endif
    
    const vect<int,dim> lo = iface.lower();
    const vect<int,dim> up = iface.upper();
    const vect<int,dim> str = iface.stride();
    
    {
      // prune boxes on the left
      list<bbox<int,dim> >::iterator ibb1 = bbl1.begin();
      while (ibb1 != bbl1.end()) {
	// is this bbox just to the left of the interface?
// 	const vect<int,dim> lo1 = ibb1->lower();
	const vect<int,dim> up1 = ibb1->upper();
	const vect<int,dim> str1 = ibb1->stride();
	assert (up1[dir]+str1[dir] <= lo[dir]);
	assert (all(str1 == str));
	if (up1[dir]+str1[dir] < lo[dir]) {
	  // no: forget it
	  bbl.push_back (*ibb1);
	  ibb1 = bbl1.erase(ibb1);
	  continue;
	}
	++ibb1;
      }	// while
    }
    
    {
      // prune boxes on the right
      list<bbox<int,dim> >::iterator ibb2 = bbl2.begin();
      while (ibb2 != bbl2.end()) {
	// is this bbox just to the right of the interface?
	const vect<int,dim> lo2 = ibb2->lower();
// 	const vect<int,dim> up2 = ibb2->upper();
	const vect<int,dim> str2 = ibb2->stride();
	assert (up[dir] <= lo2[dir]);
	assert (all(str2 == str));
	if (up[dir] < lo2[dir]) {
	  // no: forget it
	  bbl.push_back (*ibb2);
	  ibb2 = bbl2.erase(ibb2);
	  continue;
	}
	++ibb2;
      }	// while
    }
    
    {
      // walk all boxes on the left
      list<bbox<int,dim> >::iterator ibb1 = bbl1.begin();
      while (ibb1 != bbl1.end()) {
	vect<int,dim> lo1 = ibb1->lower();
	vect<int,dim> up1 = ibb1->upper();
	vect<int,dim> str1 = ibb1->stride();
	assert (up1[dir]+str1[dir] == lo[dir]);
	lo1[dir] = lo[dir];
	up1[dir] = up[dir];
	const bbox<int,dim> iface1 (lo1,up1,str1);
	
	{
	  // walk all boxes on the right
	  list<bbox<int,dim> >::iterator ibb2 = bbl2.begin();
	  while (ibb2 != bbl2.end()) {
	    vect<int,dim> lo2 = ibb2->lower();
	    vect<int,dim> up2 = ibb2->upper();
	    vect<int,dim> str2 = ibb2->stride();
	    assert (lo2[dir] == up[dir]);
	    lo2[dir] = lo[dir];
	    up2[dir] = up[dir];
	    const bbox<int,dim> iface2 (lo2,up2,str2);
	    
	    // check for a match
	    if (iface1 == iface2) {
	      const bbox<int,dim> combined (ibb1->lower(), ibb2->upper(), str);
	      bbl.push_back (combined);
	      ibb1 = bbl1.erase(ibb1);
	      ibb2 = bbl2.erase(ibb2);
	      ++numcombinedboxes;
// 	      cout << "MRA: Combining along " << "xyz"[dir] << " to " << bbl.back() << " size " << bbl.back().num_points() << endl;
	      goto continue_search;
	    }
	    
	    ++ibb2;
	  } // while
	}
	
	++ibb1;
      continue_search:;
      }	// while
    }
    
    bbl.splice (bbl.end(), bbl1);
    bbl.splice (bbl.end(), bbl2);
    
    assert (bbl1.empty() && bbl2.empty());
    
    const int newnumboxes = bbl.size();
    assert (newnumboxes + numcombinedboxes == oldnumboxes);
    
    int newnumpoints = 0;
    for (list<bbox<int,dim> >::const_iterator ibb = bbl.begin(); ibb != bbl.end(); ++ibb) {
      newnumpoints += ibb->num_points();
    }
    assert (newnumpoints == oldnumpoints);
    
#if 0
    // find new bounding boxes
    bboxset<int,dim> newboxes;
    for (list<bbox<int,dim> >::const_iterator ibb = bbl.begin(); ibb != bbl.end(); ++ibb) {
      newboxes += *ibb;
    }
    
    // Check that they are equal
    assert (newboxes.size() <= oldboxes.size());
    assert ((newboxes.size()==0) == (oldboxes.size()==0));
    assert (oldboxes == newboxes);
#endif
#if 0
    cout << "MakeRegions_Adaptively_Recombine: final list:" << endl;
    cout << bbl << endl;
    cout << endl;
#endif
  }
  
  static void
  MakeRegions_Adaptively (const cGH* const cctkGH,
			  const int minwidth,
			  const CCTK_REAL minfraction,
			  const CCTK_REAL maxerror,
			  const data<CCTK_REAL,dim>& error,
			  list<bbox<int,dim> >& bbl,
			  const bbox<int,dim>& region)
  {
    // Just to be sure
    assert (! region.empty());
    
    // Count grid points that need to be refined
    // (this doesn't work yet on multiple processors)
    assert (CCTK_nProcs(cctkGH)==1);
    int cnt = 0;
    for (bbox<int,dim>::iterator it=region.begin(); it!=region.end(); ++it) {
      if (error[*it] > maxerror) ++cnt;
    }
    const CCTK_REAL fraction = (CCTK_REAL)cnt / region.num_points();
    const int width = maxval(region.shape() / region.stride());
    
    if (cnt == 0) {
      // Don't refine
    } else if (width < 2*minwidth || fraction >= minfraction) {
      // Refine the whole region
      const vect<int,dim> lo(region.lower());
      const vect<int,dim> up(region.upper());
      const vect<int,dim> str(region.stride());
      bbl.push_back (bbox<int,dim>(lo,up+str-str/reffact,str/reffact));
//       cout << "MRA: Refining to " << bbl.back() << " size " << bbl.back().num_points() << " fraction " << fraction << endl;
    } else {
      // Split the region and check recursively
      const int dir = maxloc(region.shape());
      const vect<int,dim> lo(region.lower());
      const vect<int,dim> up(region.upper());
      const vect<int,dim> str(region.stride());
      vect<int,dim> lo1(lo), lo2(lo);
      vect<int,dim> up1(up), up2(up);
      const int mgstr = ipow(hh->mgfact, mglevels); // honour multigrid factors
      const int step = str[dir]*mgstr;
      lo2[dir] = ((up[dir]+lo[dir])/2/step)*step;
      up1[dir] = lo2[dir]-str[dir];
      const bbox<int,dim> region1(lo1,up1,str);
      const bbox<int,dim> region2(lo2,up2,str);
      assert (region1.is_contained_in(region));
      assert (region2.is_contained_in(region));
      assert ((region1 & region2).empty());
      assert (region1 + region2 == region);
      list<bbox<int,dim> > bbl1, bbl2;
      MakeRegions_Adaptively (cctkGH, minwidth, minfraction, maxerror,
			      error, bbl1, region1);
      MakeRegions_Adaptively (cctkGH, minwidth, minfraction, maxerror,
			      error, bbl2, region2);
      // Combine regions if possible
      up2 += str-str/reffact;
      up2[dir] = lo2[dir];
      const bbox<int,dim> iface(lo2,up2,str/reffact);
      MakeRegions_Adaptively_Recombine (bbl1, bbl2, bbl, iface, dir);
    }
    
  }
  
  void
  MakeRegions_Adaptively (const cGH* const cctkGH,
			  const int minwidth,
			  const CCTK_REAL minfraction,
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
    
//     cout << endl << "MRA: Choosing regions to refine in " << hh->components(rl) << " components" << endl;
    
    for (int c=0; c<hh->components(rl); ++c) {
      const bbox<int,dim> region = hh->extents[rl][c][ml];
      assert (! region.empty());
      
      const data<CCTK_REAL,dim>& errdata = *error(tl,rl,c,ml);
      
      MakeRegions_Adaptively (cctkGH, minwidth, minfraction, maxerror,
			      errdata, bbl, region);
    }
    
    int numpoints = 0;
    for (list<bbox<int,dim> >::const_iterator ibb = bbl.begin(); ibb != bbl.end(); ++ibb) {
      numpoints += ibb->num_points();
    }
//     cout << "MRA: Chose " << bbl.size() << " regions with a total size of " << numpoints << " to refine." << endl << endl;
  }
  
} // namespace CarpetRegrid
