#include <assert.h>
#include <math.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>

#include <list>
#include <sstream>
#include <string>
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

static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetRegrid/src/regrid.cc,v 1.13 2002/03/26 13:22:31 schnetter Exp $";

CCTK_FILEVERSION(CarpetRegrid_regrid_cc)



#undef AUTOMATIC_BOUNDARIES



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
			  gh<dim>::rbnds& obss,
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
    
    // Return if we want to regrid during initial data only, and this
    // is not the time for initial data
    if (regrid_every == 0 && cctkGH->cctk_iteration != 0) return 0;
    
    // Return if we want to regrid regularly, but not at this time
    if (regrid_every > 0 && cctkGH->cctk_iteration != 0
	&& cctkGH->cctk_iteration % regrid_every != 0) {
      return 0;
    }
    
    list<ibbox> bbl;
    list<bvect> obl;
    
    if (CCTK_EQUALS(refined_regions, "none")) {
      
      MakeRegions_BaseLevel (cctkGH, bbl, obl);
      
    } else if (CCTK_EQUALS(refined_regions, "centre")) {
      
      MakeRegions_RefineCentre (cctkGH, refinement_levels, bbl, obl);
      
    } else if (CCTK_EQUALS(refined_regions, "manual-gridpoints")) {
      
      if (refinement_levels > 4) {
	CCTK_WARN (0, "Cannot currently specify manual refinement regions for more than 4 refinement levels");
      }
      assert (refinement_levels<=4);
      vector<ivect> lower(3), upper(3);
      lower[0] = ivect (l1ixmin, l1iymin, l1izmin);
      upper[0] = ivect (l1ixmax, l1iymax, l1izmax);
      lower[1] = ivect (l2ixmin, l2iymin, l2izmin);
      upper[1] = ivect (l2ixmax, l2iymax, l2izmax);
      lower[2] = ivect (l3ixmin, l3iymin, l3izmin);
      upper[2] = ivect (l3ixmax, l3iymax, l3izmax);
      MakeRegions_AsSpecified (cctkGH, refinement_levels, lower, upper,
			       bbl, obl);
      
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
      MakeRegions_AsSpecified (cctkGH, refinement_levels, lower, upper,
			       bbl, obl);
      
    } else if (CCTK_EQUALS(refined_regions, "manual-gridpoint-list")) {
      
      vector<vector<ibbox> > bbss;
      if (strcmp(gridpoints, "") !=0 ) {
	istringstream gp_str(gridpoints);
	gp_str >> bbss;
      }
      
      vector<vector<vect<vect<bool,2>,dim> > > obss;
      if (strcmp(outerbounds, "") !=0 ) {
	istringstream ob_str (outerbounds);
	ob_str >> obss;
      } else {
	obss.resize(bbss.size());
	for (int rl=0; rl<(int)obss.size(); ++rl) {
	  obss[rl].resize(bbss[rl].size());
	  for (int c=0; c<(int)obss[rl].size(); ++c) {
	    obss[rl][c] = vect<vect<bool,2>,dim>(vect<bool,2>(false));
	  }
	}
      }
      
      MakeRegions_AsSpecified (cctkGH, refinement_levels, bbss, obss,
			       bbl, obl);
      
    } else if (CCTK_EQUALS(refined_regions, "manual-coordinate-list")) {
      
      abort ();
      
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
			      bbl, obl);
      
    } else {
      abort();
    }
    
#ifdef AUTOMATIC_BOUNDARIES
    assert (obl.size() == 0);
#else
    assert (bbl.size() == obl.size());
#endif
    
    // transform bbox list into bbox vector
    vector<ibbox> bbs;
    bbs.reserve (bbl.size());
    for (list<ibbox>::const_iterator it = bbl.begin();
	 it != bbl.end();
	 ++it) {
      bbs.push_back (*it);
    }
#ifdef AUTOMATIC_BOUNDARIES
    vector<bvect> obs;
    obs.resize (bbs.size());
    for (int c=0; c<(int)bbs.size(); ++c) {
      ivect extlo = hh->baseextent.lower();
      ivect extup = hh->baseextent.upper();
      ivect str = bbs[c].stride();
      for (int d=0; d<dim; ++d) {
	// Decide what is an outer boundary.
	// Allow it to be one grid point to the interior.
	assert (bbs[c].lower()[d] >= extlo[d]);
	obs[c][d][0] = bbs[c].lower()[d] <= extlo[d] + str[d];
	assert (bbs[c].upper()[d] <= extup[d]);
	obs[c][d][1] = bbs[c].upper()[d] >= extup[d] - str[d];
      }
    }
#else
    vector<bvect> obs;
    obs.reserve (obl.size());
    for (list<bvect>::const_iterator it = obl.begin();
	 it != obl.end();
	 ++it) {
      obs.push_back (*it);
    }
#endif
    
    
    
    // TODO: ensure nesting properties
    
    // make multiprocessor aware
    SplitRegions (cctkGH, bbs, obs);
    
    // make multigrid aware
    gh<dim>::cexts bbss;
    bbss = hh->make_reflevel_multigrid_boxes(bbs, mglevels);
    
    // insert into grid hierarchy
    bbsss = hh->extents;
    obss = hh->outer_boundaries;
    if (bbss.size() == 0) {
      // remove all finer levels
      // TODO: this might not work
      bbsss.resize(reflevel+1);
      obss.resize(reflevel+1);
    } else {
      assert (reflevel < (int)bbsss.size());
      if (reflevel+1 == (int)bbsss.size()) {
	// add a finer level
	bbsss.push_back (bbss);
	obss.push_back (obs);
      } else {
	// change a finer level
	bbsss[reflevel+1] = bbss;
	obss[reflevel+1] = obs;
      }
    }
    
    MakeProcessors (cctkGH, bbsss, pss);
    
    return 1;			// do recompose
  }
  
  
  
  void MakeRegions_BaseLevel (const cGH* cctkGH,
			      list<ibbox>& bbl, list<bvect>& obl)
  {
    assert (bbl.empty());
    assert (obl.empty());
  }
  
  
  
  // This is a helpful helper routine. The user can use it to define
  // how the hierarchy should be refined. But the result of this
  // routine is rather arbitrary.
  void MakeRegions_RefineCentre (const cGH* cctkGH, const int reflevels,
				 list<ibbox>& bbl, list<bvect>& obl)
  {
    assert (bbl.empty());
    assert (obl.empty());
    
    if (reflevel+1 >= reflevels) return;
      
    // note: what this routine calls "ub" is "ub+str" elsewhere
    ivect rstr = hh->baseextent.stride();
    ivect rlb  = hh->baseextent.lower();
    ivect rub  = hh->baseextent.upper() + rstr;
    
    for (int rl=0; rl<reflevel+1; ++rl) {
      // save old values
      const ivect oldrlb = rlb;
      const ivect oldrub = rub;
      // calculate extent and centre
      const ivect rextent = rub - rlb;
      const ivect rcentre = rlb + (rextent / 2 / rstr) * rstr;
      // calculate new extent
      assert (all(rextent % hh->reffact == 0));
      const ivect newrextent = rextent / hh->reffact;
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
    
    bbl.push_back (ibbox(rlb, rub-rstr, rstr));
    obl.push_back (bvect(vect<bool,2>(false)));
  }
  
  
  
  static void
  MakeRegions_AsSpecified_OneLevel (const cGH* cctkGH, const int reflevels,
				    const ivect ilower,
				    const ivect iupper,
				    const bvect obound,
				    list<ibbox>& bbl, list<bvect>& obl)
  {
    if (reflevel+1 >= reflevels) return;
    
    const ivect rstr = hh->baseextent.stride();
    const ivect rlb  = hh->baseextent.lower();
    const ivect rub  = hh->baseextent.upper();
    
    const int rl = reflevel+1;
    
    const int levfac = ipow(hh->reffact, rl);
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
    
    bbl.push_back (ibbox(lb, ub, str));
    obl.push_back (obound);
  }
  
  
  
  void MakeRegions_AsSpecified (const cGH* cctkGH, const int reflevels,
				const vector<ivect> lower,
				const vector<ivect> upper,
				list<ibbox>& bbl, list<bvect>& obl)
  {
    assert (lower.size() == upper.size());
    
    if (reflevel+1 >= reflevels) return;
    
    const int rl = reflevel+1;
    
    const ivect ilower = lower[rl-1];
    const ivect iupper = upper[rl-1];
    const bvect obound = bvect(vect<bool,2>(false));
    
    MakeRegions_AsSpecified_OneLevel (cctkGH, reflevels,
				      ilower, iupper, obound,
				      bbl, obl);
  }
  
  
  
  void MakeRegions_AsSpecified (const cGH* cctkGH, const int reflevels,
				const vector<vect<CCTK_REAL,dim> > lower,
				const vector<vect<CCTK_REAL,dim> > upper,
				list<ibbox>& bbl, list<bvect>& obl)
  {
    assert (lower.size() == upper.size());
    
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
    const ivect global_extent (hh->baseextent.upper() - hh->baseextent.lower() + hh->baseextent.stride() * (dd->lghosts + dd->ughosts));
    
    const int rl = reflevel+1;
    
    const vect<CCTK_REAL,dim> scale = vect<CCTK_REAL,dim>(global_extent) / (global_upper - global_lower);
    const vect<CCTK_REAL,dim> rlower = (lower[rl-1] - global_lower) * scale;
    const vect<CCTK_REAL,dim> rupper = (upper[rl-1] - global_lower) * scale;
    const ivect ilower = ivect(map((CCTK_REAL(*)(CCTK_REAL))floor, rlower + 0.5));
    const ivect iupper = ivect(map((CCTK_REAL(*)(CCTK_REAL))floor, rupper + 0.5));
    const bvect obound = bvect(vect<bool,2>(false));
    
    MakeRegions_AsSpecified_OneLevel (cctkGH, reflevels,
				      ilower, iupper, obound,
				      bbl, obl);
  }
  
  
  
  void MakeRegions_AsSpecified (const cGH* cctkGH, const int reflevels,
				const vector<vector<ibbox> > bbss,
				const vector<vector<bvect> > obss,
				list<ibbox>& bbl, list<bvect>& obl)
  {
    if (reflevel+1 >= reflevels) return;
    
    const int rl = reflevel+1;
    
    for (int c=0; c<(int)bbss.size(); ++c) {
      
      const ivect ilower = bbss[rl-1][c].lower();
      const ivect iupper = bbss[rl-1][c].upper();
      const bvect obound = obss[rl-1][c];
      
      MakeRegions_AsSpecified_OneLevel (cctkGH, reflevels,
					ilower, iupper, obound,
					bbl, obl);
      
    }
  }
  
  
  
  static void
  MakeRegions_Adaptively_Recombine (list<ibbox>& bbl1,
				    list<ibbox>& bbl2,
				    list<ibbox>& bbl,
				    const ibbox& iface,
				    const int dir)
  {
    assert (!iface.empty());
    assert (iface.lower()[dir] == iface.upper()[dir]);
    
    const int oldnumboxes = bbl.size() + bbl1.size() + bbl2.size();
    int numcombinedboxes = 0;
    
    int oldnumpoints = 0;
    for (list<ibbox>::const_iterator ibb = bbl.begin(); ibb != bbl.end(); ++ibb) {
      oldnumpoints += ibb->num_points();
    }
    for (list<ibbox>::const_iterator ibb1 = bbl1.begin(); ibb1 != bbl1.end(); ++ibb1) {
      oldnumpoints += ibb1->num_points();
    }
    for (list<ibbox>::const_iterator ibb2 = bbl2.begin(); ibb2 != bbl2.end(); ++ibb2) {
      oldnumpoints += ibb2->num_points();
    }
    
#if 0
    // remember old bounding boxes
    bboxset<int,dim> oldboxes;
    for (list<ibbox>::const_iterator ibb1 = bbl1.begin(); ibb1 != bbl1.end(); ++ibb1) {
      oldboxes += *ibb1;
    }
    for (list<ibbox>::const_iterator ibb2 = bbl2.begin(); ibb2 != bbl2.end(); ++ibb2) {
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
    
    const ivect lo = iface.lower();
    const ivect up = iface.upper();
    const ivect str = iface.stride();
    
    {
      // prune boxes on the left
      list<ibbox>::iterator ibb1 = bbl1.begin();
      while (ibb1 != bbl1.end()) {
	// is this bbox just to the left of the interface?
// 	const ivect lo1 = ibb1->lower();
	const ivect up1 = ibb1->upper();
	const ivect str1 = ibb1->stride();
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
      list<ibbox>::iterator ibb2 = bbl2.begin();
      while (ibb2 != bbl2.end()) {
	// is this bbox just to the right of the interface?
	const ivect lo2 = ibb2->lower();
// 	const ivect up2 = ibb2->upper();
	const ivect str2 = ibb2->stride();
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
      list<ibbox>::iterator ibb1 = bbl1.begin();
      while (ibb1 != bbl1.end()) {
	ivect lo1 = ibb1->lower();
	ivect up1 = ibb1->upper();
	ivect str1 = ibb1->stride();
	assert (up1[dir]+str1[dir] == lo[dir]);
	lo1[dir] = lo[dir];
	up1[dir] = up[dir];
	const ibbox iface1 (lo1,up1,str1);
	
	{
	  // walk all boxes on the right
	  list<ibbox>::iterator ibb2 = bbl2.begin();
	  while (ibb2 != bbl2.end()) {
	    ivect lo2 = ibb2->lower();
	    ivect up2 = ibb2->upper();
	    ivect str2 = ibb2->stride();
	    assert (lo2[dir] == up[dir]);
	    lo2[dir] = lo[dir];
	    up2[dir] = up[dir];
	    const ibbox iface2 (lo2,up2,str2);
	    
	    // check for a match
	    if (iface1 == iface2) {
	      const ibbox combined (ibb1->lower(), ibb2->upper(), str);
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
    for (list<ibbox>::const_iterator ibb = bbl.begin(); ibb != bbl.end(); ++ibb) {
      newnumpoints += ibb->num_points();
    }
    assert (newnumpoints == oldnumpoints);
    
#if 0
    // find new bounding boxes
    bboxset<int,dim> newboxes;
    for (list<ibbox>::const_iterator ibb = bbl.begin(); ibb != bbl.end(); ++ibb) {
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
			  list<ibbox>& bbl,
			  const ibbox& region)
  {
    // Just to be sure
    assert (! region.empty());
    
    // Count grid points that need to be refined
    // (this doesn't work yet on multiple processors)
    assert (CCTK_nProcs(cctkGH)==1);
    int cnt = 0;
    for (ibbox::iterator it=region.begin(); it!=region.end(); ++it) {
      if (error[*it] > maxerror) ++cnt;
    }
    const CCTK_REAL fraction = (CCTK_REAL)cnt / region.num_points();
    const int width = maxval(region.shape() / region.stride());
    
    if (cnt == 0) {
      // Don't refine
    } else if (width < 2*minwidth || fraction >= minfraction) {
      // Refine the whole region
      const ivect lo(region.lower());
      const ivect up(region.upper());
      const ivect str(region.stride());
      bbl.push_back (ibbox(lo,up+str-str/reffact,str/reffact));
//       cout << "MRA: Refining to " << bbl.back() << " size " << bbl.back().num_points() << " fraction " << fraction << endl;
    } else {
      // Split the region and check recursively
      const int dir = maxloc(region.shape());
      const ivect lo(region.lower());
      const ivect up(region.upper());
      const ivect str(region.stride());
      ivect lo1(lo), lo2(lo);
      ivect up1(up), up2(up);
      const int mgstr = ipow(hh->mgfact, mglevels); // honour multigrid factors
      const int step = str[dir]*mgstr;
      lo2[dir] = ((up[dir]+lo[dir])/2/step)*step;
      up1[dir] = lo2[dir]-str[dir];
      const ibbox region1(lo1,up1,str);
      const ibbox region2(lo2,up2,str);
      assert (region1.is_contained_in(region));
      assert (region2.is_contained_in(region));
      assert ((region1 & region2).empty());
      assert (region1 + region2 == region);
      list<ibbox> bbl1, bbl2;
      MakeRegions_Adaptively (cctkGH, minwidth, minfraction, maxerror,
			      error, bbl1, region1);
      MakeRegions_Adaptively (cctkGH, minwidth, minfraction, maxerror,
			      error, bbl2, region2);
      // Combine regions if possible
      up2 += str-str/reffact;
      up2[dir] = lo2[dir];
      const ibbox iface(lo2,up2,str/reffact);
      MakeRegions_Adaptively_Recombine (bbl1, bbl2, bbl, iface, dir);
    }
    
  }
  
  void
  MakeRegions_Adaptively (const cGH* const cctkGH,
			  const int minwidth,
			  const CCTK_REAL minfraction,
			  const CCTK_REAL maxerror,
			  const gf<CCTK_REAL,dim>& error,
			  list<ibbox>& bbl, list<bvect>& obl)
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
      const ibbox region = hh->extents[rl][c][ml];
      assert (! region.empty());
      
      const data<CCTK_REAL,dim>& errdata = *error(tl,rl,c,ml);
      
      MakeRegions_Adaptively (cctkGH, minwidth, minfraction, maxerror,
			      errdata, bbl, region);
    }
    
//     int numpoints = 0;
//     for (list<ibbox>::const_iterator ibb = bbl.begin(); ibb != bbl.end(); ++ibb) {
//       numpoints += ibb->num_points();
//     }
//     cout << "MRA: Chose " << bbl.size() << " regions with a total size of " << numpoints << " to refine." << endl << endl;
    
    // TODO: remove grid points outside the outer boundary
    
    // TODO: create obl depending on the boundary position
    // (at or close to the boundary)
    
  }
  
} // namespace CarpetRegrid
