#include <assert.h>
#include <math.h>
#include <stdlib.h>

#include <vector>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "Carpet/CarpetLib/src/bbox.hh"
#include "Carpet/CarpetLib/src/gh.hh"
#include "Carpet/CarpetLib/src/vect.hh"

#include "carpet.hh"

static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/Recompose.cc,v 1.10 2001/12/05 19:04:45 schnetter Exp $";



namespace Carpet {
  
  using namespace std;
  
  
  
  static bool do_recompose = false;
  static gh<dim>::rexts next_bbsss;
  static gh<dim>::rprocs next_pss;
  
  
  
  static void Adapt (const cGH* cgh, int reflevels, gh<dim>* hh);
  static void Output (const cGH* cgh, const gh<dim>* hh, const char* descr);
  
  
  
  void RegisterRecomposeRegions (const gh<dim>::rexts& bbsss,
				 const gh<dim>::rprocs& pss)
  {
    // save the region information for the next regridding
    next_bbsss = bbsss;
    next_pss = pss;
    do_recompose = true;
  }
  
  
  
  void Recompose (const cGH* cgh)
  {
    DECLARE_CCTK_PARAMETERS;
    
    assert (component == -1);
    Checkpoint ("%*sRecompose", 2*reflevel, "");
    
    // Check whether to recompose
    if (!do_recompose) return;
    
    // Recompose
    hh->recompose (next_bbsss, next_pss);
    Output (cgh, hh, 0);
    
    // Don't recompose to these regions any more
    do_recompose = false;
    
    // Adapt grid scalars
    Adapt (cgh, hh->reflevels(), hh0);
    Output (cgh, hh0, "");
    
    // Adapt grid arrays
    for (int group=0; group<CCTK_NumGroups(); ++group) {
      switch (CCTK_GroupTypeI(group)) {
      case CCTK_SCALAR:
	break;
      case CCTK_ARRAY:
	Adapt (cgh, hh->reflevels(), arrdata[group].hh);
	Output (cgh, arrdata[group].hh, CCTK_GroupName(group));
      case CCTK_GF:
	break;
      default:
	abort();
      }	// switch
    } // for
  }
  
  
  
  // This routine is a leftover. It determines "automatically" how
  // scalars and arrays should be refined. The user really should have
  // a possibility to define how arrays are to be refined.
  static void Adapt (const cGH* cgh, const int reflevels, gh<dim>* hh)
  {
    DECLARE_CCTK_PARAMETERS;
    
    const int nprocs    = CCTK_nProcs(cgh);
    const int mglevels  = 1;	// for now
    vector<vector<bbox<int,dim> > > bbss(reflevels);
    // note: what this routine calls "ub" is "ub+str" elsewhere
    vect<int,dim> rstr = hh->baseextent.stride();
    vect<int,dim> rlb  = hh->baseextent.lower();
    vect<int,dim> rub  = hh->baseextent.upper() + rstr;
    for (int rl=0; rl<reflevels; ++rl) {
      if (rl>0) {
	// save old values
	const vect<int,dim> oldrlb = rlb;
	const vect<int,dim> oldrub = rub;
	// calculate extent
	const vect<int,dim> rextent = rub - rlb;
	// calculate new extent
	assert (all(rextent % hh->reffact == 0));
	const vect<int,dim> newrextent = rextent / hh->reffact;
	// refined boxes have smaller stride
	assert (all(rstr%hh->reffact == 0));
	rstr /= hh->reffact;
 	// refine around the lower boundary only
 	rlb = rlb;
	rub = rlb + newrextent;
	// require rub<oldrub because we really want rub-rstr<=oldrub-oldstr
	assert (all(rlb >= oldrlb && rub < oldrub));
      }
      vector<bbox<int,dim> > bbs(nprocs);
      for (int c=0; c<nprocs; ++c) {
	vect<int,dim> cstr = rstr;
	vect<int,dim> clb = rlb;
	vect<int,dim> cub = rub;
	const int glonpz = (rub[dim-1] - rlb[dim-1]) / cstr[dim-1];
	const int locnpz = (glonpz + nprocs - 1) / nprocs;
	const int zstep = locnpz * cstr[dim-1];
	clb[dim-1] = rlb[dim-1] + zstep *  c;
	cub[dim-1] = rlb[dim-1] + zstep * (c+1);
 	if (clb[dim-1] > rub[dim-1]) clb[dim-1] = rub[dim-1];
 	if (cub[dim-1] > rub[dim-1]) cub[dim-1] = rub[dim-1];
	assert (clb[dim-1] <= cub[dim-1]);
	assert (cub[dim-1] <= rub[dim-1]);
	bbs[c] = bbox<int,dim>(clb, cub-cstr, cstr);
      }
      bbss[rl] = bbs;
    }
    vector<vector<vector<bbox<int,dim> > > > bbsss
      = hh->make_multigrid_boxes(bbss, mglevels);
    
    vector<vector<int> > pss(bbss.size());
    for (int rl=0; rl<reflevels; ++rl) {
      pss[rl] = vector<int>(bbss[rl].size());
      // make sure all processors have the same number of components
      assert (bbss[rl].size() % nprocs == 0);
      for (int c=0; c<(int)bbss[rl].size(); ++c) {
	pss[rl][c] = c % nprocs; // distribute among processors
      }
    }
    
    hh->recompose(bbsss, pss);
  }
  
  
  
  // This is a helpful helper routine. The user can use it to define
  // how the hierarchy should be refined. But the result of this
  // routine is rather arbitrary.
  void MakeRegions_RefineCentre (const cGH* cgh, const int reflevels,
				 gh<dim>::rexts& bbsss)
  {
    DECLARE_CCTK_PARAMETERS;
    
    const int mglevels = 1;	// arbitrary value
    
    vector<vector<bbox<int,dim> > > bbss(reflevels);
    // note: what this routine calls "ub" is "ub+str" elsewhere
    vect<int,dim> rstr = hh->baseextent.stride();
    vect<int,dim> rlb  = hh->baseextent.lower();
    vect<int,dim> rub  = hh->baseextent.upper() + rstr;
    
    for (int rl=0; rl<reflevels; ++rl) {
      if (rl>0) {
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
// 	// refine (arbitrarily) around the lower boundary only
// 	rlb = rlb;
	rub = rlb + newrextent;
	// require rub<oldrub because we really want rub-rstr<=oldrub-oldstr
	assert (all(rlb >= oldrlb && rub < oldrub));
      }
      vector<bbox<int,dim> > bbs(1);
      bbs[0] = bbox<int,dim>(rlb, rub-rstr, rstr);
      bbss[rl] = bbs;
    }
    bbsss = hh->make_multigrid_boxes(bbss, mglevels);
  }
  
  
  
  void MakeRegions_AsSpecified (const cGH* cgh, const int reflevels,
				gh<dim>::rexts& bbsss)
  {
    DECLARE_CCTK_PARAMETERS;
    
    const int mglevels = 1;	// arbitrary value
    
    const vect<int,dim> rstr = hh->baseextent.stride();
    const vect<int,dim> rlb  = hh->baseextent.lower();
    const vect<int,dim> rub  = hh->baseextent.upper();
    
    if (reflevels>4) {
      CCTK_WARN (0, "Cannot currently specify refinement regions for more than 4 refinement levels");
    }
    
    assert (reflevels<4);
    vector<vect<int,dim> > lower(4), upper(4);
    lower[0] = rlb;
    upper[0] = rub;
    lower[1] = vect<int,dim> (l1xmin, l1ymin, l1zmin);
    upper[1] = vect<int,dim> (l1xmax, l1ymax, l1zmax);
    lower[2] = vect<int,dim> (l2xmin, l2ymin, l2zmin);
    upper[2] = vect<int,dim> (l2xmax, l2ymax, l2zmax);
    lower[3] = vect<int,dim> (l3xmin, l3ymin, l3zmin);
    upper[3] = vect<int,dim> (l3xmax, l3ymax, l3zmax);
    
    vector<vector<bbox<int,dim> > > bbss(reflevels);
    
    for (int rl=0; rl<reflevels; ++rl) {
      const int levfac = floor(pow(hh->reffact, rl) + 0.5);
      assert (all (rstr % levfac == 0));
      const vect<int,dim> str (rstr / levfac);
      const vect<int,dim> lb  (lower[rl]);
      const vect<int,dim> ub  (upper[rl]);
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
      vector<bbox<int,dim> > bbs(1);
      bbs[0] = bbox<int,dim>(lb, ub, str);
      bbss[rl] = bbs;
    }
    
    if (reflevels>0) {
      assert (bbss[0].size() > 0);
      assert (bbss[0][0] == hh->baseextent);
    }
    
    bbsss = hh->make_multigrid_boxes(bbss, mglevels);
  }
  
  
  
  // This is a helpful helper routine. The user can use it to define
  // how the hierarchy should be refined. But the result of this
  // routine is rather arbitrary.
  void SplitRegions_AlongZ (const cGH* cgh, gh<dim>::rexts& bbsss)
  {
    DECLARE_CCTK_PARAMETERS;
    
    const int nprocs    = CCTK_nProcs(cgh);
    const int mglevels  = 1;	// arbitrary value
    
    vector<vector<bbox<int,dim> > > bbss(bbsss.size());
    
    for (int rl=0; rl<(int)bbsss.size(); ++rl) {
      assert (bbsss[rl].size() == 1);
      assert (bbsss[rl][0].size() == 1);
      
      const vect<int,dim> rstr = bbsss[rl][0][0].stride();
      const vect<int,dim> rlb  = bbsss[rl][0][0].lower();
      const vect<int,dim> rub  = bbsss[rl][0][0].upper() + rstr;
      
      vector<bbox<int,dim> > bbs(nprocs);
      for (int c=0; c<nprocs; ++c) {
	vect<int,dim> cstr = rstr;
	vect<int,dim> clb = rlb;
	vect<int,dim> cub = rub;
	const int glonpz = (rub[dim-1] - rlb[dim-1]) / cstr[dim-1];
	const int locnpz = (glonpz + nprocs - 1) / nprocs;
	const int zstep = locnpz * cstr[dim-1];
	clb[dim-1] = rlb[dim-1] + zstep *  c;
	cub[dim-1] = rlb[dim-1] + zstep * (c+1);
 	if (clb[dim-1] > rub[dim-1]) clb[dim-1] = rub[dim-1];
 	if (cub[dim-1] > rub[dim-1]) cub[dim-1] = rub[dim-1];
	assert (clb[dim-1] <= cub[dim-1]);
	assert (cub[dim-1] <= rub[dim-1]);
	bbs[c] = bbox<int,dim>(clb, cub-cstr, cstr);
      }
      bbss[rl] = bbs;
    }
    bbsss = hh->make_multigrid_boxes(bbss, mglevels);
  }
  
  
  
  void SplitRegions_AsSpecified (const cGH* cgh, gh<dim>::rexts& bbsss)
  {
    DECLARE_CCTK_PARAMETERS;
    
    const int nprocs    = CCTK_nProcs(cgh);
    const int mglevels  = 1;	// arbitrary value
    
    vector<vector<bbox<int,dim> > > bbss(bbsss.size());
    
    for (int rl=0; rl<(int)bbsss.size(); ++rl) {
      assert (bbsss[rl].size() == 1);
      assert (bbsss[rl][0].size() == 1);
      
      const vect<int,dim> rstr = bbsss[rl][0][0].stride();
      const vect<int,dim> rlb  = bbsss[rl][0][0].lower();
      const vect<int,dim> rub  = bbsss[rl][0][0].upper() + rstr;
      
      const vect<int,dim> nprocs_dir
	(processor_topology_3d_x, processor_topology_3d_y,
	 processor_topology_3d_z);
      assert (all (nprocs_dir > 0));
      if (prod(nprocs_dir) != nprocs) {
	CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
		    "The specified processor topology [%d,%d,%d] does not fit the number of processors, which is %d", nprocs_dir[0], nprocs_dir[1], nprocs_dir[2], nprocs);
      }
      assert (prod(nprocs_dir) == nprocs);
      
      vector<bbox<int,dim> > bbs(nprocs);
      assert (dim==3);
      for (int k=0; k<nprocs_dir[2]; ++k) {
	for (int j=0; j<nprocs_dir[1]; ++j) {
	  for (int i=0; i<nprocs_dir[0]; ++i) {
	    const int c = i + nprocs_dir[0] * (j + nprocs_dir[1] * k);
	    const vect<int,dim> ipos (i, j, k);
	    vect<int,dim> cstr = rstr;
	    vect<int,dim> clb = rlb;
	    vect<int,dim> cub = rub;
	    const vect<int,dim> glonp = (rub - rlb) / cstr;
	    const vect<int,dim> locnp = (glonp + nprocs_dir - 1) / nprocs_dir;
	    const vect<int,dim> step = locnp * cstr;
	    clb = rlb + step *  ipos;
	    cub = rlb + step * (ipos+1);
	    clb = min (clb, rub);
	    cub = min (cub, rub);
	    assert (all (clb <= cub));
	    assert (all (cub <= rub));
	    bbs[c] = bbox<int,dim>(clb, cub-cstr, cstr);
	  }
	}
      }
      bbss[rl] = bbs;
    }
    bbsss = hh->make_multigrid_boxes(bbss, mglevels);
  }
  
  
  
  // This is a helpful helper routine. The user can use it to define
  // how the hierarchy should be refined. But the result of this
  // routine is rather arbitrary.
  void MakeProcessors_RoundRobin (const cGH* cgh, const gh<dim>::rexts& bbsss,
				  gh<dim>::rprocs& pss)
  {
    DECLARE_CCTK_PARAMETERS;
    
    const int nprocs    = CCTK_nProcs(cgh);
    
    pss.resize(bbsss.size());
    for (int rl=0; rl<(int)bbsss.size(); ++rl) {
      pss[rl] = vector<int>(bbsss[rl].size());
      // make sure all processors have the same number of components
      assert (bbsss[rl].size() % nprocs == 0);
      for (int c=0; c<(int)bbsss[rl].size(); ++c) {
	pss[rl][c] = c % nprocs; // distribute among processors
      }
    }
  }
  
  
  
  static void Output (const cGH* cgh, const gh<dim>* hh, const char* descr)
  {
    DECLARE_CCTK_PARAMETERS;
    
    if (verbose && CCTK_MyProc(cgh)==0) {
      cout << endl;
      cout << "New bounding boxes";
      if (descr) {
	if (strlen(descr)) {
	  cout << " for group " << descr;
	} else {
	  cout << " for scalars";
	}
      }
      cout  << ":" << endl;
      for (int rl=0; rl<hh->reflevels(); ++rl) {
	for (int c=0; c<hh->components(rl); ++c) {
	  for (int ml=0; ml<hh->mglevels(rl,c); ++ml) {
	    cout << "   rl " << rl << "   c " << c << "   ml " << ml
		 << "   bbox " << hh->extents[rl][c][ml] << endl;
	  }
	}
      }
      cout << endl;
      cout << "New processor distribution";
      if (descr) {
	if (strlen(descr)) {
	  cout << " for group " << descr;
	} else {
	  cout << " for scalars";
	}
      }
      cout  << ":" << endl;
      for (int rl=0; rl<hh->reflevels(); ++rl) {
	for (int c=0; c<hh->components(rl); ++c) {
	  cout << "   rl " << rl << "   c " << c
	       << "   processor " << hh->processors[rl][c] << endl;
	}
      }
      cout << endl;
    }
  }
  
} // namespace Carpet
