#include <assert.h>
#include <stdlib.h>

#include <list>
#include <vector>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "Carpet/CarpetLib/src/bbox.hh"
#include "Carpet/CarpetLib/src/bboxset.hh"
#include "Carpet/CarpetLib/src/defs.hh"
#include "Carpet/CarpetLib/src/gh.hh"
#include "Carpet/CarpetLib/src/vect.hh"

#include "carpet.hh"

static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/Recompose.cc,v 1.18 2002/01/11 17:19:45 schnetter Exp $";



namespace Carpet {
  
  using namespace std;
  
  
  
  static int (*regrid_routine) (const cGH * cckgGH,
				gh<dim>::rexts& bbsss,
				gh<dim>::rprocs& pss) = 0;
  
  
  
  static void SplitRegions_AlongZ       (const cGH* cgh,
					 vector<bbox<int,dim> >& bbs);
  static void SplitRegions_AsSpecified  (const cGH* cgh,
					 vector<bbox<int,dim> >& bbs);
  
  static void MakeProcessors_RoundRobin (const cGH* cgh,
					 const gh<dim>::rexts& bbss,
					 gh<dim>::rprocs& pss);
  
  
  
  static void CheckRegions (const gh<dim>::rexts& bbsss,
			    const gh<dim>::rprocs& pss);
  static void Adapt (const cGH* cgh, int reflevels, gh<dim>* hh);
  static void Output (const cGH* cgh, const gh<dim>* hh, const char* descr);
  
  
  
  void CheckRegions (const gh<dim>::rexts& bbsss, const gh<dim>::rprocs& pss)
  {
    // At least one level
    assert (bbsss.size() > 0);
    for (int rl=0; rl<(int)bbsss.size(); ++rl) {
      // No empty levels
      assert (bbsss[rl].size() > 0);
      for (int c=0; c<(int)bbsss[rl].size(); ++c) {
	// At least one multigrid level
	assert (bbsss[rl][c].size() > 0);
	for (int ml=0; ml<(int)bbsss[rl][c].size(); ++ml) {
	  // Check sizes
	  assert (all(bbsss[rl][c][ml].lower() <= bbsss[rl][c][ml].upper()));
	  // Check strides
	  const int str = ipow(reffact, maxreflevels-rl-1) * ipow(mgfact, ml);
	  assert (all(bbsss[rl][c][ml].stride() == str));
	  // Check alignments
	  assert (all(bbsss[rl][c][ml].lower() % str == 0));
	  assert (all(bbsss[rl][c][ml].upper() % str == 0));
	}
      }
    }
    
    assert (pss.size() == bbsss.size());
    for (int rl=0; rl<(int)bbsss.size(); ++rl) {
      assert (pss[rl].size() == bbsss[rl].size());
    }
  }
  
  
  
  void RegisterRegridRoutine (int (*routine)(const cGH * cckgGH,
					     gh<dim>::rexts& bbsss,
					     gh<dim>::rprocs& pss))
  {
    assert (!regrid_routine);
    regrid_routine = routine;
  }
  
  
  
  void Regrid (const cGH* cgh)
  {
    assert (mglevel == -1);
    assert (component == -1);
    
    // Check whether to recompose
    gh<dim>::rexts bbsss;
    gh<dim>::rprocs pss;
    int do_recompose = (*regrid_routine) (cgh, bbsss, pss);
    assert (do_recompose >= 0);
    if (do_recompose == 0) return;
    Recompose (cgh, bbsss, pss);
  }
  
  
  
  void Recompose (const cGH* const cgh,
		  const gh<dim>::rexts& bbsss,
		  const gh<dim>::rprocs& pss)
  {
    assert (mglevel == -1);
    assert (component == -1);
    
    // Check the regions
    CheckRegions (bbsss, pss);
    
    // Recompose
    hh->recompose (bbsss, pss);
    Output (cgh, hh, 0);
    
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
  
  
  
  // This routine is a leftover.  It determines "automatically" how
  // scalars and arrays should be refined.  The user really should
  // have a possibility to define how arrays are to be refined.
  static void Adapt (const cGH* cgh, const int reflevels, gh<dim>* hh)
  {
    const int nprocs   = CCTK_nProcs(cgh);
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
//       bbss[rl] = bbs;
      bbss[rl].clear();
      assert (Carpet::hh->components(rl) % nprocs == 0);
      for (int cc=0; cc<(int)Carpet::hh->components(rl)/nprocs; ++cc) {
	bbss[rl].insert(bbss[rl].end(), bbs.begin(), bbs.end());
      }
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
  
  
  
  static void Output (const cGH* cgh, const gh<dim>* hh, const char* descr)
  {
    DECLARE_CCTK_PARAMETERS;
    
    if (verbose) {
      cout << endl;
      cout << "New bounding boxes";
      if (!descr) {
	cout << " for grid functions";
      } else {
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
      if (!descr) {
	cout << " for grid functions";
      } else {
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
  
  
  
  void SplitRegions (const cGH* cgh, vector<bbox<int,dim> >& bbs)
  {
    DECLARE_CCTK_PARAMETERS;
    
    if (CCTK_EQUALS (processor_topology, "automatic")) {
      SplitRegions_AlongZ (cgh, bbs);
    } else if (CCTK_EQUALS (processor_topology, "manual")) {
      SplitRegions_AsSpecified (cgh, bbs);
    } else {
      abort();
    }
  }
  
  
  
  void SplitRegions_AlongZ (const cGH* cgh, vector<bbox<int,dim> >& bbs)
  {
    // Something to do?
    if (bbs.size() == 0) return;
    
    const int nprocs = CCTK_nProcs(cgh);
    
    if (nprocs==1) return;
    
    assert (bbs.size() == 1);
    
    const vect<int,dim> rstr = bbs[0].stride();
    const vect<int,dim> rlb  = bbs[0].lower();
    const vect<int,dim> rub  = bbs[0].upper() + rstr;
    
    bbs.resize(nprocs);
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
  }
  
  
  
  void SplitRegions_AsSpecified (const cGH* cgh, vector<bbox<int,dim> >& bbs)
  {
    DECLARE_CCTK_PARAMETERS;
    
    // Something to do?
    if (bbs.size() == 0) return;
    
    const int nprocs = CCTK_nProcs(cgh);
    
    if (nprocs==1) return;
    
    assert (bbs.size() == 1);
    
    const vect<int,dim> rstr = bbs[0].stride();
    const vect<int,dim> rlb  = bbs[0].lower();
    const vect<int,dim> rub  = bbs[0].upper() + rstr;
    
    const vect<int,dim> nprocs_dir
      (processor_topology_3d_x, processor_topology_3d_y,
       processor_topology_3d_z);
    assert (all (nprocs_dir > 0));
    if (prod(nprocs_dir) != nprocs) {
      CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
		  "The specified processor topology [%d,%d,%d] does not fit the number of processors, which is %d", nprocs_dir[0], nprocs_dir[1], nprocs_dir[2], nprocs);
    }
    assert (prod(nprocs_dir) == nprocs);
    
    bbs.resize(nprocs);
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
  }
  
  
  
  void MakeProcessors (const cGH* cgh, const gh<dim>::rexts& bbsss,
		       gh<dim>::rprocs& pss)
  {
    MakeProcessors_RoundRobin (cgh, bbsss, pss);
  }
  
  
  
  // This is a helpful helper routine. The user can use it to define
  // how the hierarchy should be refined.  But the result of this
  // routine is rather arbitrary.
  void MakeProcessors_RoundRobin (const cGH* cgh, const gh<dim>::rexts& bbsss,
				  gh<dim>::rprocs& pss)
  {
    const int nprocs = CCTK_nProcs(cgh);
    
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
  
} // namespace Carpet
