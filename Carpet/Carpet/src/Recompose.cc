#include <assert.h>
#include <stdlib.h>

#include <vector>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "Carpet/CarpetLib/src/bbox.hh"
#include "Carpet/CarpetLib/src/gh.hh"
#include "Carpet/CarpetLib/src/vect.hh"

#include "carpet.hh"

static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/Recompose.cc,v 1.6 2001/10/29 08:36:45 schnetter Exp $";



namespace Carpet {
  
  using namespace std;
  
  
  
  static bool do_recompose = false;
  static gh<dim>::rexts next_bbsss;
  static gh<dim>::rprocs next_pss;
  
  
  
  static void Adapt (cGH* cgh, int reflevels, gh<dim>* hh);
  
  
  
  void RegisterRecomposeRegions (const gh<dim>::rexts& bbsss,
				 const gh<dim>::rprocs& pss)
  {
    // save the region information for the next regridding
    next_bbsss = bbsss;
    next_pss = pss;
    do_recompose = true;
  }
  
  
  
  void Recompose (cGH* cgh)
  {
    DECLARE_CCTK_PARAMETERS;
    
    assert (component == -1);
    Checkpoint ("%*sRecompose", 2*reflevel, "");
    
    // Check whether to recompose
    if (!do_recompose) return;
    
    // Recompose
    hh->recompose (next_bbsss, next_pss);
    
    if (verbose && CCTK_MyProc(cgh)==0) {
      cout << endl;
      cout << "New bounding boxes:" << endl;
      for (int rl=0; rl<hh->reflevels(); ++rl) {
	for (int c=0; c<hh->components(rl); ++c) {
	  for (int ml=0; ml<hh->mglevels(rl,c); ++ml) {
	    cout << "   rl " << rl << "   c " << c << "   ml " << ml
		 << "   bbox " << hh->extents[rl][c][ml] << endl;
	  }
	}
      }
      cout << endl;
      cout << "New processor distribution:" << endl;
      for (int rl=0; rl<hh->reflevels(); ++rl) {
	for (int c=0; c<hh->components(rl); ++c) {
	  cout << "   rl " << rl << "   c " << c
	       << "   processor " << hh->processors[rl][c] << endl;
	}
      }
      cout << endl;
    }
    
    // Don't recompose to these regions any more
    do_recompose = false;
    
    // Adapt grid scalars
    Adapt (cgh, hh->reflevels(), hh0);
    
    // Adapt grid arrays
    for (int group=0; group<CCTK_NumGroups(); ++group) {
      switch (CCTK_GroupTypeI(group)) {
      case CCTK_SCALAR:
	break;
      case CCTK_ARRAY:
	Adapt (cgh, hh->reflevels(), arrdata[group].hh);
      case CCTK_GF:
	break;
      default:
	abort();
      }	// switch
    } // for
  }
  
  
  
  static void Adapt (cGH* cgh, const int reflevels, gh<dim>* hh)
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
	if (cub[dim-1] > rub[dim-1]) cub[dim-1] = rub[dim-1];
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
  
  
  
  void MakeRegions_RefineCentre (cGH* cgh, const int reflevels,
				 gh<dim>::rexts& bbsss, gh<dim>::rprocs& pss)
  {
    DECLARE_CCTK_PARAMETERS;
    
    const int nprocs    = CCTK_nProcs(cgh);
    const int mglevels  = 1;	// arbitrary value
    
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
	if (c == nprocs-1) cub[dim-1] = rub[dim-1];
	assert (cub[dim-1] <= rub[dim-1]);
	bbs[c] = bbox<int,dim>(clb, cub-cstr, cstr);
      }
      bbss[rl] = bbs;
    }
    bbsss = hh->make_multigrid_boxes(bbss, mglevels);
    
    pss.resize(bbss.size());
    for (int rl=0; rl<reflevels; ++rl) {
      pss[rl] = vector<int>(bbss[rl].size());
      // make sure all processors have the same number of components
      assert (bbss[rl].size() % nprocs == 0);
      for (int c=0; c<(int)bbss[rl].size(); ++c) {
	pss[rl][c] = c % nprocs; // distribute among processors
      }
    }
  }
  
} // namespace Carpet
