#include <assert.h>
#include <stdlib.h>

#include <vector>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "Carpet/CarpetLib/src/bbox.hh"
#include "Carpet/CarpetLib/src/gh.hh"
#include "Carpet/CarpetLib/src/vect.hh"

#include "carpet.hh"

static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/Recompose.cc,v 1.2 2001/07/09 09:00:09 schnetter Exp $";



namespace Carpet {
  
  using namespace std;
  
  
  
  static void Recompose_gh (cGH* cgh, gh<dim>* hh);
  
  
  
  void Recompose (cGH* cgh)
  {
    DECLARE_CCTK_PARAMETERS;
    
    assert (component == -1);
    Checkpoint ("%*sRecompose", 2*reflevel, "");
    
    Recompose_gh (cgh, hh);
    Recompose_gh (cgh, hh0);
    
    for (int group=0; group<CCTK_NumGroups(); ++group) {
      switch (CCTK_GroupTypeI(group)) {
      case CCTK_SCALAR:
	break;
      case CCTK_ARRAY:
	Recompose_gh (cgh, arrdata[group].hh);
      case CCTK_GF:
	break;
      default:
	abort();
      }	// switch
    } // for
  }
  
  static void Recompose_gh (cGH* cgh, gh<dim>* hh)
  {
    DECLARE_CCTK_PARAMETERS;
    
    const int nprocs    = CCTK_nProcs(cgh);
    const int reflevels = max_refinement_levels; // arbitrary value
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
	const vect<int,dim> rextent = rub - rlb + rstr;
	const vect<int,dim> rcentre = rlb + (rextent / 2 / rstr) * rstr;
	// calculate new extent
	assert (all(rextent % hh->reffact == 0));
	const vect<int,dim> newrextent = rextent / hh->reffact;
	// refined boxes have smaller stride
	assert (all(rstr%hh->reffact == 0));
	rstr /= hh->reffact;
	// refine (arbitrarily) around the center only
	rlb = rcentre - (newrextent/2 / rstr) * rstr;
	rub = rlb + newrextent - rstr;
// 	// refine (arbitrarily) around the lower boundary only
// 	rlb = rlb;
// 	rub = rlb + newrextent - rstr;
	assert (all(rlb >= oldrlb && rub <= oldrub));
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
  
} // namespace Carpet
