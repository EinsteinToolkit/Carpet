#include <assert.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <fstream>
#include <iomanip>
#include <list>
#include <vector>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "CactusBase/IOUtil/src/ioutil_CheckpointRecovery.h"

#include "Carpet/CarpetLib/src/bbox.hh"
#include "Carpet/CarpetLib/src/bboxset.hh"
#include "Carpet/CarpetLib/src/defs.hh"
#include "Carpet/CarpetLib/src/gh.hh"
#include "Carpet/CarpetLib/src/vect.hh"

#include "carpet.hh"

static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/Recompose.cc,v 1.22 2002/03/23 20:20:54 schnetter Exp $";



namespace Carpet {
  
  using namespace std;
  
  
  
  typedef vect<int,dim> ivect;
  typedef bbox<int,dim> ibbox;
  
  typedef vect<vect<bool,2>,dim> bvect;
  
  
  
  static int (*regrid_routine) (const cGH * cckgGH,
				gh<dim>::rexts& bbsss,
				gh<dim>::rbnds& obss,
				gh<dim>::rprocs& pss) = 0;
  
  
  
  static void CheckRegions (const gh<dim>::rexts& bbsss,
			    const gh<dim>::rbnds& obss,
			    const gh<dim>::rprocs& pss);
  static void Adapt (const cGH* cgh, int reflevels, gh<dim>* hh);
  
  static void Output (const cGH* cgh, const gh<dim>* hh, const char* descr);
  
  static void OutputGridStructure (const cGH *cgh,
				   const gh<dim>::rexts& bbsss,
				   const gh<dim>::rbnds& obss,
				   const gh<dim>::rprocs& pss);
  
  
  
  static void SplitRegions_AlongZ       (const cGH* cgh,
					 vector<ibbox>& bbs,
					 vector<bvect>& obs);
  static void SplitRegions_AsSpecified  (const cGH* cgh,
					 vector<ibbox>& bbs,
					 vector<bvect>& obs);
  
  static void MakeProcessors_RoundRobin (const cGH* cgh,
					 const gh<dim>::rexts& bbss,
					 gh<dim>::rprocs& pss);
  
  
  
  void CheckRegions (const gh<dim>::rexts& bbsss, const gh<dim>::rbnds& obss,
		     const gh<dim>::rprocs& pss)
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
    assert (obss.size() == bbsss.size());
    for (int rl=0; rl<(int)bbsss.size(); ++rl) {
      assert (obss[rl].size() == bbsss[rl].size());
      assert (pss[rl].size() == bbsss[rl].size());
    }
  }
  
  
  
  void RegisterRegridRoutine (int (*routine)(const cGH * cckgGH,
					     gh<dim>::rexts& bbsss,
					     gh<dim>::rbnds& obss,
					     gh<dim>::rprocs& pss))
  {
    assert (!regrid_routine);
    regrid_routine = routine;
  }
  
  
  
  void Regrid (const cGH* cgh)
  {
    assert (mglevel == -1);
    assert (component == -1);
    
    if (!regrid_routine) {
      static bool didtell = false;
      if (!didtell) {
	CCTK_WARN (1, "No regridding routine has been registered.  There will be no regridding.  (Maybe you forgot to activate the regridding thorn?)");
	didtell = true;
      }
      return;
    }
    
    // Check whether to recompose
    gh<dim>::rexts bbsss;
    gh<dim>::rbnds obss;
    gh<dim>::rprocs pss;
    int do_recompose = (*regrid_routine) (cgh, bbsss, obss, pss);
    assert (do_recompose >= 0);
    if (do_recompose == 0) return;
    Recompose (cgh, bbsss, obss, pss);
  }
  
  
  
  void Recompose (const cGH* const cgh,
		  const gh<dim>::rexts& bbsss,
		  const gh<dim>::rbnds& obss,
		  const gh<dim>::rprocs& pss)
  {
    assert (mglevel == -1);
    assert (component == -1);
    
    // Check the regions
    CheckRegions (bbsss, obss, pss);
    
    // Write grid structure to file
    OutputGridStructure (cgh, bbsss, obss, pss);
    
    // Recompose
    hh->recompose (bbsss, obss, pss);
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
	break;
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
    vector<vector<ibbox> > bbss(reflevels);
    // note: what this routine calls "ub" is "ub+str" elsewhere
    ivect rstr = hh->baseextent.stride();
    ivect rlb  = hh->baseextent.lower();
    ivect rub  = hh->baseextent.upper() + rstr;
    for (int rl=0; rl<reflevels; ++rl) {
      if (rl>0) {
	// save old values
	const ivect oldrlb = rlb;
	const ivect oldrub = rub;
	// calculate extent
	const ivect rextent = rub - rlb;
	// calculate new extent
	assert (all(rextent % hh->reffact == 0));
	const ivect newrextent = rextent / hh->reffact;
	// refined boxes have smaller stride
	assert (all(rstr%hh->reffact == 0));
	rstr /= hh->reffact;
 	// refine around the lower boundary only
 	rlb = rlb;
	rub = rlb + newrextent;
	// require rub<oldrub because we really want rub-rstr<=oldrub-oldstr
	assert (all(rlb >= oldrlb && rub < oldrub));
      }
      vector<ibbox> bbs(nprocs);
      for (int c=0; c<nprocs; ++c) {
	ivect cstr = rstr;
	ivect clb = rlb;
	ivect cub = rub;
	// split the components along the z axis
	const int glonpz = (rub[dim-1] - rlb[dim-1]) / cstr[dim-1];
	const int locnpz = (glonpz + nprocs - 1) / nprocs;
	const int zstep = locnpz * cstr[dim-1];
	clb[dim-1] = rlb[dim-1] + zstep *  c;
	cub[dim-1] = rlb[dim-1] + zstep * (c+1);
 	if (clb[dim-1] > rub[dim-1]) clb[dim-1] = rub[dim-1];
 	if (cub[dim-1] > rub[dim-1]) cub[dim-1] = rub[dim-1];
	assert (clb[dim-1] <= cub[dim-1]);
	assert (cub[dim-1] <= rub[dim-1]);
	bbs[c] = ibbox(clb, cub-cstr, cstr);
      }
//       bbss[rl] = bbs;
      bbss[rl].clear();
      assert (Carpet::hh->components(rl) % nprocs == 0);
      for (int cc=0; cc<(int)Carpet::hh->components(rl)/nprocs; ++cc) {
	bbss[rl].insert(bbss[rl].end(), bbs.begin(), bbs.end());
      }
    }
    
    vector<vector<vector<ibbox> > > bbsss
      = hh->make_multigrid_boxes(bbss, mglevels);
    
    vector<vector<int> > pss(bbss.size());
    vector<vector<bvect> > obss(bbss.size());
    for (int rl=0; rl<reflevels; ++rl) {
      pss[rl] = vector<int>(bbss[rl].size());
      obss[rl] = vector<bvect>(bbss[rl].size());
      // make sure all processors have the same number of components
      assert (bbss[rl].size() % nprocs == 0);
      for (int c=0; c<(int)bbss[rl].size(); ++c) {
	// distribute among processors
	pss[rl][c] = c % nprocs;
	for (int d=0; d<dim; ++d) {
	  // assume the components are split along the z axis
	  obss[rl][c][d][0] = d<dim-1 || c==0;
	  obss[rl][c][d][1] = d<dim-1 || c==(int)bbss[rl].size()-1;
	}
      }
    }
    
    hh->recompose(bbsss, obss, pss);
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
  
  
  
  static void OutputGridStructure (const cGH * const cgh,
				   const gh<dim>::rexts& bbsss,
				   const gh<dim>::rbnds& obss,
				   const gh<dim>::rprocs& pss)
  {
    DECLARE_CCTK_PARAMETERS;
    
    // Output only on the root processor
    if (CCTK_MyProc(cgh) != 0) return;
    
    // Output only if output is desired
    if (strcmp(grid_structure_filename, "") == 0) return;
    
    ofstream file;
    char filename[1000];
    snprintf (filename, sizeof(filename), "%s/%s",
	      outdir, grid_structure_filename);
    
    static bool do_truncate = true;
    
    if (do_truncate) {
      do_truncate = false;
      struct stat fileinfo;
      if (! IOUtil_RestartFromRecovery(cgh)
	  || stat(filename, &fileinfo)!=0) {
	file.open (filename, ios::out | ios::trunc);
	assert (file.good());
	file << "# grid structure" << endl
	     << "# format: reflevel component mglevel   processor bounding-box is-outer-boundary" << endl;
	assert (file.good());
      }
    }
    if (! file.is_open()) {
      file.open (filename, ios::app);
      assert (file.good());
    }
    
    file << "iteration " << cgh->cctk_iteration << endl;
    file << "reflevels " << bbsss.size() << endl;
    for (int rl=0; rl<(int)bbsss.size(); ++rl) {
      file << rl << " components " << bbsss[rl].size() << endl;
      for (int c=0; c<(int)bbsss[rl].size(); ++c) {
	file << rl << " " << c << " mglevels " << bbsss[rl][c].size() << endl;
	for (int ml=0; ml<(int)bbsss[rl][c].size(); ++ml) {
	  file << rl << " " << c << " " << ml << "   " << pss[rl][c] << " " << bbsss[rl][c][ml] << obss[rl][c] << endl;
	}
      }
    }
    file << endl;
    
    file.close();
    assert (file.good());
  }
  
  
  
  void SplitRegions (const cGH* cgh, vector<ibbox>& bbs, vector<bvect>& obs)
  {
    DECLARE_CCTK_PARAMETERS;
    
    if (CCTK_EQUALS (processor_topology, "automatic")) {
      SplitRegions_AlongZ (cgh, bbs, obs);
    } else if (CCTK_EQUALS (processor_topology, "manual")) {
      SplitRegions_AsSpecified (cgh, bbs, obs);
    } else {
      abort();
    }
  }
  
  
  
  void SplitRegions_AlongZ (const cGH* cgh, vector<ibbox>& bbs,
			    vector<bvect>& obs)
  {
    // Something to do?
    if (bbs.size() == 0) return;
    
    const int nprocs = CCTK_nProcs(cgh);
    
    if (nprocs==1) return;
    
    assert (bbs.size() == 1);
    
    const ivect rstr = bbs[0].stride();
    const ivect rlb  = bbs[0].lower();
    const ivect rub  = bbs[0].upper() + rstr;
    const bvect obnd = obs[0];
    
    bbs.resize(nprocs);
    obs.resize(nprocs);
    for (int c=0; c<nprocs; ++c) {
      ivect cstr = rstr;
      ivect clb = rlb;
      ivect cub = rub;
      const int glonpz = (rub[dim-1] - rlb[dim-1]) / cstr[dim-1];
      const int locnpz = (glonpz + nprocs - 1) / nprocs;
      const int zstep = locnpz * cstr[dim-1];
      clb[dim-1] = rlb[dim-1] + zstep *  c;
      cub[dim-1] = rlb[dim-1] + zstep * (c+1);
      if (clb[dim-1] > rub[dim-1]) clb[dim-1] = rub[dim-1];
      if (cub[dim-1] > rub[dim-1]) cub[dim-1] = rub[dim-1];
      assert (clb[dim-1] <= cub[dim-1]);
      assert (cub[dim-1] <= rub[dim-1]);
      bbs[c] = ibbox(clb, cub-cstr, cstr);
      obs[c] = obnd;
      if (c>0)        obs[c][dim-1][0] = false;
      if (c<nprocs-1) obs[c][dim-1][1] = false;
    }
  }
  
  
  
  void SplitRegions_AsSpecified (const cGH* cgh, vector<ibbox>& bbs,
				 vector<bvect>& obs)
  {
    DECLARE_CCTK_PARAMETERS;
    
    // Something to do?
    if (bbs.size() == 0) return;
    
    const int nprocs = CCTK_nProcs(cgh);
    
    if (nprocs==1) return;
    
    assert (bbs.size() == 1);
    
    const ivect rstr = bbs[0].stride();
    const ivect rlb  = bbs[0].lower();
    const ivect rub  = bbs[0].upper() + rstr;
    const bvect obnd = obs[0];
    
    const ivect nprocs_dir
      (processor_topology_3d_x, processor_topology_3d_y,
       processor_topology_3d_z);
    assert (all (nprocs_dir > 0));
    if (prod(nprocs_dir) != nprocs) {
      CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
		  "The specified processor topology [%d,%d,%d] does not fit the number of processors, which is %d", nprocs_dir[0], nprocs_dir[1], nprocs_dir[2], nprocs);
    }
    assert (prod(nprocs_dir) == nprocs);
    
    bbs.resize(nprocs);
    obs.resize(nprocs);
    assert (dim==3);
    for (int k=0; k<nprocs_dir[2]; ++k) {
      for (int j=0; j<nprocs_dir[1]; ++j) {
	for (int i=0; i<nprocs_dir[0]; ++i) {
	  const int c = i + nprocs_dir[0] * (j + nprocs_dir[1] * k);
	  const ivect ipos (i, j, k);
	  ivect cstr = rstr;
	  ivect clb = rlb;
	  ivect cub = rub;
	  const ivect glonp = (rub - rlb) / cstr;
	  const ivect locnp = (glonp + nprocs_dir - 1) / nprocs_dir;
	  const ivect step = locnp * cstr;
	  clb = rlb + step *  ipos;
	  cub = rlb + step * (ipos+1);
	  clb = min (clb, rub);
	  cub = min (cub, rub);
	  assert (all (clb <= cub));
	  assert (all (cub <= rub));
	  bbs[c] = ibbox(clb, cub-cstr, cstr);
	  obs[c] = obnd;
	  if (i>0) obs[c][0][0] = false;
	  if (j>0) obs[c][1][0] = false;
	  if (k>0) obs[c][2][0] = false;
	  if (i<nprocs_dir[0]-1) obs[c][0][1] = false;
	  if (j<nprocs_dir[1]-1) obs[c][1][1] = false;
	  if (k<nprocs_dir[2]-1) obs[c][2][1] = false;
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
