#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <list>
#include <sstream>
#include <string>
#include <vector>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "CactusBase/IOUtil/src/ioGH.h"

#include "Carpet/CarpetLib/src/bbox.hh"
#include "Carpet/CarpetLib/src/bboxset.hh"
#include "Carpet/CarpetLib/src/defs.hh"
#include "Carpet/CarpetLib/src/gh.hh"
#include "Carpet/CarpetLib/src/vect.hh"

#include "carpet.hh"

extern "C" {
  static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/Recompose.cc,v 1.39 2003/05/07 10:03:21 schnetter Exp $";
  CCTK_FILEVERSION(Carpet_Carpet_Recompose_cc);
}



namespace Carpet {
  
  using namespace std;
  
  
  
  typedef vect<bool,dim> bvect;
  
  typedef vect<int,dim> ivect;
  typedef bbox<int,dim> ibbox;
  
  typedef vect<double,dim> dvect;
  
  typedef vect<vect<bool,2>,dim> bbvect;
  
  
  
  static int (*regrid_routine) (const cGH * cckgGH,
				gh<dim>::rexts& bbsss,
				gh<dim>::rbnds& obss,
				gh<dim>::rprocs& pss) = 0;
  
  
  
  static void CheckRegions (const gh<dim>::rexts& bbsss,
			    const gh<dim>::rbnds& obss,
			    const gh<dim>::rprocs& pss);
  static void Adapt (const cGH* cgh, int reflevels, gh<dim>* hh);
  
  static void Output (const cGH* cgh, const gh<dim>* hh);
  
  static void OutputGridStructure (const cGH *cgh,
				   const gh<dim>::rexts& bbsss,
				   const gh<dim>::rbnds& obss,
				   const gh<dim>::rprocs& pss);
  
  
  
  void SplitRegions_AlongZ              (const cGH* cgh,
					 vector<ibbox>& bbs,
					 vector<bbvect>& obs);
  static void SplitRegions_Automatic_Recursively (bvect const & dims,
                                                  int const nprocs,
                                                  dvect const dshape,
                                                  ibbox const & bb,
                                                  bbvect const & ob,
                                                  vector<ibbox> & bbs,
                                                  vector<bbvect> & obs);
  static void SplitRegions_Automatic    (const cGH* cgh,
					 vector<ibbox>& bbs,
					 vector<bbvect>& obs);
  static void SplitRegions_AsSpecified  (const cGH* cgh,
					 vector<ibbox>& bbs,
					 vector<bbvect>& obs);
  
  static void MakeProcessors_RoundRobin (const cGH* cgh,
					 const gh<dim>::rexts& bbss,
					 gh<dim>::rprocs& pss);
  
  
  
  void CheckRegions (const gh<dim>::rexts& bbsss, const gh<dim>::rbnds& obss,
		     const gh<dim>::rprocs& pss)
  {
    // At least one level
    if (bbsss.size() == 0) {
      CCTK_WARN (0, "I cannot set up a grid hierarchy with zero refinement levels.");
    }
    assert (bbsss.size() > 0);
    // At most maxreflevels levels
    if ((int)bbsss.size() > maxreflevels) {
      CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
		  "I cannot set up a grid hierarchy with more than Carpet::max_refinement_levels refinement levels.  I found Carpet::max_refinement_levels=%d, while %d levels were requested.",
		  (int)maxreflevels, (int)bbsss.size());
    }
    assert ((int)bbsss.size() <= maxreflevels);
    for (int rl=0; rl<(int)bbsss.size(); ++rl) {
      // No empty levels
      assert (bbsss[rl].size() > 0);
      for (int c=0; c<(int)bbsss[rl].size(); ++c) {
	// At least one multigrid level
	assert (bbsss[rl][c].size() > 0);
	for (int ml=0; ml<(int)bbsss[rl][c].size(); ++ml) {
	  // Check sizes
	  // Do allow processors with zero grid points
// 	  assert (all(bbsss[rl][c][ml].lower() <= bbsss[rl][c][ml].upper()));
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
  
  
  
  void Regrid (const cGH* cgh, const int initialise_upto)
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
    Recompose (cgh, bbsss, obss, pss, initialise_upto);
  }
  
  
  
  void Recompose (const cGH* const cgh,
		  const gh<dim>::rexts& bbsss,
		  const gh<dim>::rbnds& obss,
		  const gh<dim>::rprocs& pss,
                  const int initialise_upto)
  {
    assert (mglevel == -1);
    assert (component == -1);
    
    // Check the regions
    CheckRegions (bbsss, obss, pss);
    
    // Write grid structure to file
    OutputGridStructure (cgh, bbsss, obss, pss);
    
    // Recompose
    hh->recompose (bbsss, obss, pss, initialise_upto);
    Output (cgh, hh);
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
    vector<vector<bbvect> > obss(bbss.size());
    for (int rl=0; rl<reflevels; ++rl) {
      pss[rl] = vector<int>(bbss[rl].size());
      obss[rl] = vector<bbvect>(bbss[rl].size());
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
  
  
  
  static void Output (const cGH* cgh, const gh<dim>* hh)
  {
    DECLARE_CCTK_PARAMETERS;
    
    if (verbose) {
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
    
    // Get grid hierarchy extentsion from IOUtil
    const ioGH * const iogh = (const ioGH *)CCTK_GHExtension (cgh, "IO");
    
    // Output only if IO exists and has been initialised
    if (! iogh) return;
    
    assert (iogh);
    
    // Create the output directory
    CCTK_CreateDirectory (0755, out_dir);
    
    ostringstream filenamebuf;
    filenamebuf << out_dir << "/" << grid_structure_filename;
    // we need a persistent temporary here
    string filenamestr = filenamebuf.str();
    const char * filename = filenamestr.c_str();
    
    ofstream file;
    
    static bool do_truncate = true;
    
    if (do_truncate) {
      do_truncate = false;
      struct stat fileinfo;
      if (! iogh->recovered
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
  
  
  
  // TODO: this routine should go into CarpetRegrid (except maybe
  // SplitRegions_AlongZ for grid arrays)
  void SplitRegions (const cGH* cgh, vector<ibbox>& bbs, vector<bbvect>& obs)
  {
    DECLARE_CCTK_PARAMETERS;
    
    if (CCTK_EQUALS (processor_topology, "along-z")) {
      SplitRegions_AlongZ (cgh, bbs, obs);
    } else if (CCTK_EQUALS (processor_topology, "automatic")) {
      SplitRegions_Automatic (cgh, bbs, obs);
    } else if (CCTK_EQUALS (processor_topology, "manual")) {
      SplitRegions_AsSpecified (cgh, bbs, obs);
    } else {
      assert (0);
    }
  }
  
  
  
  void SplitRegions_AlongZ (const cGH* cgh, vector<ibbox>& bbs,
			    vector<bbvect>& obs)
  {
    // Something to do?
    if (bbs.size() == 0) return;
    
    const int nprocs = CCTK_nProcs(cgh);
    
    if (nprocs==1) return;
    
    assert (bbs.size() == 1);
    
    const ivect rstr = bbs[0].stride();
    const ivect rlb  = bbs[0].lower();
    const ivect rub  = bbs[0].upper() + rstr;
    const bbvect obnd = obs[0];
    
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
  
  
  
  void SplitRegions_Automatic_Recursively (bvect const & dims,
                                           int const nprocs,
                                           dvect const dshape,
                                           ibbox const & bb,
                                           bbvect const & ob,
                                           vector<ibbox> & bbs,
                                           vector<bbvect> & obs)
  {
    // check preconditions
    assert (nprocs >= 1);
    
    // are we done?
    if (all(dims)) {
      
      // check precondition
      assert (nprocs == 1);
      
      // return arguments
      bbs.assign (1, bb);
      obs.assign (1, ob);
      
      // return
      return;
    }
    
    // choose a direction
    int mydim = -1;
    double mysize = 0;
    int alldims = 0;
    double allsizes = 1;
    for (int d=0; d<dim; ++d) {
      if (! dims[d]) {
        ++ alldims;
        allsizes *= dshape[d];
        if (dshape[d] > mysize) {
          mydim = d;
          mysize = dshape[d];
        }
      }
    }
    assert (mydim>=0 && mydim<dim);
    assert (mysize>0);
    
    // mark this direction as done
    assert (! dims[mydim]);
    bvect const newdims = dims.replace(mydim, true);
    
    // choose a number of slices for this direction
    int const nslices = min(nprocs, (int)floor(mysize * pow(nprocs/allsizes, 1.0/alldims) + 0.5));
    assert (nslices <= nprocs);
    
    // split the remaining processors
    vector<int> mynprocs(nslices);
    int const mynprocs_base = nprocs / nslices;
    int const mynprocs_left = nprocs - nslices * mynprocs_base;
    for (int n=0; n<nslices; ++n) {
      mynprocs[n] = n < mynprocs_left ? mynprocs_base+1 : mynprocs_base;
    }
    int sum_mynprocs = 0;
    for (int n=0; n<nslices; ++n) {
      sum_mynprocs += mynprocs[n];
    }
    assert (sum_mynprocs == nprocs);
    
    // split the region
    vector<int> myslice(nslices);
    int slice_left = ((bb.upper() - bb.lower()) / bb.stride())[mydim] + 1;
    int nprocs_left = nprocs;
    for (int n=0; n<nslices; ++n) {
      if (n == nslices-1) {
        myslice[n] = slice_left;
      } else {
        myslice[n] = (int)floor(1.0 * slice_left * mynprocs[n] / nprocs_left + 0.5);
      }
      slice_left -= myslice[n];
      nprocs_left -= mynprocs[n];
    }
    assert (slice_left == 0);
    assert (nprocs_left == 0);
    
    // create the bboxes and recurse
    bbs.clear();
    obs.clear();
    bbs.reserve(nprocs);
    obs.reserve(nprocs);
    for (int n=0; n<nslices; ++n) {
      
      // create a new bbox
      ivect lo = bb.lower();
      ivect up = bb.upper();
      ivect str = bb.stride();
      bbvect newob = ob;
      if (n > 0) {
        lo = lo.replace(mydim, bbs[n-1].upper()[mydim] + str[mydim]);
        newob[mydim][0] = false;
      }
      if (n < nslices-1) {
        up = up.replace(mydim, lo[mydim] + (myslice[n]-1) * str[mydim]);
        newob[mydim][1] = false;
      }
      ibbox newbb(lo, up, str);
      
      // recurse
      vector<ibbox> newbbs;
      vector<bbvect> newobs;
      SplitRegions_Automatic_Recursively
        (newdims, mynprocs[n], dshape, newbb, newob, newbbs, newobs);
      
      // store
      assert (newbbs.size() == mynprocs[n]);
      assert (newobs.size() == mynprocs[n]);
      bbs.insert (bbs.end(), newbbs.begin(), newbbs.end());
      obs.insert (obs.end(), newobs.begin(), newobs.end());
    }
    
    // check postconditions
    assert (bbs.size() == nprocs);
    assert (obs.size() == nprocs);
  }
  
  void SplitRegions_Automatic (const cGH* cgh, vector<ibbox>& bbs,
                               vector<bbvect>& obs)
  {
    // Something to do?
    if (bbs.size() == 0) return;
    
    const int nprocs = CCTK_nProcs(cgh);
    
//     if (nprocs==1) return;
    
    // split the processors
    int const nslices = bbs.size();
    int const ncomps = (nslices + nprocs - 1) / nprocs;
    assert (ncomps > 0);
    vector<int> mysize(nslices);
    for (int c=0; c<nslices; ++c) {
      mysize[c] = bbs[c].num_points();
    }
    vector<int> mynprocs(nslices);
    {
      int ncomps_left = nprocs * ncomps;
      for (int c=0; c<nslices; ++c) {
        mynprocs[c] = 1;
        -- ncomps_left;
      }
      while (ncomps_left > 0) {
        int maxc = -1;
        double maxratio = 0;
        for (int c=0; c<nslices; ++c) {
          double const ratio = (double)mysize[c] / mynprocs[c];
          if (ratio > maxratio) { maxc=c; maxratio=ratio; }
        }
        assert (maxc>=0 && maxc<nslices);
        ++ mynprocs[maxc];
        -- ncomps_left;
      }
      assert (ncomps_left == 0);
    }
    
    vector<ibbox> allbbs;
    vector<bbvect> allobs;
    
    for (int c=0; c<nslices; ++c) {
      
      const ibbox bb = bbs[c];
      const bbvect ob = obs[c];
      
      const ivect rstr = bb.stride();
      const ivect rlb  = bb.lower();
      const ivect rub  = bb.upper() + rstr;
    
      // calculate real shape factors
      dvect dshape;
      for (int d=0; d<dim; ++d) {
        dshape[d] = (double)(rub[d]-rlb[d]) / (rub[0]-rlb[0]);
      }
      const double dfact = pow(nprocs / prod(dshape), 1.0/dim);
      dshape *= dfact;
      assert (abs(prod(dshape) - nprocs) < 1e-6);
      
      bvect const dims = false;
      
      vector<ibbox> thebbs;
      vector<bbvect> theobs;

      SplitRegions_Automatic_Recursively
        (dims, mynprocs[c], dshape, bb, ob, thebbs, theobs);
      
      allbbs.insert(allbbs.end(), thebbs.begin(), thebbs.end());
      allobs.insert(allobs.end(), theobs.begin(), theobs.end());
      
    } // for c
    
    bbs = allbbs;
    obs = allobs;
  }
  
  
  
  void SplitRegions_AsSpecified (const cGH* cgh, vector<ibbox>& bbs,
				 vector<bbvect>& obs)
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
    const bbvect obnd = obs[0];
    
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
