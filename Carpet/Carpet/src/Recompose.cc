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

#include "bbox.hh"
#include "bboxset.hh"
#include "defs.hh"
#include "gh.hh"
#include "vect.hh"

#include "carpet.hh"
#include "modes.hh"

extern "C" {
  static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/Recompose.cc,v 1.60 2004/03/23 19:32:59 schnetter Exp $";
  CCTK_FILEVERSION(Carpet_Carpet_Recompose_cc);
}

#define DEBUG false             // false or true



namespace Carpet {
  
  using namespace std;
  
  
  
  // Reduction operator
  template<typename iter, typename func>
  static typename func::result_type
  reduce (iter const first, iter const last,
          typename func::result_type const & init)
  {
    typename func::result_type res (init);
    for (iter it (first); it != last; ++it) {
      res = func::operator() (res, *it);
    }
    return res;
  }
  
  
  
  static void SplitRegions_Automatic_Recursively (bvect const & dims,
                                                  int const nprocs,
                                                  rvect const rshape,
                                                  ibbox const & bb,
                                                  bbvect const & ob,
                                                  int const & p,
                                                  vector<ibbox> & bbs,
                                                  vector<bbvect> & obs,
                                                  vector<int> & ps);
  static void SplitRegions_AsSpecified (const cGH* cgh,
                                        vector<ibbox>& bbs,
                                        vector<bbvect>& obs,
                                        vector<int>& ps);
  
  
  
  void CheckRegions (const gh<dim>::rexts & bbsss,
                     const gh<dim>::rbnds & obss,
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
      assert (bbsss.at(rl).size() > 0);
      for (int c=0; c<(int)bbsss.at(rl).size(); ++c) {
        // At least one multigrid level
        assert (bbsss.at(rl).at(c).size() > 0);
        for (int ml=0; ml<(int)bbsss.at(rl).at(c).size(); ++ml) {
          // Check sizes
          // Do allow processors with zero grid points
//           assert (all(bbsss.at(rl).at(c).at(ml).lower() <= bbsssi.at(rl).at(c).at(ml).upper()));
          // Check strides
          const int str = ipow(reffact, maxreflevels-rl-1) * ipow(mgfact, ml);
          assert (all(bbsss.at(rl).at(c).at(ml).stride() == str));
          // Check alignments
          assert (all(bbsss.at(rl).at(c).at(ml).lower() % str == 0));
          assert (all(bbsss.at(rl).at(c).at(ml).upper() % str == 0));
        }
      }
    }
    
    assert (pss.size() == bbsss.size());
    assert (obss.size() == bbsss.size());
    for (int rl=0; rl<(int)bbsss.size(); ++rl) {
      assert (obss.at(rl).size() == bbsss.at(rl).size());
      assert (pss.at(rl).size() == bbsss.at(rl).size());
    }
    
  }
  
  
  
  void Regrid (const cGH* cgh, const int rl,
               const int initialise_from, const bool do_prolongate)
  {
    assert (is_meta_mode());
    
    if (! CCTK_IsFunctionAliased ("Carpet_Regrid")) {
      static bool didtell = false;
      if (!didtell) {
	CCTK_WARN (1, "No regridding routine has been provided.  There will be no regridding.  Maybe you forgot to activate a regridding thorn?");
	didtell = true;
      }
      return;
    }
    
    for (int m=0; m<maps; ++m) {
      
      jjvect nboundaryzones, is_internal, is_staggered, shiftout;
      CCTK_INT const ierr = GetBoundarySpecification
        (2*dim, &nboundaryzones[0][0], &is_internal[0][0],
         &is_staggered[0][0], &shiftout[0][0]);
      assert (!ierr);
      
      gh<dim>::rexts  bbsss = vhh.at(m)->extents;
      gh<dim>::rbnds  obss  = vhh.at(m)->outer_boundaries;
      gh<dim>::rprocs pss   = vhh.at(m)->processors;
      
      // Check whether to recompose
      CCTK_INT const do_recompose = Carpet_Regrid
        (cgh, rl, m,
         2*dim, &nboundaryzones[0][0], &is_internal[0][0],
         &is_staggered[0][0], &shiftout[0][0],
         &bbsss, &obss, &pss);
      assert (do_recompose >= 0);
      
      if (do_recompose) {
        
        // Check the regions
        CheckRegions (bbsss, obss, pss);
        // TODO: check also that the current and all coarser levels
        // did not change
        
        // Write grid structure to file
        OutputGridStructure (cgh, m, bbsss, obss, pss);
        
        // Recompose
        vhh.at(m)->recompose (bbsss, obss, pss,
                              initialise_from, do_prolongate);
        
        OutputGrids (cgh, m, *vhh.at(m));
        
      }
      
    } // for m
    
    // Calculate new number of levels
    reflevels = vhh.at(0)->reflevels();
    for (int m=0; m<maps; ++m) {
      assert (vhh.at(m)->reflevels() == reflevels);
    }
    
    // One cannot switch off the current level
    assert (reflevels>rl);
  }
  
  
  
  void OutputGrids (const cGH* cgh, const int m, const gh<dim>& hh)
  {
    DECLARE_CCTK_PARAMETERS;
    
    if (verbose) {
      CCTK_INFO ("New bounding boxes:");
      for (int rl=0; rl<hh.reflevels(); ++rl) {
        for (int c=0; c<hh.components(rl); ++c) {
          for (int ml=0; ml<hh.mglevels(rl,c); ++ml) {
            cout << "   m " << m << "   rl " << rl << "   c " << c
                 << "   ml " << ml
                 << "   bbox " << hh.extents.at(rl).at(c).at(ml)
                 << endl;
          }
	}
      }
      CCTK_INFO ("New processor distribution:");
      for (int rl=0; rl<hh.reflevels(); ++rl) {
        for (int c=0; c<hh.components(rl); ++c) {
          cout << "   m " << m << "   rl " << rl << "   c " << c
               << "   processor " << hh.processors.at(rl).at(c) << endl;
	}
      }
      cout << endl;
    }
    
    CCTK_INFO ("New grid structure (grid points):");
    for (int rl=0; rl<hh.reflevels(); ++rl) {
      for (int c=0; c<hh.components(rl); ++c) {
        for (int ml=0; ml<hh.mglevels(rl,c); ++ml) {
          const int convfact = ipow(mgfact, ml);
          const int levfact = ipow(reffact, rl);
          const ivect lower = hh.extents.at(rl).at(c).at(ml).lower();
          const ivect upper = hh.extents.at(rl).at(c).at(ml).upper();
          assert (all(lower * levfact % maxreflevelfact == 0));
          assert (all(upper * levfact % maxreflevelfact == 0));
          assert (all(((upper - lower) * levfact / maxreflevelfact)
                      % convfact == 0));
          cout << "   [" << ml << "][" << rl << "][" << m << "][" << c << "]"
               << "   exterior extent: "
               << lower * levfact / maxreflevelfact
               << " : "
               << upper * levfact / maxreflevelfact
               << "   ("
               << (upper - lower) * levfact / maxreflevelfact / convfact + 1
               << ")" << endl;
        }
      }
    }
    
    CCTK_INFO ("New grid structure (coordinates):");
    for (int rl=0; rl<hh.reflevels(); ++rl) {
      for (int c=0; c<hh.components(rl); ++c) {
        for (int ml=0; ml<hh.mglevels(rl,c); ++ml) {
          const rvect origin = origin_space.at(0);
          const rvect delta = delta_space;
          const ivect lower = hh.extents.at(rl).at(c).at(ml).lower();
          const ivect upper = hh.extents.at(rl).at(c).at(ml).upper();
          const int convfact = ipow(mgfact, ml);
          const int levfact = ipow(reffact, rl);
          cout << "   [" << ml << "][" << rl << "][" << m << "][" << c << "]"
               << "   exterior extent: "
               << origin + delta * lower / maxreflevelfact
               << " : "
               << origin + delta * upper / maxreflevelfact
               << " : "
               << delta * convfact / levfact << endl;
        }
      }
    }
  }
  
  
  
  void OutputGridStructure (const cGH * const cgh,
                            const int m,
                            const gh<dim>::rexts & bbsss,
                            const gh<dim>::rbnds & obss,
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
	     << "# format: map reflevel component mglevel   processor bounding-box is-outer-boundary" << endl;
	assert (file.good());
      }
    }
    if (! file.is_open()) {
      file.open (filename, ios::app);
      assert (file.good());
    }
    
    file << "iteration " << cgh->cctk_iteration << endl;
    file << "maps " << maps << endl;
    file << m << " reflevels " << bbsss.size() << endl;
    for (size_t rl=0; rl<bbsss.size(); ++rl) {
      file << m << " " << rl << " components " << bbsss.at(rl).size() << endl;
      for (size_t c=0; c<bbsss.at(rl).size(); ++c) {
        file << m << " " << rl << " " << c << " mglevels " << bbsss.at(rl).at(c).size() << endl;
        for (size_t ml=0; ml<bbsss.at(rl).at(c).size(); ++ml) {
          file << m << " " << rl << " " << c << " " << ml << "   " << pss.at(rl).at(c) << " " << bbsss.at(rl).at(c).at(ml) << obss.at(rl).at(c) << endl;
        }
      }
    }
    file << endl;
    
    file.close();
    assert (file.good());
  }
  
  
  
  // TODO: this routine should go into CarpetRegrid (except maybe
  // SplitRegions_AlongZ for grid arrays)
  void SplitRegions (const cGH* cgh, vector<ibbox>& bbs, vector<bbvect>& obs,
                     vector<int>& ps)
  {
    DECLARE_CCTK_PARAMETERS;
    
    if (CCTK_EQUALS (processor_topology, "along-z")) {
      SplitRegions_AlongZ (cgh, bbs, obs, ps);
    } else if (CCTK_EQUALS (processor_topology, "along-dir")) {
      SplitRegions_AlongDir (cgh, bbs, obs, ps, split_direction);
    } else if (CCTK_EQUALS (processor_topology, "automatic")) {
      SplitRegions_Automatic (cgh, bbs, obs, ps);
    } else if (CCTK_EQUALS (processor_topology, "manual")) {
      SplitRegions_AsSpecified (cgh, bbs, obs, ps);
    } else {
      assert (0);
    }
  }
  
  
  
  void SplitRegions_AlongZ (const cGH* cgh, vector<ibbox>& bbs,
			    vector<bbvect>& obs, vector<int>& ps)
  {
    SplitRegions_AlongDir (cgh, bbs, obs, ps, 2);
  }
  
  
  
  void SplitRegions_AlongDir (const cGH* cgh, vector<ibbox>& bbs,
                              vector<bbvect>& obs, vector<int>& ps,
                              const int dir)
  {
    // Something to do?
    if (bbs.size() == 0) {
      ps.resize(0);
      return;
    }
    
    const int nprocs = CCTK_nProcs(cgh);
    
    if (nprocs==1) {
      ps.resize(1);
      ps.at(0) = 0;
      return;
    }
    
    assert (bbs.size() == 1);
    
    assert (dir>=0 && dir<dim);
    
    const ivect rstr = bbs.at(0).stride();
    const ivect rlb  = bbs.at(0).lower();
    const ivect rub  = bbs.at(0).upper() + rstr;
    const bbvect obnd = obs.at(0);
    
    bbs.resize(nprocs);
    obs.resize(nprocs);
    ps.resize(nprocs);
    for (int c=0; c<nprocs; ++c) {
      ivect cstr = rstr;
      ivect clb = rlb;
      ivect cub = rub;
      const int glonpz = (rub[dir] - rlb[dir]) / cstr[dir];
      const int locnpz = (glonpz + nprocs - 1) / nprocs;
      const int zstep = locnpz * cstr[dir];
      clb[dir] = rlb[dir] + zstep *  c;
      cub[dir] = rlb[dir] + zstep * (c+1);
      if (clb[dir] > rub[dir]) clb[dir] = rub[dir];
      if (cub[dir] > rub[dir]) cub[dir] = rub[dir];
      assert (clb[dir] <= cub[dir]);
      assert (cub[dir] <= rub[dir]);
      bbs.at(c) = ibbox(clb, cub-cstr, cstr);
      obs.at(c) = obnd;
      ps.at(c) = c;
      if (c>0)        obs.at(c)[dir][0] = false;
      if (c<nprocs-1) obs.at(c)[dir][1] = false;
    }
    
    for (size_t n=0; n<ps.size(); ++n) {
      assert (ps.at(n) == n);
    }
  }
  
  
  
  static void SplitRegions_Automatic_Recursively (bvect const & dims,
                                                  int const nprocs,
                                                  rvect const rshape,
                                                  ibbox const & bb,
                                                  bbvect const & ob,
                                                  int const & p,
                                                  vector<ibbox> & bbs,
                                                  vector<bbvect> & obs,
                                                  vector<int> & ps)
  {
    if (DEBUG) cout << "SRAR enter" << endl;
    // check preconditions
    assert (nprocs >= 1);
    
    // are we done?
    if (all(dims)) {
      if (DEBUG) cout << "SRAR bottom" << endl;
      
      // check precondition
      assert (nprocs == 1);
      
      // return arguments
      bbs.assign (1, bb);
      obs.assign (1, ob);
      ps.assign (1, p);
      
      // return
      if (DEBUG) cout << "SRAR exit" << endl;
      return;
    }
    
    // choose a direction
    int mydim = -1;
    CCTK_REAL mysize = 0;
    int alldims = 0;
    CCTK_REAL allsizes = 1;
    for (int d=0; d<dim; ++d) {
      if (! dims[d]) {
        ++ alldims;
        allsizes *= rshape[d];
        if (rshape[d] >= mysize) {
          mydim = d;
          mysize = rshape[d];
        }
      }
    }
    assert (mydim>=0 && mydim<dim);
    assert (mysize>=0);
    if (DEBUG) cout << "SRAR mydim " << mydim << endl;
    if (DEBUG) cout << "SRAR mysize " << mysize << endl;
    
    if (mysize == 0) {
      // the bbox is empty
      if (DEBUG) cout << "SRAR empty" << endl;
      
      // create the bboxes
      bbs.clear();
      obs.clear();
      ps.clear();
      bbs.reserve(nprocs);
      obs.reserve(nprocs);
      ps.reserve(nprocs);
      
      // create a new bbox
      assert (bb.empty());
      bbvect const newob (false);
      ibbox const newbb (bb);
      int const newp (p);
      if (DEBUG) cout << "SRAR " << mydim << " newbb " << newbb << endl;
      if (DEBUG) cout << "SRAR " << mydim << " newob " << newob << endl;
      if (DEBUG) cout << "SRAR " << mydim << " newp " << newp << endl;
      
      // store
      bbs.insert (bbs.end(), nprocs, newbb);
      obs.insert (obs.end(), nprocs, newob);
      for (int p=0; p<nprocs; ++p) ps.insert (ps.end(), 1, newp+p);
      
      // check postconditions
      assert (bbs.size() == nprocs);
      assert (obs.size() == nprocs);
      assert (ps.size() == nprocs);
      if (DEBUG) cout << "SRAR exit" << endl;
      return;
    }
    
    // mark this direction as done
    assert (! dims[mydim]);
    bvect const newdims = dims.replace(mydim, true);
    
    // choose a number of slices for this direction
    int const nslices = min(nprocs, (int)floor(mysize * pow(nprocs/allsizes, 1.0/alldims) + 0.5));
    assert (nslices <= nprocs);
    if (DEBUG) cout << "SRAR " << mydim << " nprocs " << nprocs << endl;
    if (DEBUG) cout << "SRAR " << mydim << " nslices " << nslices << endl;
    
    // split the remaining processors
    vector<int> mynprocs(nslices);
    int const mynprocs_base = nprocs / nslices;
    int const mynprocs_left = nprocs - nslices * mynprocs_base;
    for (int n=0; n<nslices; ++n) {
      mynprocs.at(n) = n < mynprocs_left ? mynprocs_base+1 : mynprocs_base;
    }
    int sum_mynprocs = 0;
    for (int n=0; n<nslices; ++n) {
      sum_mynprocs += mynprocs.at(n);
    }
    assert (sum_mynprocs == nprocs);
    if (DEBUG) cout << "SRAR " << mydim << " mynprocs " << mynprocs << endl;
    
    // split the region
    vector<int> myslice(nslices);
    int slice_left = ((bb.upper() - bb.lower()) / bb.stride())[mydim] + 1;
    int nprocs_left = nprocs;
    for (int n=0; n<nslices; ++n) {
      if (n == nslices-1) {
        myslice.at(n) = slice_left;
      } else {
        myslice.at(n) = (int)floor(1.0 * slice_left * mynprocs.at(n) / nprocs_left + 0.5);
      }
      assert (myslice.at(n) >= 0);
      slice_left -= myslice.at(n);
      nprocs_left -= mynprocs.at(n);
    }
    assert (slice_left == 0);
    assert (nprocs_left == 0);
    if (DEBUG) cout << "SRAR " << mydim << " myslice " << myslice << endl;
    
    // create the bboxes and recurse
    if (DEBUG) cout << "SRAR " << mydim << ": create bboxes" << endl;
    bbs.clear();
    obs.clear();
    ps.clear();
    bbs.reserve(nprocs);
    obs.reserve(nprocs);
    ps.reserve(nprocs);
    ivect last_up;
    for (int n=0; n<nslices; ++n) {
      if (DEBUG) cout << "SRAR " << mydim << " n " << n << endl;
      
      // create a new bbox
      ivect lo = bb.lower();
      ivect up = bb.upper();
      ivect str = bb.stride();
      bbvect newob = ob;
      if (n > 0) {
        lo[mydim] = last_up[mydim] + str[mydim];
        newob[mydim][0] = false;
      }
      if (n < nslices-1) {
        up[mydim] = lo[mydim] + (myslice.at(n)-1) * str[mydim];
        newob[mydim][1] = false;
        last_up = up;
      }
      ibbox newbb(lo, up, str);
      int newp(p + n * mynprocs_base + (n < mynprocs_left ? n : mynprocs_left));
      if (DEBUG) cout << "SRAR " << mydim << " newbb " << newbb << endl;
      if (DEBUG) cout << "SRAR " << mydim << " newob " << newob << endl;
      if (DEBUG) cout << "SRAR " << mydim << " newp " << newp << endl;
      
      // recurse
      vector<ibbox> newbbs;
      vector<bbvect> newobs;
      vector<int> newps;
      SplitRegions_Automatic_Recursively
        (newdims, mynprocs.at(n), rshape,
         newbb, newob, newp, newbbs, newobs, newps);
      if (DEBUG) cout << "SRAR " << mydim << " newbbs " << newbbs << endl;
      if (DEBUG) cout << "SRAR " << mydim << " newobs " << newobs << endl;
      if (DEBUG) cout << "SRAR " << mydim << " newps " << newps << endl;
      
      // store
      assert (newbbs.size() == mynprocs.at(n));
      assert (newobs.size() == mynprocs.at(n));
      assert (newps.size() == mynprocs.at(n));
      bbs.insert (bbs.end(), newbbs.begin(), newbbs.end());
      obs.insert (obs.end(), newobs.begin(), newobs.end());
      ps.insert (ps.end(), newps.begin(), newps.end());
    }
    
    // check postconditions
    assert (bbs.size() == nprocs);
    assert (obs.size() == nprocs);
    assert (ps.size() == nprocs);
    for (size_t n=0; n<ps.size(); ++n) {
      assert (ps.at(n) == p+n);
    }
    if (DEBUG) cout << "SRAR exit" << endl;
  }
  
  
  
  void SplitRegions_Automatic (const cGH* cgh, vector<ibbox>& bbs,
                               vector<bbvect>& obs, vector<int>& ps)
  {
    if (DEBUG) cout << "SRA enter" << endl;
    // Something to do?
    if (bbs.size() == 0) {
      ps.resize(0);
      return;
    }
    
    const int nprocs = CCTK_nProcs(cgh);
    if (DEBUG) cout << "SRA nprocs " << nprocs << endl;
    
//     if (nprocs==1) return;
    
    // split the processors
    // nslices: number of disjoint bboxes
    int const nslices = bbs.size();
    if (DEBUG) cout << "SRA nslices " << nslices << endl;
    // ncomps: number of components per processor
    int const ncomps = (nslices + nprocs - 1) / nprocs;
    if (DEBUG) cout << "SRA ncomps " << ncomps << endl;
    assert (ncomps > 0);
    vector<int> mysize(nslices);
    for (int c=0; c<nslices; ++c) {
      mysize.at(c) = bbs.at(c).size();
    }
    vector<int> mynprocs(nslices);
    {
      if (DEBUG) cout << "SRA: distributing processors to slices" << endl;
      int ncomps_left = nprocs * ncomps;
      for (int c=0; c<nslices; ++c) {
        mynprocs.at(c) = 1;
        -- ncomps_left;
      }
      while (ncomps_left > 0) {
        if (DEBUG) cout << "SRA ncomps_left " << ncomps_left << endl;
        int maxc = -1;
        CCTK_REAL maxratio = -1;
        for (int c=0; c<nslices; ++c) {
          CCTK_REAL const ratio = (CCTK_REAL)mysize.at(c) / mynprocs.at(c);
          if (ratio > maxratio) { maxc=c; maxratio=ratio; }
        }
        assert (maxc>=0 && maxc<nslices);
        ++ mynprocs.at(maxc);
        if (DEBUG) cout << "SRA maxc " << maxc << endl;
        if (DEBUG) cout << "SRA mynprocs[maxc] " << mynprocs.at(maxc) << endl;
        -- ncomps_left;
      }
      assert (ncomps_left == 0);
    }
    if (DEBUG) cout << "SRA mynprocs " << mynprocs << endl;
    
    vector<ibbox> allbbs;
    vector<bbvect> allobs;
    vector<int> allps;
    
    if (DEBUG) cout << "SRA: splitting regions" << endl;
    for (int c=0; c<nslices; ++c) {
      
      const ibbox bb = bbs.at(c);
      const bbvect ob = obs.at(c);
      int p = 0;
      for (int cc=0; cc<c; ++cc) {
        p += mynprocs.at(cc);
      }
      assert (p>=0 && p<=nprocs);
      if (DEBUG) cout << "SRA c " << c << endl;
      if (DEBUG) cout << "SRA bb " << bb << endl;
      if (DEBUG) cout << "SRA ob " << ob << endl;
      
      const ivect rstr = bb.stride();
      const ivect rlb  = bb.lower();
      const ivect rub  = bb.upper() + rstr;
    
      // calculate real shape factors
      rvect rshape;
      if (any(rub == rlb)) {
        // the bbox is empty
        rshape = 0.0;
      } else {
        for (int d=0; d<dim; ++d) {
          rshape[d] = (CCTK_REAL)(rub[d]-rlb[d]) / (rub[0]-rlb[0]);
        }
        const CCTK_REAL rfact = pow(nprocs / prod(rshape), 1.0/dim);
        rshape *= rfact;
        assert (abs(prod(rshape) - nprocs) < 1e-6);
      }
      if (DEBUG) cout << "SRA shapes " << rshape << endl;
      
      bvect const dims = false;
      
      vector<ibbox> thebbs;
      vector<bbvect> theobs;
      vector<int> theps;
      
      SplitRegions_Automatic_Recursively
        (dims, mynprocs.at(c), rshape, bb, ob, p, thebbs, theobs, theps);
      if (DEBUG) cout << "SRA thebbs " << thebbs << endl;
      if (DEBUG) cout << "SRA theobs " << theobs << endl;
      if (DEBUG) cout << "SRA theps " << theps << endl;
      
      allbbs.insert(allbbs.end(), thebbs.begin(), thebbs.end());
      allobs.insert(allobs.end(), theobs.begin(), theobs.end());
      allps.insert(allps.end(), theps.begin(), theps.end());
      
    } // for c
    
    bbs = allbbs;
    obs = allobs;
    ps = allps;
    for (size_t n=0; n<ps.size(); ++n) {
      assert (ps.at(n) == n);
    }
    
    if (DEBUG) cout << "SRA exit" << endl;
  }
  
  
  
  static void SplitRegions_AsSpecified (const cGH* cgh,
                                        vector<ibbox>& bbs,
                                        vector<bbvect>& obs,
                                        vector<int>& ps)
  {
    DECLARE_CCTK_PARAMETERS;
    
    // Something to do?
    if (bbs.size() == 0) {
      ps.resize(0);
      return;
    }
    
    const int nprocs = CCTK_nProcs(cgh);
    
    assert (bbs.size() == 1);
    
    const ivect rstr = bbs.at(0).stride();
    const ivect rlb  = bbs.at(0).lower();
    const ivect rub  = bbs.at(0).upper() + rstr;
    const bbvect obnd = obs.at(0);
    
    const ivect nprocs_dir
      (processor_topology_3d_x, processor_topology_3d_y,
       processor_topology_3d_z);
    assert (all (nprocs_dir > 0));
    if (prod(nprocs_dir) != nprocs) {
      CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
		  "The specified processor topology [%d,%d,%d] requires %d processors, but there are %d processors", nprocs_dir[0], nprocs_dir[1], nprocs_dir[2], prod(nprocs_dir), nprocs);
    }
    assert (prod(nprocs_dir) == nprocs);
    
    bbs.resize(nprocs);
    obs.resize(nprocs);
    ps.resize(nprocs);
    const ivect cstr = rstr;
    const ivect glonp = (rub - rlb) / cstr;
//     const ivect locnp = (glonp + nprocs_dir - 1) / nprocs_dir;
    const ivect locnp = glonp / nprocs_dir;
    const ivect rem = glonp % nprocs_dir;
    const ivect step = locnp * cstr;
    assert (dim==3);
    for (int k=0; k<nprocs_dir[2]; ++k) {
      for (int j=0; j<nprocs_dir[1]; ++j) {
	for (int i=0; i<nprocs_dir[0]; ++i) {
	  const int c = i + nprocs_dir[0] * (j + nprocs_dir[1] * k);
	  const ivect ipos (i, j, k);
	  ivect clb = rlb + step *  ipos;
	  ivect cub = rlb + step * (ipos+1);
// 	  clb = min (clb, rub);
// 	  cub = min (cub, rub);
          for (int d=0; d<dim; ++d) {
            if (ipos[d]<rem[d]) {
              clb[d] += cstr[d] * ipos[d];
              cub[d] += cstr[d] * (ipos[d]+1);
            } else {
              clb[d] += cstr[d] * rem[d];
              cub[d] += cstr[d] * rem[d];
            }
          }
	  assert (all (clb >= 0));
	  assert (all (clb <= cub));
	  assert (all (cub <= rub));
          assert (all (! (ipos==0) || clb==rlb));
          assert (all (! (ipos==nprocs_dir-1) || cub==rub));
	  bbs.at(c) = ibbox(clb, cub-cstr, cstr);
	  obs.at(c) = obnd;
          ps.at(c) = c;
	  if (i>0) obs.at(c)[0][0] = false;
	  if (j>0) obs.at(c)[1][0] = false;
	  if (k>0) obs.at(c)[2][0] = false;
	  if (i<nprocs_dir[0]-1) obs.at(c)[0][1] = false;
	  if (j<nprocs_dir[1]-1) obs.at(c)[1][1] = false;
	  if (k<nprocs_dir[2]-1) obs.at(c)[2][1] = false;
	}
      }
    }
    
    for (size_t n=0; n<ps.size(); ++n) {
      assert (ps.at(n) == n);
    }
  }
  
  
  
  static void MakeMultigridBoxes (const cGH* cgh,
                                  int const size,
                                  jjvect const & nboundaryzones,
                                  jjvect const & is_internal,
                                  jjvect const & is_staggered,
                                  jjvect const & shiftout,
                                  ibbox const & bb,
                                  bbvect const & ob,
                                  vector<ibbox>& bbs)
  {
    bbs.resize (mglevels);
    bbs.at(0) = bb;
    // boundary offsets
    assert (size==2*dim);
    // (distance in grid points between the exterior and the physical boundary)
    iivect offset;
    for (int d=0; d<dim; ++d) {
      for (int f=0; f<2; ++f) {
        assert (! is_staggered[d][f]);
        offset[d][f] = (+ (is_internal[d][f] ? 0 : nboundaryzones[d][f] - 1)
                        + shiftout[d][f]);
      }
    }
    for (int ml=1; ml<mglevels; ++ml) {
      // next finer grid
      ivect const flo = bbs.at(ml-1).lower();
      ivect const fhi = bbs.at(ml-1).upper();
      ivect const fstr = bbs.at(ml-1).stride();
      // this grid
      ivect const str = fstr * mgfact;
      ivect const modoffset = (xpose(offset)[0] - ivect(mgfact) * xpose(offset)[0] + ivect(mgfact) * xpose(offset)[0] * str) % str;
      ivect const lo = flo + xpose(ob)[0].ifthen (    xpose(offset)[0] - ivect(mgfact) * xpose(offset)[0] , modoffset);
      ivect const hi = fhi + xpose(ob)[1].ifthen ( - (xpose(offset)[1] - ivect(mgfact) * xpose(offset)[1]), modoffset + fstr - str);
      bbs[ml] = ibbox(lo,hi,str);
    }
  }
  
  void MakeMultigridBoxes (const cGH* cgh,
                           int const size,
                           jjvect const & nboundaryzones,
                           jjvect const & is_internal,
                           jjvect const & is_staggered,
                           jjvect const & shiftout,
                           vector<ibbox> const & bbs,
                           vector<bbvect> const & obs,
                           vector<vector<ibbox> >& bbss)
  {
    assert (bbs.size() == obs.size());
    bbss.resize(bbs.size());
    for (size_t c=0; c<bbs.size(); ++c) {
      MakeMultigridBoxes
        (cgh,
         size, nboundaryzones, is_internal, is_staggered, shiftout,
         bbs.at(c), obs.at(c), bbss.at(c));
    }
  }
  
} // namespace Carpet
