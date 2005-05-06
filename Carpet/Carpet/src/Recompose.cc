#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <list>
#include <sstream>
#include <string>
#include <vector>

#include <sys/stat.h>
#include <sys/types.h>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "bbox.hh"
#include "bboxset.hh"
#include "defs.hh"
#include "gh.hh"
#include "vect.hh"

#include "carpet.hh"
#include "modes.hh"

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
  
  
  
  void CheckRegions (const gh::mexts & bbsss,
                     const gh::rbnds & obss,
                     const gh::rprocs& pss)
  {
    // At least one multigrid level
    if (bbsss.size() == 0) {
      CCTK_WARN (0, "I cannot set up a grid hierarchy with zero multigrid levels.");
    }
    assert (bbsss.size() > 0);
    // At least one refinement level
    if (bbsss.at(0).size() == 0) {
      CCTK_WARN (0, "I cannot set up a grid hierarchy with zero refinement levels.");
    }
    assert (bbsss.at(0).size() > 0);
    // At most maxreflevels levels
    if ((int)bbsss.at(0).size() > maxreflevels) {
      CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "I cannot set up a grid hierarchy with more than Carpet::max_refinement_levels refinement levels.  I found Carpet::max_refinement_levels=%d, while %d levels were requested.",
                  (int)maxreflevels, (int)bbsss.at(0).size());
    }
    assert ((int)bbsss.at(0).size() <= maxreflevels);
    for (int ml=0; ml<(int)bbsss.size(); ++ml) {
      for (int rl=0; rl<(int)bbsss.at(0).size(); ++rl) {
        // No empty levels
        assert (bbsss.at(ml).at(rl).size() > 0);
        for (int c=0; c<(int)bbsss.at(ml).at(rl).size(); ++c) {
          // Check sizes
          // Do allow processors with zero grid points
//           assert (all(bbsss.at(rl).at(c).at(ml).lower() <= bbsssi.at(rl).at(c).at(ml).upper()));
          // Check strides
          const ivect str
            = (maxspacereflevelfact / spacereffacts.at(rl) * ipow(mgfact, ml));
          assert (all(bbsss.at(ml).at(rl).at(c).stride() == str));
          // Check alignments
          assert (all(bbsss.at(ml).at(rl).at(c).lower() % str == 0));
          assert (all(bbsss.at(ml).at(rl).at(c).upper() % str == 0));
        }
      }
    }
    
    assert (pss.size() == bbsss.at(0).size());
    assert (obss.size() == bbsss.at(0).size());
    for (int rl=0; rl<(int)bbsss.at(0).size(); ++rl) {
      assert (obss.at(rl).size() == bbsss.at(0).at(rl).size());
      assert (pss.at(rl).size() == bbsss.at(0).at(rl).size());
    }
    
  }
  
  
 
  bool Regrid (const cGH* cgh, const bool force_recompose, const bool do_init)
  {
    DECLARE_CCTK_PARAMETERS;
    
    assert (is_level_mode());
    
    if (! CCTK_IsFunctionAliased ("Carpet_Regrid")) {
      static bool didtell = false;
      if (maxreflevels > 1 and ! didtell) {
	CCTK_WARN (2, "No regridding routine has been provided.  There will be no regridding.  Maybe you forgot to activate a regridding thorn?");
	didtell = true;
      }
      return false;
    }
    
    bool did_change = false;
    BEGIN_MAP_LOOP(cgh, CCTK_GF) {
      
      gh::mexts  bbsss = vhh.at(map)->extents();
      gh::rbnds  obss  = vhh.at(map)->outer_boundaries();
      gh::rprocs pss   = vhh.at(map)->processors();
      
      // Check whether to recompose
      CCTK_INT const do_recompose
        = Carpet_Regrid (cgh, &bbsss, &obss, &pss, force_recompose);
      assert (do_recompose >= 0);
      did_change = did_change or do_recompose;
      
      if (do_recompose) {
        
        CCTK_INFO ("Recomposing the grid hierarchy");
        
        // Check the regions
        CheckRegions (bbsss, obss, pss);
        // TODO: check also that the current and all coarser levels
        // did not change
        
        // Write grid structure to file
        OutputGridStructure (cgh, map, bbsss, obss, pss);
        
        // Recompose
        vhh.at(map)->recompose (bbsss, obss, pss, do_init);
        
        if (verbose) OutputGrids (cgh, map, *vhh.at(map));
        
      }

    } END_MAP_LOOP;
    
    // Calculate new number of levels
    reflevels = vhh.at(0)->reflevels();
    for (int m=0; m<maps; ++m) {
      assert (vhh.at(m)->reflevels() == reflevels);
    }
    
    // One cannot switch off the current level
    assert (reflevels > reflevel);
    
    return did_change;
  }
  
  
  
  void OutputGrids (const cGH* cgh, const int m, const gh& hh)
  {
    CCTK_INFO ("New grid structure (grid points):");
    for (int ml=0; ml<hh.mglevels(); ++ml) {
      for (int rl=0; rl<hh.reflevels(); ++rl) {
        for (int c=0; c<hh.components(rl); ++c) {
          const int convfact = ipow(mgfact, ml);
          const ivect levfact = spacereffacts.at(rl);
          const ivect lower = hh.extents().at(ml).at(rl).at(c).lower();
          const ivect upper = hh.extents().at(ml).at(rl).at(c).upper();
          const ivect stride = hh.extents().at(ml).at(rl).at(c).stride();
          assert (all(lower * levfact % maxspacereflevelfact == 0));
          assert (all(upper * levfact % maxspacereflevelfact == 0));
          assert (all(((upper - lower) * levfact / maxspacereflevelfact)
                      % convfact == 0));
          cout << "   [" << ml << "][" << rl << "][" << m << "][" << c << "]"
               << "   exterior: "
               << "proc "
               << hh.processors().at(rl).at(c)
               << "   "
//                << lower * levfact / maxreflevelfact
               << lower / stride
               << " : "
//                << upper * levfact / maxreflevelfact
               << upper / stride
               << "   ("
//                << (upper - lower) * levfact / maxreflevelfact / convfact + 1
               << (upper - lower) / stride + 1
               << ")"
               << endl;
        }
      }
    }
    
    CCTK_INFO ("New grid structure (coordinates):");
    for (int ml=0; ml<hh.mglevels(); ++ml) {
      for (int rl=0; rl<hh.reflevels(); ++rl) {
        for (int c=0; c<hh.components(rl); ++c) {
          const rvect origin = origin_space.at(0);
          const rvect delta = delta_space;
          const ivect lower = hh.extents().at(ml).at(rl).at(c).lower();
          const ivect upper = hh.extents().at(ml).at(rl).at(c).upper();
          const int convfact = ipow(mgfact, ml);
          const ivect levfact = spacereffacts.at(rl);
          cout << "   [" << ml << "][" << rl << "][" << m << "][" << c << "]"
               << "   exterior: "
               << origin + delta * lower / maxspacereflevelfact
               << " : "
               << origin + delta * upper / maxspacereflevelfact
               << " : "
               << delta * convfact / levfact << endl;
        }
      }
    }
    
    CCTK_INFO ("New grid statistics:");
    cout.setf (ios::fixed);
    const int oldprecision = cout.precision();
    for (int ml=0; ml<hh.mglevels(); ++ml) {
      CCTK_REAL coarsevolume = 0;
      for (int rl=0; rl<hh.reflevels(); ++rl) {
        
        const CCTK_REAL basevolume
          = prod (rvect (hh.baseextent.shape()) / hh.baseextent.stride());
        CCTK_REAL countvolume = 0;
        CCTK_REAL totalvolume = 0;
        CCTK_REAL totalvolume2 = 0;
        
        for (int c=0; c<hh.components(rl); ++c) {
          const ivect shape = hh.extents().at(ml).at(rl).at(c).shape();
          const ivect stride = hh.extents().at(ml).at(rl).at(c).stride();
          const CCTK_REAL volume = prod (rvect (shape) / stride);
          ++ countvolume;
          totalvolume += volume;
          totalvolume2 += ipow(volume, 2);
        }
        
        const CCTK_REAL avgvolume = totalvolume / countvolume;
        const CCTK_REAL stddevvolume
          = sqrt (max ((CCTK_REAL) 0.0,
                       totalvolume2 / countvolume - ipow (avgvolume, 2)));
        
        for (int c=0; c<hh.components(rl); ++c) {
          const ivect shape = hh.extents().at(ml).at(rl).at(c).shape();
          const ivect stride = hh.extents().at(ml).at(rl).at(c).stride();
          const CCTK_REAL volume = prod(rvect (shape) / stride);
          cout << "   [" << ml << "][" << rl << "][" << m << "][" << c << "]"
               << "   volume: " << setprecision(0) << volume
               << "   of parent: " << setprecision(1) << 100 * volume / totalvolume << "%"
               << "   of domain: " << setprecision(1) << 100 * volume / basevolume << "%"
               << endl;
        }
        
        cout << "   [" << ml << "][" << rl << "][" << m << "]"
             << "   average volume: " << setprecision(0) << avgvolume
             << "   of parent: " << setprecision(1) << 100 * avgvolume / totalvolume << "%"
             << "   of domain: " << setprecision(1) << 100 * avgvolume / basevolume << "%"
             << endl;
        cout << "   [" << ml << "][" << rl << "][" << m << "]"
             << "   standard deviation: " << setprecision(0) << stddevvolume
             << "   of parent: " << setprecision(1) << 100 * stddevvolume / totalvolume << "%"
             << "   of domain: " << setprecision(1) << 100 * stddevvolume / basevolume << "%"
             << endl;
        
        // TODO: ghost points vs. interior points (and boundary
        // points)
        
        CCTK_REAL countquadrupole = 0;
        CCTK_REAL minquadrupole = 1;
        CCTK_REAL totalquadrupole = 0;
        CCTK_REAL totalquadrupole2 = 0;
        
        for (int c=0; c<hh.components(rl); ++c) {
          const ivect shape = hh.extents().at(ml).at(rl).at(c).shape();
          const ivect stride = hh.extents().at(ml).at(rl).at(c).stride();
          const CCTK_REAL minlength = minval (rvect (shape) / stride);
          const CCTK_REAL maxlength = maxval (rvect (shape) / stride);
          const CCTK_REAL quadrupole = minlength / maxlength;
          ++ countquadrupole;
          minquadrupole = min (minquadrupole, quadrupole);
          totalquadrupole += quadrupole;
          totalquadrupole2 += ipow (quadrupole, 2);
          cout << "   [" << ml << "][" << rl << "][" << m << "][" << c << "]"
               << "   length ratio: " << setprecision(3) << quadrupole
               << endl;
        }
        
        const CCTK_REAL avgquadrupole = totalquadrupole / countquadrupole;
        const CCTK_REAL stddevquadrupole
          = sqrt (max ((CCTK_REAL) 0.0,
                       (totalquadrupole2 / countquadrupole
                        - ipow (avgquadrupole, 2))));
        
        cout << "   [" << ml << "][" << rl << "][" << m << "]"
             << "   average length ratio: " << setprecision(3) << avgquadrupole
             << "   standard deviation: " << setprecision(3) << stddevquadrupole
             << endl;
        
        // TODO: processor distribution, average load, std deviation
        
        coarsevolume = totalvolume * prod (rvect (spacereflevelfact));
      }
    }
    cout.precision (oldprecision);
    cout.unsetf (ios::fixed);
    
  }
  
  
  
  void OutputGridStructure (const cGH * const cgh,
                            const int m,
                            const gh::mexts & bbsss,
                            const gh::rbnds & obss,
                            const gh::rprocs& pss)
  {
    DECLARE_CCTK_PARAMETERS;
    
    // Output only on the root processor
    if (CCTK_MyProc(cgh) != 0) return;
    
    // Output only if output is desired
    if (strcmp(grid_structure_filename, "") == 0) return;
    
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
      if (IO_TruncateOutputFiles (cgh)
	  or stat(filename, &fileinfo)!=0) {
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
    file << m << " mglevels " << bbsss.size() << endl;
    for (int ml=0; ml<(int)bbsss.size(); ++ml) {
      file << m << " " << ml << " reflevels " << bbsss.at(ml).size() << endl;
      for (int rl=0; rl<(int)bbsss.at(ml).size(); ++rl) {
        file << m << " " << ml << " " << rl << " components " << bbsss.at(ml).at(rl).size() << endl;
        for (int c=0; c<(int)bbsss.at(ml).at(rl).size(); ++c) {
          file << m << " " << ml << " " << rl << " " << c << "   " << pss.at(rl).at(c) << " " << bbsss.at(ml).at(rl).at(c) << obss.at(rl).at(c) << endl;
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
    
    assert (dir>=0 and dir<dim);
    
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
    
    for (int n=0; n<(int)ps.size(); ++n) {
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
    assert (mydim>=0 and mydim<dim);
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
      for (int pp=0; pp<nprocs; ++pp) ps.insert (ps.end(), 1, newp+pp);
      
      // check postconditions
      assert ((int)bbs.size() == nprocs);
      assert ((int)obs.size() == nprocs);
      assert ((int)ps.size() == nprocs);
      if (DEBUG) cout << "SRAR exit" << endl;
      return;
    }
    
    // mark this direction as done
    assert (! dims[mydim]);
    bvect const newdims = dims.replace(mydim, true);
    
    // choose a number of slices for this direction
    int const nslices
      = min(nprocs,
            (int)floor(mysize * pow(nprocs/allsizes, (CCTK_REAL)1/alldims)
                       + 0.5));
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
      assert ((int)newbbs.size() == mynprocs.at(n));
      assert ((int)newobs.size() == mynprocs.at(n));
      assert ((int)newps.size() == mynprocs.at(n));
      bbs.insert (bbs.end(), newbbs.begin(), newbbs.end());
      obs.insert (obs.end(), newobs.begin(), newobs.end());
      ps.insert (ps.end(), newps.begin(), newps.end());
    }
    
    // check postconditions
    assert ((int)bbs.size() == nprocs);
    assert ((int)obs.size() == nprocs);
    assert ((int)ps.size() == nprocs);
    for (int n=0; n<(int)ps.size(); ++n) {
      assert ((int)ps.at(n) == p+n);
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
        assert (maxc>=0 and maxc<nslices);
        ++ mynprocs.at(maxc);
        if (DEBUG) cout << "SRA maxc " << maxc << endl;
        if (DEBUG) cout << "SRA mynprocs[maxc] " << mynprocs.at(maxc) << endl;
        -- ncomps_left;
      }
      assert (ncomps_left == 0);
      int sum_nprocs = 0;
      for (int c=0; c<nslices; ++c) {
        sum_nprocs += mynprocs.at(c);
      }
      assert (sum_nprocs == nprocs * ncomps);
    }
    if (DEBUG) cout << "SRA mynprocs " << mynprocs << endl;
    
    vector<ibbox> allbbs;
    vector<bbvect> allobs;
    vector<int> allps;
    
    if (DEBUG) cout << "SRA: splitting regions" << endl;
    for (int c=0, p=0; c<nslices; p+=mynprocs.at(c), ++c) {
      
      const ibbox bb = bbs.at(c);
      const bbvect ob = obs.at(c);
      if (DEBUG) cout << "SRA c " << c << endl;
      if (DEBUG) cout << "SRA p " << p << endl;
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
        const CCTK_REAL rfact = pow(nprocs / prod(rshape), (CCTK_REAL)1/dim);
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
    for (int n=0; n<(int)ps.size(); ++n) {
      ps.at(n) /= ncomps;
      assert (ps.at(n) >= 0 and ps.at(n) < nprocs);
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
          assert (all (! (ipos==0) or clb==rlb));
          assert (all (! (ipos==nprocs_dir-1) or cub==rub));
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
    
    for (int n=0; n<(int)ps.size(); ++n) {
      assert (ps.at(n) == n);
    }
  }
  
  
  
  static void MakeMultigridBoxes (const cGH* cgh,
                                  ibbox const & base,
                                  ibbox const & bb,
                                  bbvect const & ob,
                                  vector<ibbox>& bbs)
  {
    bbs.resize (mglevels);
    bbs.at(0) = bb;
    vector<ibbox> bases(mglevels);
    assert (mglevels >= 1);
    bases.at(0) = base;
    if (mglevels > 1) {
      // boundary offsets
      jjvect nboundaryzones, is_internal, is_staggered, shiftout;
      const int ierr = GetBoundarySpecification
        (2*dim, &nboundaryzones[0][0], &is_internal[0][0],
         &is_staggered[0][0], &shiftout[0][0]);
      assert (!ierr);
      // (distance in grid points between the exterior and the physical boundary)
      iivect offset;
      for (int d=0; d<dim; ++d) {
        for (int f=0; f<2; ++f) {
          assert (! is_staggered[d][f]);
          offset[d][f] = (+ (is_internal[d][f] ? 0 : nboundaryzones[d][f] - 1)
                          + shiftout[d][f]);
        }
      }
      vector<ibbox> bases(mglevels);
      bases.at(0) = base;
      for (int ml=1; ml<mglevels; ++ml) {
        // next finer base
        ivect const fbaselo = bases.at(ml-1).lower();
        ivect const fbasehi = bases.at(ml-1).upper();
        ivect const fbasestr = bases.at(ml-1).stride();
        // this base
        ivect const basestr = fbasestr * mgfact;
        ivect const baselo = fbaselo + (xpose(offset)[0] - ivect(mgfact) * xpose(offset)[0]) * fbasestr;
        ivect const basehi = fbasehi + (xpose(offset)[1] - ivect(mgfact) * xpose(offset)[1]) * fbasestr;
        ivect const baselo1 = baselo;
        ivect const basehi1 = baselo1 + (basehi - baselo1) / basestr * basestr;
        bases.at(ml) = ibbox(baselo1, basehi1, basestr);
        // next finer grid
        ivect const flo = bbs.at(ml-1).lower();
        ivect const fhi = bbs.at(ml-1).upper();
        ivect const fstr = bbs.at(ml-1).stride();
        // this grid
        ivect const str = fstr * mgfact;
        ivect const lo = flo + xpose(ob)[0].ifthen (  (xpose(offset)[0] - ivect(mgfact) * xpose(offset)[0]) * fstr, ivect(0));
        ivect const hi = fhi + xpose(ob)[1].ifthen (- (xpose(offset)[1] - ivect(mgfact) * xpose(offset)[1]) * fstr, ivect(0));
        ivect const lo1 = baselo1 + (lo - baselo1 + str - 1) / str * str;
        ivect const hi1 = lo1 + (hi - lo1) / str * str;
        bbs.at(ml) = ibbox(lo1, hi1, str);
      }
    }
  }
  
  void MakeMultigridBoxes (const cGH* cgh,
                           vector<vector<ibbox> > const & bbss,
                           vector<vector<bbvect> > const & obss,
                           vector<vector<vector<ibbox> > > & bbsss)
  {
    assert (bbss.size() == obss.size());
    for (int rl=0; rl<(int)bbss.size(); ++rl) {
      assert (bbss.at(rl).size() == obss.at(rl).size());
    }
    
    bbsss.resize (mglevels);
    for (int ml=0; ml<mglevels; ++ml) {
      bbsss.at(ml).resize (bbss.size());
      for (int rl=0; rl<(int)bbss.size(); ++rl) {
        bbsss.at(ml).at(rl).resize (bbss.at(rl).size());
      }
    }
    
    for (int rl=0; rl<(int)bbss.size(); ++rl) {
      
      ibbox base;
      for (int c=0; c<(int)bbss.at(rl).size(); ++c) {
        base = base.expanded_containing(bbss.at(rl).at(c));
      }
      
      for (int c=0; c<(int)bbss.at(rl).size(); ++c) {
        vector<ibbox> mg_bbs;
        MakeMultigridBoxes
          (cgh, base, bbss.at(rl).at(c), obss.at(rl).at(c), mg_bbs);
        assert ((int)mg_bbs.size() == mglevels);
        for (int ml=0; ml<mglevels; ++ml) {
          bbsss.at(ml).at(rl).at(c) = mg_bbs.at(ml);
        }
        
      } // for c
      
    } // for rl
  }
  
} // namespace Carpet
