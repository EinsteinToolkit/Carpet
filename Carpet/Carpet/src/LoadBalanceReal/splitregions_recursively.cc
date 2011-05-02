#include <cctk.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <list>
#include <sstream>
#include <string>
#include <vector>

#include <sys/stat.h>
#include <sys/types.h>

#ifdef CCTK_MPI
#  include <mpi.h>
#else
#  include "nompi.h"
#endif

#include <loopcontrol.h>

#include <bbox.hh>
#include <bboxset.hh>
#include <defs.hh>
#include <dh.hh>
#include <gh.hh>
#include <region.hh>
#include <vect.hh>

#include <carpet.hh>
#include <modes.hh>
#include <variables.hh>
#include <Timers.hh>



namespace Carpet {
  
  using namespace std;
  
  
  
  // Helper routines for spliting regions automatically
  
  // The cost for a region, assuming a cost of 1 per interior point
  static rvect
  cost (region_t const & reg)
  {
    DECLARE_CCTK_PARAMETERS;
    static rvect costfactor;
    static bool initialised = false;
    if (not initialised) {
      costfactor = rvect(1.0);
      if (dim > 0) costfactor[0] = 1.0 / aspect_ratio_x;
      if (dim > 1) costfactor[1] = 1.0 / aspect_ratio_y;
      if (dim > 2) costfactor[2] = 1.0 / aspect_ratio_z;
    }
    if (reg.extent.empty()) return rvect(0);
    return rvect (reg.extent.shape() / reg.extent.stride()) * costfactor;
  }
  
  
  
  struct f_range {
    int lower, upper, stride;
  };
  
  struct f_bbox {
    f_range dim[3];
    
    f_bbox () { }
    f_bbox (ibbox const& box)
    {
      assert (::dim == 3);
      for (int d=0; d<3; ++d) {
        dim[d].lower  = box.lower()[d];
        dim[d].upper  = box.upper()[d];
        dim[d].stride = box.stride()[d];
      }
    }
    /*explicit*/ operator ibbox () const
    {
      ivect lower, upper, stride;
      assert (::dim == 3);
      for (int d=0; d<3; ++d) {
        lower [d] = dim[d].lower;
        upper [d] = dim[d].upper;
        stride[d] = dim[d].stride;
      }
      return ibbox (lower, upper, stride);
    }
  };
  
  struct f_boundary {
    int obound[2][3];
    
    /*explicit*/ operator b2vect () const
    {
      b2vect ob;
      assert (::dim == 3);
      for (int d=0; d<3; ++d) {
        for (int f=0; f<2; ++f) {
          ob[f][d] = obound[f][d];
        }
      }
      return ob;
    }
  };
  
  struct f_superregion2slim {
    f_bbox     extent;
    f_boundary outer_boundaries;
    int        map;
    int        processor;
    
    /*explicit*/ operator region_t () const
    {
      region_t reg;
      reg.extent           = ibbox(extent);
      reg.outer_boundaries = b2vect(outer_boundaries);
      reg.map              = map;
      reg.processor        = processor;
      reg.processors       = NULL;
      return reg;
    }
    
    /*explicit*/ operator pseudoregion_t () const
    {
      pseudoregion_t preg;
      preg.extent    = ibbox(extent);
      preg.component = processor;
      return preg;
    }
  };
  
  
  
  extern "C"
  CCTK_FCALL void
  CCTK_FNAME(splitregions_recursively) (CCTK_POINTER const& cxx_superregs,
                                        int const& nsuperregs,
                                        CCTK_POINTER const& cxx_regs,
                                        int const& nprocs);
  
  void
  SplitRegionsMaps_Recursively (cGH const * const cctkGH,
                                vector<vector<region_t> > & superregss,
                                vector<vector<region_t> > & regss)
  {
    DECLARE_CCTK_PARAMETERS;
    
    if (recompose_verbose) cout << "SRMR enter" << endl;
    
    int const nmaps = superregss.size();
    int map_offset = 1000000000;
    for (int m=0; m<nmaps; ++m) {
      for (int r=0; r<int(superregss.AT(m).size()); ++r) {
        map_offset = min (map_offset, superregss.AT(m).AT(r).map);
      }
    }
    
    int nsuperregs = 0;
    for (int m=0; m<nmaps; ++m) {
      nsuperregs += superregss.AT(m).size();
    }
    if (recompose_verbose) cout << "SRMR nsuperregs " << nsuperregs << endl;
    
    // Something to do?
    if (nsuperregs == 0) return;
    
    // Collect slices
    vector<region_t> superregs;
    {
      for (int m=0; m<nmaps; ++m) {
        combine_regions (superregss.AT(m), superregs);
      }
      nsuperregs = superregs.size();
      
      // If the last region was removed, add a new empty region again.
      // A set of regions (corresponding to a refinement level or a
      // grid array) cannot be empty.
      if (nsuperregs == 0) {
        assert (nmaps == 1);    // we should only be here for grid
                                // arrays
        region_t reg;
        reg.extent           = ibbox (ivect (0), ivect (-1), ivect (1));
        reg.outer_boundaries = b2vect (bvect (true), bvect (true));
        reg.map              = 0;
        superregs.push_back (reg);
        nsuperregs = superregs.size();
      }
    }
    
    int const real_nprocs = CCTK_nProcs (cctkGH);
    if (recompose_verbose) cout << "SRMR real_nprocs " << real_nprocs << endl;
    
    // Deactivate some processors if there are too many
    int nprocs;
    if (min_points_per_proc == 0) {
      nprocs = real_nprocs;
    } else {
      CCTK_REAL mycost = 0;
      for (int r=0; r<nsuperregs; ++r) {
        mycost += prod (cost (superregs.AT(r)));
      }
      int const goodnprocs = int (floor (mycost / min_points_per_proc));
      nprocs = max (1, min (real_nprocs, goodnprocs));
    }
    if (recompose_verbose) cout << "SRMR nprocs " << nprocs << endl;
    
    // Distribute load
    CCTK_POINTER const cxx_superregs = & superregs;
    assert ((int)superregs.size() == nsuperregs);
    vector<region_t> regs;
    // regs.reserve (...);
    CCTK_POINTER const cxx_regs = & regs;
    CCTK_FNAME(splitregions_recursively)
      (cxx_superregs, nsuperregs, cxx_regs, nprocs);
    int const nregs = regs.size();
    
    // Allocate regions
    assert ((int)regss.size() == nmaps);
    for (int m=0; m<nmaps; ++m) {
      assert (regss.AT(m).empty());
      // regss.AT(m).reserve (...);
      superregss.AT(m).clear();
      // superregss.AT(m).reserve (...);
    }
    // Assign regions
    for (int r=0; r<nsuperregs; ++r) {
      int const m = superregs.AT(r).map - map_offset;
      assert (m>=0 and m<nmaps);
      superregss.AT(m).push_back (superregs.AT(r));
    }
    for (int r=0; r<nregs; ++r) {
      int const m = regs.AT(r).map - map_offset;
      assert (m>=0 and m<nmaps);
      regss.AT(m).push_back (regs.AT(r));
    }
    // Output regions
    if (recompose_verbose) {
      cout << "SRMR superregss " << superregss << endl;
      cout << "SRMR regss " << regss << endl;
    }
    
    if (recompose_verbose) cout << "SRMR exit" << endl;
  }
  
  
  
  extern "C"
  CCTK_FCALL void
  CCTK_FNAME(carpet_get_region) (CCTK_POINTER& cxx_superregs,
                                 int const& i,
                                 CCTK_POINTER& cxx_superreg)
  {
    vector<region_t>& superregs = *(vector<region_t>*)cxx_superregs;
    region_t& superreg = superregs.AT(i);
    cxx_superreg = &superreg;
  }
  
  extern "C"
  CCTK_FCALL void
  CCTK_FNAME(carpet_get_bbox) (CCTK_POINTER& cxx_superreg,
                               f_bbox& box)
  {
    region_t& superreg = *(region_t*)cxx_superreg;
    box = f_bbox(superreg.extent);
  }
  
  extern "C"
  CCTK_FCALL void
  CCTK_FNAME(carpet_insert_region) (CCTK_POINTER& cxx_regs,
                                    f_superregion2slim const& reg)
  {
    vector<region_t>& regs = *(vector<region_t>*)cxx_regs;
    regs.push_back (region_t (reg));
  }
  
  extern "C"
  CCTK_FCALL void
  CCTK_FNAME(carpet_create_tree_branch) (int const& nch,
                                         int const& dir,
                                         int const fbounds[],
                                         CCTK_POINTER cxx_subtrees[],
                                         CCTK_POINTER& cxx_tree)
  {
    vector<int> bounds(nch+1);
    vector<ipfulltree*> subtrees(nch);
    for (int i=0; i<nch+1; ++i) {
      bounds.AT(i) = fbounds[i];
    }
    for (int i=0; i<nch; ++i) {
      ipfulltree* const tree = (ipfulltree*)cxx_subtrees[i];
      assert (tree->invariant());
      subtrees.AT(i) = tree;
    }
    cxx_tree = new ipfulltree (dir, bounds, subtrees);
  }
  
  extern "C"
  CCTK_FCALL void
  CCTK_FNAME(carpet_create_tree_leaf) (f_superregion2slim const& sreg,
                                       CCTK_POINTER& cxx_tree)
  {
    cxx_tree = new ipfulltree (pseudoregion_t (sreg));
  }
  
  extern "C"
  CCTK_FCALL void
  CCTK_FNAME(carpet_set_tree) (CCTK_POINTER& cxx_superreg,
                               CCTK_POINTER& cxx_tree)
  {
    region_t& superreg = *(region_t*)cxx_superreg;
    ipfulltree* tree = (ipfulltree*)cxx_tree;
    assert (not superreg.processors);
    superreg.processors = tree;
  }
  
} // namespace Carpet
