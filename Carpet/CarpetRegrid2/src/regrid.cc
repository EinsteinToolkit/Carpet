#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "bbox.hh"
#include "bboxset.hh"
#include "defs.hh"
#include "dh.hh"
#include "gh.hh"
#include "vect.hh"

#include "carpet.hh"



namespace CarpetRegrid2 {
  
  using namespace std;
  using namespace Carpet;
  
  
  
  typedef bboxset <int, dim> ibboxset;
  typedef vect <CCTK_REAL, dim> rvect;
  typedef bbox <CCTK_REAL, dim> rbbox;
  
  
  
  struct centre_description {
    int                num_levels;
    rvect              position;
    vector <CCTK_REAL> radius;
    
    static int num_centres ();
    centre_description (int n);
  };
  
  
  
  int
  centre_description::
  num_centres ()
  {
    DECLARE_CCTK_PARAMETERS;
    return num_centres;
  }
  
  
  
  centre_description::
  centre_description (int const n)
  {
    DECLARE_CCTK_PARAMETERS;
    
    assert (n >= 0 and n < centre_description::num_centres ());
    switch (n) {
      
    case 0:
      num_levels = num_levels_1;
      position = rvect (position_x_1, position_y_1, position_z_1);
      radius.resize (num_levels);
      for (int rl = 0; rl < num_levels; ++ rl) {
        radius.at(rl) = radius_1[rl];
      }
      break;
      
    case 1:
      num_levels = num_levels_2;
      position = rvect (position_x_2, position_y_2, position_z_2);
      radius.resize (num_levels);
      for (int rl = 0; rl < num_levels; ++ rl) {
        radius.at(rl) = radius_2[rl];
      }
      break;
      
    case 2:
      num_levels = num_levels_3;
      position = rvect (position_x_3, position_y_3, position_z_3);
      radius.resize (num_levels);
      for (int rl = 0; rl < num_levels; ++ rl) {
        radius.at(rl) = radius_3[rl];
      }
      break;
      
    default:
      assert (0);
    }
    
    assert (num_levels <= maxreflevels);
  }
  
  
  
  // Convert a coordinate location to an index location
  ivect
  rpos2ipos (rvect const & rpos,
             rvect const & origin, rvect const & scale,
             gh const & hh, int const rl)
  {
    ivect const levfac = hh.reffacts.at(rl);
    assert (all (hh.baseextent.stride() % levfac == 0));
    ivect const istride = hh.baseextent.stride() / levfac;
    
    return (ivect (floor ((rpos - origin) * scale / rvect(istride) +
                          static_cast<CCTK_REAL> (0.5))) *
            istride);
  }
  
  // Convert a coordinate location to an index location, rounding in
  // the opposite manner as rpos2ipos
  ivect
  rpos2ipos1 (rvect const & rpos,
              rvect const & origin, rvect const & scale,
              gh const & hh, int const rl)
  {
    ivect const levfac = hh.reffacts.at(rl);
    assert (all (hh.baseextent.stride() % levfac == 0));
    ivect const istride = hh.baseextent.stride() / levfac;
    
    return (ivect (ceil ((rpos - origin) * scale / rvect(istride) -
                         static_cast<CCTK_REAL> (0.5))) *
            istride);
  }
  
  // Convert an index location to a coordinate location
  rvect
  ipos2rpos (ivect const & ipos,
             rvect const & origin, rvect const & scale,
             gh const & hh, int const rl)
  {
    return rvect(ipos) / scale + origin;
  }
  
  // Convert an index bbox to a coordinate bbox
  rbbox
  ibbox2rbbox (ibbox const & ib,
               rvect const & origin, rvect const & scale,
               gh const & hh, int const rl)
  {
    rvect const zero (0);
    return rbbox (ipos2rpos (ib.lower() , origin, scale, hh, rl),
                  ipos2rpos (ib.upper() , origin, scale, hh, rl),
                  ipos2rpos (ib.stride(), zero  , scale, hh, rl));
  }
  
  
  
  extern "C" {
    CCTK_INT
    CarpetRegrid2_Regrid (CCTK_POINTER_TO_CONST const cctkGH_,
                          CCTK_POINTER          const bbsss_,
                          CCTK_POINTER          const obss_,
                          CCTK_POINTER          const pss_,
                          CCTK_INT              const force);
    
    CCTK_INT
    CarpetRegrid2_RegridMaps (CCTK_POINTER_TO_CONST const cctkGH_,
                              CCTK_POINTER          const bbssss_,
                              CCTK_POINTER          const obsss_,
                              CCTK_POINTER          const psss_,
                              CCTK_INT              const force);
  }
  
  
  
  void
  Regrid (cGH const * const cctkGH,
          gh::rexts & bbss,
          gh::rbnds & obss)
  {
    DECLARE_CCTK_PARAMETERS;
    
    assert (Carpet::is_singlemap_mode());
    gh const & hh = * Carpet::vhh.at (Carpet::map);
    dh const & dd = * Carpet::vdd.at (Carpet::map);
    
    //
    // Find extent of domain
    //
    
    // This requires that CoordBase is used
#warning "TODO: check this (for Carpet, and maybe also for CartGrid3D)"
#warning "TODO: (the check that these two are consistent should be in Carpet)"
    
    typedef vect<vect<CCTK_INT,2>,3> jjvect;
    jjvect nboundaryzones;
    jjvect is_internal;
    jjvect is_staggered;
    jjvect shiftout;
    {
      CCTK_INT const ierr = GetBoundarySpecification
        (2*dim,
         & nboundaryzones[0][0],
         & is_internal[0][0],
         & is_staggered[0][0],
         & shiftout[0][0]);
      assert (not ierr);
    }
    
    rvect physical_lower, physical_upper;
    rvect interior_lower, interior_upper;
    rvect exterior_lower, exterior_upper;
    rvect spacing;
    {
      CCTK_INT const ierr = GetDomainSpecification
        (dim,
         & physical_lower[0], & physical_upper[0],
         & interior_lower[0], & interior_upper[0],
         & exterior_lower[0], & exterior_upper[0],
         & spacing[0]);
      assert (not ierr);
    }
    
    // Adapt spacing for convergence level
    spacing *= ipow ((CCTK_REAL) Carpet::mgfact, Carpet::basemglevel);
    
    {
      CCTK_INT const ierr =
        ConvertFromPhysicalBoundary
        (dim,
         & physical_lower[0], & physical_upper[0],
         & interior_lower[0], & interior_upper[0],
         & exterior_lower[0], & exterior_upper[0],
         & spacing[0]);
      assert (not ierr);
    }
    
    //
    // Calculate the union of the bounding boxes for all levels
    //
    
    rvect const origin (exterior_lower);
    rvect const scale (rvect (hh.baseextent.stride()) / spacing);
    
    ivect const physical_ilower =
      rpos2ipos (physical_lower, origin, scale, hh, 0);
    ivect const physical_iupper =
      rpos2ipos1 (physical_upper, origin, scale, hh, 0);
    
    // The set of refined regions
    vector <ibboxset> regions (1);
    
    // Loop over all centres
    for (int n = 0; n < centre_description::num_centres(); ++ n) {
      centre_description centre (n);
      
      // Loop over all levels for this centre
      for (int rl = 1; rl < centre.num_levels; ++ rl) {
        
        // Calculate a bbox for this region
        rvect const rmin = centre.position - centre.radius.at(rl);
        rvect const rmax = centre.position + centre.radius.at(rl);
        
        // Convert to an integer bbox
        ivect const imin =
          rpos2ipos (rmin, origin, scale, hh, rl);
        ivect const imax =
          rpos2ipos1 (rmax, origin, scale, hh, rl);
        
        ivect const levfac = hh.reffacts.at(rl);
        assert (all (hh.baseextent.stride() % levfac == 0));
        ivect const istride = hh.baseextent.stride() / levfac;
        
        ibbox const region (imin, imax, istride);
        
        // Add this region to the list of regions
        if (static_cast <int> (regions.size()) < rl+1) regions.resize (rl+1);
        regions.at(rl) |= region;
        
      } // for rl
      
    } // for n
    
      // Ensure that the coarser grids contain the finer ones
    for (size_t rl = regions.size() - 1; rl >= 2; -- rl) {
      ibbox coarse = * regions.at(rl-1).begin();
      
      regions.at(rl).normalize();
      ibboxset coarsified;
      for (ibboxset::const_iterator ibb = regions.at(rl).begin();
           ibb != regions.at(rl).end();
           ++ ibb)
      {
        ivect const distance (min_distance);
        ibbox const fbb = * ibb;
        ibbox const cbb = fbb.expanded_for(coarse);
        ibbox const ebb = cbb.expand (distance, distance);
        coarsified |= ebb;
      }
      
      regions.at(rl-1) |= coarsified;
    }
    
    
    
    //
    // Clip at the outer boundary, and convert to (bbss, obss) pair
    //
    
    bbss.resize (regions.size());
    obss.resize (regions.size());
    
    for (size_t rl = 1; rl < regions.size(); ++ rl) {
      
      // Find the location of the outer boundary
      rvect const level_physical_lower = physical_lower;
      rvect const level_physical_upper = physical_upper;
      rvect const level_spacing = spacing / rvect (hh.reffacts.at(rl));
      rvect level_interior_lower, level_interior_upper;
      rvect level_exterior_lower, level_exterior_upper;
      {
        CCTK_INT const ierr =
          ConvertFromPhysicalBoundary
          (dim,
           & level_physical_lower[0], & level_physical_upper[0],
           & level_interior_lower[0], & level_interior_upper[0],
           & level_exterior_lower[0], & level_exterior_upper[0],
           & level_spacing[0]);
        assert (not ierr);
      }
      
      ivect const level_exterior_ilower =
        rpos2ipos (level_exterior_lower, origin, scale, hh, rl);
      ivect const level_exterior_iupper =
        rpos2ipos1 (level_exterior_upper, origin, scale, hh, rl);
      
      // Find the minimum necessary distance to the outer boundary
      // due to buffer zones.  This is in terms of grid points.
      i2vect const min_bnd_dist = dd.buffers;
      
      // Clip at the outer boundary
      regions.at(rl).normalize();
      ibboxset clipped;
      for (ibboxset::const_iterator ibb = regions.at(rl).begin();
           ibb != regions.at(rl).end();
           ++ ibb)
      {
        ibbox const bb = * ibb;
        
        // Clip boxes that extend outside the boundary.  Enlarge boxes
        // that are inside but too close to the outer boundary.
        bvect const lower_is_outer =
          bb.lower() - min_bnd_dist[0] <= physical_ilower;
        bvect const upper_is_outer =
          bb.upper() + min_bnd_dist[1] >= physical_iupper;
        
        ibbox const clipped_bb
          (either (lower_is_outer, level_exterior_ilower, bb.lower()),
           either (upper_is_outer, level_exterior_iupper, bb.upper()),
           bb.stride());
        
        clipped += clipped_bb;
      }
      
      regions.at(rl) = clipped;
      
      // Simplify the set of refined regions
      regions.at(rl).normalize();
      
      // Create a vector of bboxes for this level
      gh::cexts bbs;
      gh::cbnds obs;
      bbs.reserve (regions.at(rl).setsize());
      obs.reserve (regions.at(rl).setsize());
      for (ibboxset::const_iterator ibb = regions.at(rl).begin();
           ibb != regions.at(rl).end();
           ++ ibb)
      {
        ibbox bb = * ibb;
        
        bvect const lower_is_outer = bb.lower() <= physical_ilower;
        bvect const upper_is_outer = bb.upper() >= physical_iupper;
        
        bbvect ob;
        for (int d = 0; d < dim; ++ d) {
          ob[d][0] = lower_is_outer[d];
          ob[d][1] = upper_is_outer[d];
        }
        bbs.push_back (bb);
        obs.push_back (ob);
      }
      
      bbss.at(rl) = bbs;
      obss.at(rl) = obs;
      
    } // for rl
  }
  
  
  
  CCTK_INT
  CarpetRegrid2_Regrid (CCTK_POINTER_TO_CONST const cctkGH_,
                        CCTK_POINTER          const bbsss_,
                        CCTK_POINTER          const obss_,
                        CCTK_POINTER          const pss_,
                        CCTK_INT              const force)
  {
    DECLARE_CCTK_PARAMETERS;
    
    cGH const * const cctkGH = static_cast <cGH const *> (cctkGH_);
    
    assert (Carpet::is_singlemap_mode());
    
    // Decide whether to change the grid hierarchy
    // (We always do)
    bool do_recompose;
    if (force) {
      do_recompose = true;
    } else {
      if (regrid_every == -1) {
        do_recompose = false;
      } else if (regrid_every == 0) {
        do_recompose = cctkGH->cctk_iteration == 0;
      } else {
        do_recompose =
          reflevel == 0 and
          (cctkGH->cctk_iteration == 0 or
           (cctkGH->cctk_iteration > 0 and
            (cctkGH->cctk_iteration - 1) % regrid_every == 0));
      }
    }
    
    if (do_recompose) {
      
      gh::mexts  & bbsss = * static_cast <gh::mexts  *> (bbsss_);
      gh::rbnds  & obss  = * static_cast <gh::rbnds  *> (obss_ );
      gh::rprocs & pss   = * static_cast <gh::rprocs *> (pss_  );
      
      // Make multigrid unaware
      vector <vector <ibbox> > bbss = bbsss.at(0);
      
      Regrid (cctkGH, bbss, obss);
      
      // Make multiprocessor aware
      pss.resize (bbss.size());
      for (size_t rl = 0; rl < pss.size(); ++ rl) {
        Carpet::SplitRegions (cctkGH, bbss.at(rl), obss.at(rl), pss.at(rl));
      } // for rl
      
      // Make multigrid aware
      Carpet::MakeMultigridBoxes (cctkGH, bbss, obss, bbsss);
      
    } // if do_recompose
    
    return do_recompose;
  }
  
  
  
  CCTK_INT
  CarpetRegrid2_RegridMaps (CCTK_POINTER_TO_CONST const cctkGH_,
                            CCTK_POINTER          const bbssss_,
                            CCTK_POINTER          const obsss_,
                            CCTK_POINTER          const psss_,
                            CCTK_INT              const force)
  {
    DECLARE_CCTK_PARAMETERS;
    
    cGH const * const cctkGH = static_cast <cGH const *> (cctkGH_);
    
    assert (Carpet::is_level_mode());
    
    // Decide whether to change the grid hierarchy
    // (We always do)
    bool do_recompose;
    if (force) {
      do_recompose = true;
    } else {
      if (regrid_every == -1) {
        do_recompose = false;
      } else if (regrid_every == 0) {
        do_recompose = cctkGH->cctk_iteration == 0;
      } else {
        do_recompose =
          reflevel == 0 and
          (cctkGH->cctk_iteration == 0 or
           (cctkGH->cctk_iteration > 0 and
            (cctkGH->cctk_iteration - 1) % regrid_every == 0));
      }
    }
    
    if (do_recompose) {
      vector <gh::mexts>  & bbssss =
        * static_cast <vector <gh::mexts>  *> (bbssss_);
      vector <gh::rbnds>  & obsss  =
        * static_cast <vector <gh::rbnds>  *> (obsss_ );
      vector <gh::rprocs> & psss   =
        * static_cast <vector <gh::rprocs> *> (psss_  );
      
      // Make multigrid unaware
      vector< vector <vector <ibbox> > > bbsss (maps);
      for (int m = 0; m < maps; ++ m) {
        bbsss.at(m) = bbssss.at(m).at(0);
      }
      
      for (int m = 0; m < maps; ++ m) {
        Regrid (cctkGH, bbsss.at(m), obsss.at(m));
      }
      
      // Count levels
      vector <int> rls (maps);
      for (int m = 0; m < maps; ++ m) {
        rls.at(m) = bbsss.at(m).size();
      }
      int maxrl = 0;
      for (int m = 0; m < maps; ++ m) {
        maxrl = max (maxrl, rls.at(m));
      }
      
      // Make multiprocessor aware
      for (int m = 0; m < maps; ++ m) {
        psss.at(m).resize (bbsss.at(m).size());
      }
      for (int rl = 0; rl < maxrl; ++ rl) {
        vector <vector <ibbox > > bbss (maps);
        vector <vector <bbvect> > obss (maps);
        vector <vector <int   > > pss  (maps);
        for (int m = 0; m < maps; ++ m) {
          if (rl < rls.at(m)) {
            bbss.at(m) = bbsss.at(m).at(rl);
            obss.at(m) = obsss.at(m).at(rl);
            pss .at(m) = psss .at(m).at(rl);
          }
        }
        Carpet::SplitRegionsMaps (cctkGH, bbss, obss, pss);
        for (int m = 0; m < maps; ++ m) {
          if (rl < rls.at(m)) {
            bbsss.at(m).at(rl) = bbss.at(m);
            obsss.at(m).at(rl) = obss.at(m);
            psss .at(m).at(rl) = pss .at(m);
          }
        }
      } // for rl
      
      // Make multigrid aware
      Carpet::MakeMultigridBoxesMaps (cctkGH, bbsss, obsss, bbssss);
      
    } // if do_recompose
    
    return do_recompose;
  }
  
  
  
} // namespace CarpetRegrid2
