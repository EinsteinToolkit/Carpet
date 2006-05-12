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
#include "gh.hh"
#include "vect.hh"

#include "carpet.hh"



namespace CarpetRegrid2 {
  
  using namespace std;
  
  
  
  typedef bboxset <int, dim> ibboxset;
  typedef vect <CCTK_REAL, dim> rvect;
  typedef bbox <CCTK_REAL, dim> rbbox;
  
  
  
  struct centre_description {
    size_t             num_levels;
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
      for (size_t rl = 0; rl < num_levels; ++ rl) {
        radius.at(rl) = radius_1 [rl];
      }
      break;
      
    case 1:
      num_levels = num_levels_2;
      position = rvect (position_x_2, position_y_2, position_z_2);
      radius.resize (num_levels);
      for (size_t rl = 0; rl < num_levels; ++ rl) {
        radius.at(rl) = radius_2 [rl];
      }
      break;
      
    case 2:
      num_levels = num_levels_3;
      position = rvect (position_x_3, position_y_3, position_z_3);
      radius.resize (num_levels);
      for (size_t rl = 0; rl < num_levels; ++ rl) {
        radius.at(rl) = radius_3 [rl];
      }
      break;
      
    default:
      assert (0);
    }
  }
  
  
  
  ivect
  rpos2ipos (rvect const & rpos,
             rvect const & origin, rvect const & scale,
             gh const & hh, int const rl)
  {
    ivect const levfac = hh.reffacts.at(rl);
    assert (all (hh.baseextent.stride() % levfac == 0));
    ivect const istride = hh.baseextent.stride() / levfac;
    
    return (ivect (floor ((rpos - origin) * scale / rvect(istride) + 0.5)) *
            istride);
  }
  
  rvect
  ipos2rpos (ivect const & ipos,
             rvect const & origin, rvect const & scale,
             gh const & hh, int const rl)
  {
    return rvect(ipos) / scale + origin;
  }
  
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
    gh::mexts  & bbsss = * static_cast <gh::mexts  *> (bbsss_);
    gh::rbnds  & obss  = * static_cast <gh::rbnds  *> (obss_ );
    gh::rprocs & pss   = * static_cast <gh::rprocs *> (pss_  );
    
    assert (Carpet::is_singlemap_mode());
    gh const & hh = * Carpet::vhh.at (Carpet::map);
    
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
          cctkGH->cctk_iteration == 0 or
          (cctkGH->cctk_iteration > 0 and
           (cctkGH->cctk_iteration - 1) % regrid_every == 0);
      }
    }
    
    if (do_recompose) {
      
      // Find extent of domain
      
      // This requires that CoordBase is used (but this is not
      // checked)
      
      iivect nboundaryzones;
      iivect is_internal;
      iivect is_staggered;
      iivect shiftout;
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
      
      rvect const origin (exterior_lower);
      rvect const scale (rvect (hh.baseextent.stride()) / spacing);
      
      ivect const physical_ilower =
        rpos2ipos (physical_lower, origin, scale, hh, 0);
      ivect const physical_iupper =
        rpos2ipos (physical_upper, origin, scale, hh, 0);
      
      // The set of refined regions
      vector <ibboxset> regions (1);
      
      // Loop over all centres
      for (int n = 0; n < centre_description::num_centres(); ++ n) {
        centre_description centre (n);
        
        // Loop over all levels for this centre
        for (size_t rl = 1; rl < centre.num_levels; ++ rl) {
          
          // Calculate a bbox for this region
          rvect const rmin = centre.position - centre.radius.at(rl);
          rvect const rmax = centre.position + centre.radius.at(rl);
          
          // Convert to an integer bbox
          ivect const imin =
            rpos2ipos (rmin, origin, scale, hh, rl);
          ivect const imax =
            rpos2ipos (rmax, origin, scale, hh, rl);
          
          ivect const levfac = hh.reffacts.at(rl);
          assert (all (hh.baseextent.stride() % levfac == 0));
          ivect const istride = hh.baseextent.stride() / levfac;
          
          ibbox const region (imin, imax, istride);
          
          // Add this region to the list of regions
          if (regions.size() < rl+1) regions.resize (rl+1);
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
      
      
      
      // Convert to (bbsss, obss, pss) triplet
      vector <vector <ibbox> > bbss;
      
      bbss.resize (regions.size());
      obss.resize (regions.size());
      pss .resize (regions.size());
      
      bbss.at(0) = bbsss.at(0).at(0);
      
      for (size_t rl = 1; rl < regions.size(); ++ rl) {
        
        rvect const level_physical_lower = physical_lower;
        rvect const level_physical_upper = physical_upper;
        rvect const level_spacing = spacing / hh.reffacts.at(rl);
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
          rpos2ipos (level_exterior_upper, origin, scale, hh, rl);
        
        // Clip at the outer boundary
        regions.at(rl).normalize();
        ibboxset clipped;
        for (ibboxset::const_iterator ibb = regions.at(rl).begin();
             ibb != regions.at(rl).end();
             ++ ibb)
        {
          ibbox const bb = * ibb;
          
          bvect const lower_is_outer = bb.lower() <= physical_ilower;
          bvect const upper_is_outer = bb.upper() >= physical_iupper;
          
          ibbox const clipped_bb (max (bb.lower(), level_exterior_ilower),
                                  min (bb.upper(), level_exterior_iupper),
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
        
        // Make multiprocessor aware
        gh::cprocs ps;
        Carpet::SplitRegions (cctkGH, bbs, obs, ps);
        
        bbss.at(rl) = bbs;
        obss.at(rl) = obs;
        pss .at(rl) = ps;
        
      } // for rl
      
      // Make multigrid aware
      Carpet::MakeMultigridBoxes (cctkGH, bbss, obss, bbsss);
      
    } // if do_recompose
    
    return do_recompose;
  }
  
  
  
} // namespace CarpetRegrid2
