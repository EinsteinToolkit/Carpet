#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <bbox.hh>
#include <bboxset.hh>
#include <defs.hh>
#include <dh.hh>
#include <gh.hh>
#include <vect.hh>

#include <carpet.hh>

#include "indexing.hh"



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
    
    centre_description (cGH const * cctkGH, int n);
  };
  
  
  
  centre_description::
  centre_description (cGH const * const cctkGH, int const n)
  {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;
    
    assert (n >= 0 and n < num_centres);
    
    int lsh[2];
    getvectorindex2 (cctkGH, "CarpetRegrid2::radii", lsh);
    
    this->num_levels = num_levels[n];
    this->position = rvect (position_x[n], position_y[n], position_z[n]);
    this->radius.resize (this->num_levels);
    for (int rl = 0; rl < this->num_levels; ++ rl) {
      this->radius.at(rl) = radius[index2 (lsh, rl, n)];
    }
    
    assert (this->num_levels <= maxreflevels);
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
    
    if (verbose) CCTK_INFO ("Regridding");
    
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
    for (int n = 0; n < num_centres; ++ n) {
      centre_description centre (cctkGH, n);
      
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
      ibbox const coarse = * regions.at(rl-1).begin();
      
      i2vect const fdistance =
        i2vect(ivect(0 * min_distance)) + dd.ghosts + dd.buffers;
      i2vect const cdistance =
        i2vect(ivect(min_distance + dd.inner_buffer_width +
                     dd.prolongation_stencil_size()));
      
      regions.at(rl).normalize();
      ibboxset coarsified;
      for (ibboxset::const_iterator ibb = regions.at(rl).begin();
           ibb != regions.at(rl).end();
           ++ ibb)
      {
        ibbox const fbb = * ibb;
        ibbox const efbb = fbb.expand (fdistance[0], fdistance[1]);
        ibbox const cbb = efbb.expanded_for(coarse);
        ibbox const ecbb = cbb.expand (cdistance[0], cdistance[1]);
        coarsified |= ecbb;
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
      
      // Find the minimum necessary distance away from the outer
      // boundary due to buffer and ghost zones.  This is e.g. the
      // distance that the lower boundary of a bbox has to have from
      // the lower boundary.  This is in terms of grid points.
      i2vect const min_bnd_dist_away = dd.buffers + dd.ghosts;
      // Find the minimum necessary distance from the outer boundary
      // due to buffer and ghost zones.  This is e.g. the distance
      // that the upper boundary of a bbox has to have from the lower
      // boundary.  This is in terms of grid points.
      i2vect const min_bnd_dist_incl = dd.ghosts;
      
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
        bvect const lower_is_outside_lower =
          bb.lower() - min_bnd_dist_away[0] * bb.stride() <= physical_ilower;
        // Remove bboxes that are completely outside.
        bvect const upper_is_outside_lower =
          bb.upper() < physical_ilower;
        // Enlarge bboxes that extend not far enough inwards.
        bvect const upper_is_almost_outside_lower =
          bb.upper() < physical_ilower + min_bnd_dist_incl[0] * bb.stride();
        
        // Ditto for the upper boundary.
        bvect const upper_is_outside_upper =
          bb.upper() + min_bnd_dist_away[1] * bb.stride() >= physical_iupper;
        bvect const lower_is_outside_upper =
          bb.lower() > physical_iupper;
        bvect const lower_is_almost_outside_upper =
          bb.lower() > physical_iupper - min_bnd_dist_incl[1] * bb.stride();
        
        assert (not any (lower_is_almost_outside_upper and
                         lower_is_outside_lower));
        assert (not any (upper_is_almost_outside_lower and
                         upper_is_outside_upper));
        
        if (any (upper_is_outside_lower or lower_is_outside_upper)) {
          // The box is completely outside.  Ignore it.
          continue;
        }
        
        ibbox const clipped_bb
          (either (lower_is_outside_lower,
                   level_exterior_ilower,
                   either (lower_is_almost_outside_upper,
                           physical_iupper - min_bnd_dist_incl[1] * bb.stride(),
                           bb.lower())),
           either (upper_is_outside_upper,
                   level_exterior_iupper,
                   either (upper_is_almost_outside_lower,
                           physical_ilower + min_bnd_dist_incl[0] * bb.stride(),
                           bb.upper())),
           bb.stride());
        
        clipped |= clipped_bb;
        
        if (symmetry_rotating90) {
          assert (0);           // There is also a paramwarn about this
        }
        
        if (symmetry_rotating180) {
          // Make the boxes rotating-180 symmetric
          if (lower_is_outside_lower[0]) {
            ivect const ilo = clipped_bb.lower();
            ivect const iup = clipped_bb.upper();
            ivect const istr = clipped_bb.stride();
            
            // Origin
            rvect const axis (physical_lower[0],
                              (physical_lower[1] + physical_upper[1]) / 2,
                              physical_lower[2]); // z component is unused
            ivect const iaxis0 = rpos2ipos  (axis, origin, scale, hh, rl);
            assert (all (iaxis0 % istr == 0));
            ivect const iaxis1 = rpos2ipos1 (axis, origin, scale, hh, rl);
            assert (all (iaxis1 % istr == 0));
            ivect const offset = iaxis1 - iaxis0;
            assert (all (offset % istr == 0));
            assert (all (offset >= 0 and offset < 2*istr));
            assert (all ((iaxis0 + iaxis1 - offset) % (2*istr) == 0));
            ivect const iaxis = (iaxis0 + iaxis1 - offset) / 2;
            
            // Mirror about the y axis, cutting off most of the extent
            // in the positive x direction
            ivect const new_ilo (ilo[0],
                                 (2*iaxis[1]+offset[1]) - iup[1],
                                 ilo[2]);
            ivect const new_iup ((2*iaxis[0]+offset[0]) - ilo[0],
                                 (2*iaxis[1]+offset[1]) - ilo[1],
                                 iup[2]);
            ivect const new_istr (istr);
            
            ibbox const new_bb (new_ilo, new_iup, new_istr);
            
            clipped |= new_bb;
          }
        }
      }
      
      regions.at(rl) = clipped;
      
      // Simplify the set of refined regions
      regions.at(rl).normalize();
      
      // Check whether it would be worthwhile to combine the regions
      // into a single region
      {
        ibbox single;
        for (ibboxset::const_iterator ibb = regions.at(rl).begin();
             ibb != regions.at(rl).end();
             ++ ibb)
        {
          ibbox bb = * ibb;
          single = single.expanded_containing (bb);
        }
        
        CCTK_REAL const regions_size
          = static_cast <CCTK_REAL> (regions.at(rl).size());
        CCTK_REAL const single_size
          = static_cast <CCTK_REAL> (single.size());
        
        // Would a single bbox be efficient enough?
        if (regions_size >= min_fraction * single_size) {
          // Combine the boxes
          regions.at(rl) = ibboxset (single);
        }
      }
      
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
    cGH const * const cctkGH = static_cast <cGH const *> (cctkGH_);
    
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;
    
    assert (Carpet::is_singlemap_mode());
    
    // Decide whether to change the grid hierarchy
    bool do_recompose;
    if (force) {
      do_recompose = true;
    } else {
      if (regrid_every == -1) {
        do_recompose = false;
      } else if (regrid_every == 0) {
        do_recompose = cctk_iteration == 0;
      } else {
        do_recompose =
          (cctk_iteration == 0 or
           (cctk_iteration > 0 and
            (cctk_iteration - 1) % regrid_every == 0 and
            (cctk_iteration > * last_iteration or
             (cctk_iteration == * last_iteration and
              Carpet::map > * last_map))));
      }
    }
    if (do_recompose) {
      * last_iteration = cctk_iteration;
      * last_map = Carpet::map;
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
    cGH const * const cctkGH = static_cast <cGH const *> (cctkGH_);
    
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;
    
    assert (Carpet::is_level_mode());
    
    // Decide whether to change the grid hierarchy
    assert (* last_map == -1);  // ensure this remains unused
    bool do_recompose;
    if (force) {
      do_recompose = true;
    } else {
      if (regrid_every == -1) {
        do_recompose = false;
      } else if (regrid_every == 0) {
        do_recompose = cctk_iteration == 0;
      } else {
        do_recompose =
          (cctk_iteration == 0 or
           (cctk_iteration > 0 and
            (cctk_iteration - 1) % regrid_every == 0 and
            cctk_iteration > * last_iteration));
      }
    }
    if (do_recompose) {
      * last_iteration = cctk_iteration;
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
      
      BEGIN_MAP_LOOP (cctkGH, CCTK_GF) {
        Regrid (cctkGH, bbsss.at(Carpet::map), obsss.at(Carpet::map));
      } END_MAP_LOOP;
      
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
