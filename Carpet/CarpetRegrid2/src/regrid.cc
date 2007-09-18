#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <vector>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <bbox.hh>
#include <bboxset.hh>
#include <defs.hh>
#include <dh.hh>
#include <gh.hh>
#include <region.hh>
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
    int                active;
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
    this->active = active[n];
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
    DECLARE_CCTK_PARAMETERS;
    
    ivect const istride = hh.baseextents.at(0).at(rl).stride();
    
    ivect ipos = (ivect (floor ((rpos - origin) * scale / rvect(istride) +
                                static_cast<CCTK_REAL> (0.5))) *
                  istride);
    
    if (snap_to_coarse) {
      if (rl > 0) {
        ivect const cistride = hh.baseextents.at(0).at(rl - 1).stride();
        ipos = ipos / cistride * cistride;
      }
    }
    
    if (hh.refcent == cell_centered) {
#if 0
      if (rl > 0) {
        ivect const cistride = hh.baseextents.at(0).at(rl - 1).stride();
        
        ivect const offset = (cistride / istride + 1) % 2;
        
        assert (all (istride % offset == 0));
        ipos += offset;
      }
#endif
      assert (all (istride % 2 == 0));
      ipos += istride / 2;
    }
    
    return ipos;
  }
  
  // Convert a coordinate location to an index location, rounding in
  // the opposite manner as rpos2ipos
  ivect
  rpos2ipos1 (rvect const & rpos,
              rvect const & origin, rvect const & scale,
              gh const & hh, int const rl)
  {
    DECLARE_CCTK_PARAMETERS;
    
    ivect const istride = hh.baseextents.at(0).at(rl).stride();
    
    ivect ipos = (ivect (ceil ((rpos - origin) * scale / rvect(istride) -
                               static_cast<CCTK_REAL> (0.5))) *
                  istride);
    
    if (snap_to_coarse) {
      if (rl > 0) {
        ivect const cistride = hh.baseextents.at(0).at(rl - 1).stride();
        ipos = (ipos + cistride - 1) / cistride * cistride;
      }
    }
    
    if (hh.refcent == cell_centered) {
#if 0
      if (rl > 0) {
        ivect const cistride = hh.baseextents.at(0).at(rl - 1).stride();
        
        ivect const offset = (cistride / istride + 1) % 2;
        
        assert (all (istride % offset == 0));
        ipos -= offset;
      }
#endif
      assert (all (istride % 2 == 0));
      ipos -= istride / 2;
    }
    
    return ipos;
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
                          CCTK_POINTER          const regsss_,
                          CCTK_INT              const force);
    
    CCTK_INT
    CarpetRegrid2_RegridMaps (CCTK_POINTER_TO_CONST const cctkGH_,
                              CCTK_POINTER          const regssss_,
                              CCTK_INT              const force);
  }
  
  
  
  void
  Regrid (cGH const * const cctkGH,
          gh::rregs & regss)
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
#warning "TODO: use hh.baseextents[0] [rl] instead of [0]"
    rvect const scale (rvect (hh.baseextents.at(0).at(0).stride()) / spacing);
    
    ivect const physical_ilower =
      rpos2ipos (physical_lower, origin, scale, hh, 0);
    ivect const physical_iupper =
      rpos2ipos1 (physical_upper, origin, scale, hh, 0);
    i2vect const physical_ibounds = i2vect (physical_ilower, physical_iupper);
    
    // The set of refined regions
    vector <ibboxset> regions (1);
    
    // Set up coarsest level
    for (size_t c = 0; c < regss.at(0).size(); ++ c) {
      regions.at(0) += regss.at(0).at(c).extent;
    }
    regions.at(0).normalize();
    
    // Loop over all centres
    for (int n = 0; n < num_centres; ++ n) {
      centre_description centre (cctkGH, n);
      if (centre.active) {
        
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
          
          ivect const istride = hh.baseextents.at(0).at(rl).stride();
          
          ibbox const region (imin, imax, istride);
          
          // Add this region to the list of regions
          if (static_cast <int> (regions.size()) < rl+1) regions.resize (rl+1);
          regions.at(rl) |= region;
          regions.at(rl).normalize();
          
        } // for rl
        
      } // if centre is active
    }   // for n
    
    
      
    regss.resize (regions.size());
    
    assert (regions.size() > 0);
    for (size_t rl = regions.size() - 1; rl >= 1; -- rl) {
      
      // Sanity check
      assert (not regions.at(rl).empty());
      
      
      
      //
      // Ensure that this grid contains the next finer grid
      //
      if (rl < regions.size() - 1) {
        
        ibbox const & coarse0 = * regions.at(rl).begin();
        
        i2vect const fdistance = dd.ghost_width;
        i2vect const cdistance =
          i2vect (min_distance + dd.prolongation_stencil_size());
        
        if (ensure_proper_nesting) {
          
          for (ibboxset::const_iterator ibb = regions.at(rl+1).begin();
               ibb != regions.at(rl+1).end();
               ++ ibb)
          {
            ibbox const & fbb = * ibb;
            
            bvect const lower_is_outer = fbb.lower() <= physical_ilower;
            bvect const upper_is_outer = fbb.upper() >= physical_iupper;
            b2vect const ob (lower_is_outer, upper_is_outer);
            
            ibbox const ebb = fbb.expand (i2vect (not ob) * fdistance);
            ibbox const cbb = ebb.expanded_for (coarse0);
            ibbox const ecbb = cbb.expand (i2vect (not ob) * cdistance);
            
            // Enlarge this level
            regions.at(rl) |= ecbb;
          }
          
          regions.at(rl).normalize();
          
        } // if ensure proper nesting
        
        {
          bool is_properly_nested = true;
          
          for (ibboxset::const_iterator ibb = regions.at(rl+1).begin();
               ibb != regions.at(rl+1).end();
               ++ ibb)
          {
            ibbox const & fbb = * ibb;
            
            bvect const lower_is_outer = fbb.lower() <= physical_ilower;
            bvect const upper_is_outer = fbb.upper() >= physical_iupper;
            b2vect const ob (lower_is_outer, upper_is_outer);
            
            ibbox const ebb = fbb.expand (i2vect (ob) * fdistance);
            ibbox const cbb = ebb.expanded_for (coarse0);
            ibbox const ecbb = cbb.expand (i2vect (ob) * cdistance);
            
            is_properly_nested = is_properly_nested and ecbb <= regions.at(rl);
          }
          
          if (not is_properly_nested) {
            ostringstream msg;
            msg << "Level " << rl << " of the refinement hierarchy is not properly nested.  It does not contain level " << (rl+1) << ".";
            CCTK_WARN (CCTK_WARN_ALERT, msg.str().c_str());
          }
        }
        
      }
      
      
      
      //
      // Add buffer zones
      //
      {
        ibboxset buffered;
        for (ibboxset::const_iterator ibb = regions.at(rl).begin();
             ibb != regions.at(rl).end();
             ++ ibb)
        {
          ibbox const & bb = * ibb;
          
          // bvect const lower_is_outer = bb.lower() <= physical_ilower;
          // bvect const upper_is_outer = bb.upper() >= physical_iupper;
          // 
          // b2vect const ob (lower_is_outer, upper_is_outer);
          // 
          // ibbox const bbb = bb.expand (i2vect (not ob) * dd.buffer_width);
          
          ibbox const bbb = bb.expand (dd.buffer_width);

          buffered |= bbb;
        }
        regions.at(rl) = buffered;
        regions.at(rl).normalize();
      }
      
      
      
      //
      // Check whether it would be worthwhile to combine all regions
      // into a single region
      //
      // TODO: Check this also for pairs of regions
      //
      // TODO: Check after clipping
      {
        ibbox single;
        for (ibboxset::const_iterator ibb = regions.at(rl).begin();
             ibb != regions.at(rl).end();
             ++ ibb)
        {
          ibbox const & bb = * ibb;
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
      i2vect const min_bnd_dist_away = dd.ghost_width;
      // Find the minimum necessary distance from the outer boundary
      // due to buffer and ghost zones.  This is e.g. the distance
      // that the upper boundary of a bbox has to have from the lower
      // boundary.  This is in terms of grid points.
      i2vect const min_bnd_dist_incl = dd.ghost_width;
      // TODO: The above is required only near symmetry boundaries.
      
      
      
      //
      // Clip at the outer boundary
      //
      {
        ibboxset clipped;
        for (ibboxset::const_iterator ibb = regions.at(rl).begin();
             ibb != regions.at(rl).end();
             ++ ibb)
        {
          ibbox const & bb = * ibb;
          
          // Clip boxes that extend outside the boundary.  Enlarge
          // boxes that are inside but too close to the outer
          // boundary.
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
                             (physical_iupper -
                              min_bnd_dist_incl[1] * bb.stride()),
                             bb.lower())),
             either (upper_is_outside_upper,
                     level_exterior_iupper,
                     either (upper_is_almost_outside_lower,
                             (physical_ilower +
                              min_bnd_dist_incl[0] * bb.stride()),
                             bb.upper())),
             bb.stride());
          if (not clipped_bb.is_contained_in (hh.baseextents.at(0).at(rl))) {
            ostringstream msg;
            msg << "Level " << rl << " of the refinement hierarchy is not contained in the simulation domain."
                << "  (There may be too many ghost of buffer zones.)"
                << "  One bbox is " << clipped_bb << "."
                << "  lower_is_outside_lower=" << lower_is_outside_lower
                << "  upper_is_outside_upper=" << upper_is_outside_upper
                << "  lower_is_almost_outside_upper=" << lower_is_almost_outside_upper
                << "  upper_is_almost_outside_lower=" << upper_is_almost_outside_lower
                << "  level_exterior_ilower=" << level_exterior_ilower
                << "  level_exterior_iupper=" << level_exterior_iupper
                << "  physical_ilower=" << physical_ilower
                << "  physical_iupper=" << physical_iupper
                << "  baseextent=" << hh.baseextents.at(0).at(rl);
            CCTK_WARN (CCTK_WARN_ABORT, msg.str().c_str());
          }
          assert (clipped_bb.is_contained_in (hh.baseextents.at(0).at(rl)));
          
          clipped |= clipped_bb;
          
        } // for ibb
        
        regions.at(rl) = clipped;
        regions.at(rl).normalize();
      }
      
      
      
      //
      // Make the boxes rotating-90 symmetric
      //
      if (symmetry_rotating90) {
        assert (0);            // There is also a paramwarn about this
      }
      
      
      
      //
      // Make the boxes rotating-180 symmetric
      //
      if (symmetry_rotating180) {
        
        ibboxset added;
        for (ibboxset::const_iterator ibb = regions.at(rl).begin();
             ibb != regions.at(rl).end();
             ++ ibb)
        {
          ibbox const & bb = * ibb;
          
          bvect const lower_is_outside_lower =
            bb.lower() - min_bnd_dist_away[0] * bb.stride() <= physical_ilower;
          
          if (lower_is_outside_lower[0]) {
            ivect const ilo = bb.lower();
            ivect const iup = bb.upper();
            ivect const istr = bb.stride();
            
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
            assert (new_bb.is_contained_in (hh.baseextents.at(0).at(rl)));
            
            added |= new_bb;
          }
        }
        
        regions.at(rl) |= added;
        regions.at(rl).normalize();
        
      } // if symmetry_rotating180
      
      
      
      //
      // Ensure that this grid is contained in the domain
      //
      assert (regions.at(rl) <= hh.baseextents.at(0).at(rl));
      
      
      
      //
      // Create a vector of bboxes for this level
      //
      {
        gh::cregs regs;
        regs.reserve (regions.at(rl).setsize());
        for (ibboxset::const_iterator ibb = regions.at(rl).begin();
             ibb != regions.at(rl).end();
             ++ ibb)
        {
          ibbox const & bb = * ibb;
          assert (bb.is_contained_in (hh.baseextents.at(0).at(rl)));
          
          bvect const lower_is_outer = bb.lower() <= physical_ilower;
          bvect const upper_is_outer = bb.upper() >= physical_iupper;
          
          b2vect const ob (lower_is_outer, upper_is_outer);
          
          region_t reg;
          reg.extent = bb;
          reg.map = Carpet::map;
          reg.outer_boundaries = ob;
          regs.push_back (reg);
        }
        
        regss.at(rl) = regs;
      }
      
    } // for rl
  }
  
  
  
  CCTK_INT
  CarpetRegrid2_Regrid (CCTK_POINTER_TO_CONST const cctkGH_,
                        CCTK_POINTER          const regsss_,
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
        // Regrid at most once per iteration
        do_recompose =
          (cctk_iteration == 0 or
           (cctk_iteration > 0 and
            (cctk_iteration - 1) % regrid_every == 0 and
            (cctk_iteration > * last_iteration or
             (cctk_iteration == * last_iteration and
              Carpet::map > * last_map))));
      }
    }
    
    if (verbose) {
      if (do_recompose) {
        for (int n = 0; n < num_centres; ++ n) {
          if (active[n]) {
            CCTK_VInfo (CCTK_THORNSTRING,
                        "Centre %d is at position [%g,%g,%g]",
                        n,
                        static_cast <double> (position_x[n]),
                        static_cast <double> (position_y[n]),
                        static_cast <double> (position_z[n]));
          } else {
            CCTK_VInfo (CCTK_THORNSTRING,
                        "Centre %d is not active", n);
          }
        }
      }
    }
    
    if (do_recompose and * last_iteration != -1) {
      // Regrid only if the positions have changed sufficiently
      do_recompose = false;
      for (int n = 0; n < num_centres; ++ n) {
        CCTK_REAL const dist2 =
          pow (position_x[n] - old_position_x[n], 2) +
          pow (position_y[n] - old_position_y[n], 2) +
          pow (position_z[n] - old_position_z[n], 2);
        CCTK_REAL mindist;
        switch (n) {
        case 0: mindist = movement_threshold_1; break;
        case 1: mindist = movement_threshold_2; break;
        case 2: mindist = movement_threshold_3; break;
        case 3: mindist = movement_threshold_4; break;
        case 4: mindist = movement_threshold_5; break;
        case 5: mindist = movement_threshold_6; break;
        case 6: mindist = movement_threshold_7; break;
        case 7: mindist = movement_threshold_8; break;
        case 8: mindist = movement_threshold_9; break;
        case 9: mindist = movement_threshold_10; break;
        default: assert (0);
        }
        do_recompose = dist2 >= pow (mindist, 2);
        if (do_recompose) break;
      } // for n
      if (verbose) {
        if (not do_recompose) {
          CCTK_INFO
            ("Centres have not moved sufficiently; skipping regridding");
        }
      }
    }
    
    if (do_recompose) {
      * last_iteration = cctk_iteration;
      * last_map = Carpet::map;
    }
    
    if (do_recompose) {
      
      gh::mregs & regsss = * static_cast <gh::mregs *> (regsss_);
      
      // Make multigrid unaware
      vector <vector <region_t> > regss = regsss.at(0);
      
      Regrid (cctkGH, regss);
      
      // Make multiprocessor aware
      for (size_t rl = 0; rl < regss.size(); ++ rl) {
        Carpet::SplitRegions (cctkGH, regss.at(rl));
      } // for rl
      
      // Make multigrid aware
      Carpet::MakeMultigridBoxes (cctkGH, regss, regsss);
      
      // Remember current positions
      for (int n = 0; n < num_centres; ++ n) {
        old_position_x[n] = position_x[n];
        old_position_y[n] = position_y[n];
        old_position_z[n] = position_z[n];
      }
      
    } // if do_recompose
    
    return do_recompose;
  }
  
  
  
  CCTK_INT
  CarpetRegrid2_RegridMaps (CCTK_POINTER_TO_CONST const cctkGH_,
                            CCTK_POINTER          const regssss_,
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
        // Regrid at most once per iteration
        do_recompose =
          (cctk_iteration == 0 or
           (cctk_iteration > 0 and
            (cctk_iteration - 1) % regrid_every == 0 and
            cctk_iteration > * last_iteration));
      }
    }
    
    if (verbose) {
      if (do_recompose) {
        for (int n = 0; n < num_centres; ++ n) {
          if (active[n]) {
            CCTK_VInfo (CCTK_THORNSTRING,
                        "Centre %d is at position [%g,%g,%g]",
                        n,
                        static_cast <double> (position_x[n]),
                        static_cast <double> (position_y[n]),
                        static_cast <double> (position_z[n]));
          } else {
            CCTK_VInfo (CCTK_THORNSTRING,
                        "Centre %d is not active", n);
          }
        }
      }
    }
    
    if (do_recompose and * last_iteration != -1) {
      // Regrid only if the positions have changed sufficiently
      do_recompose = false;
      for (int n = 0; n < num_centres; ++ n) {
        CCTK_REAL const dist2 =
          pow (position_x[n] - old_position_x[n], 2) +
          pow (position_y[n] - old_position_y[n], 2) +
          pow (position_z[n] - old_position_z[n], 2);
        CCTK_REAL mindist;
        switch (n) {
        case 0: mindist = movement_threshold_1; break;
        case 1: mindist = movement_threshold_2; break;
        case 2: mindist = movement_threshold_3; break;
        case 3: mindist = movement_threshold_4; break;
        case 4: mindist = movement_threshold_5; break;
        case 5: mindist = movement_threshold_6; break;
        case 6: mindist = movement_threshold_7; break;
        case 7: mindist = movement_threshold_8; break;
        case 8: mindist = movement_threshold_9; break;
        case 9: mindist = movement_threshold_10; break;
        default: assert (0);
        }
        do_recompose = dist2 >= pow (mindist, 2);
        if (do_recompose) break;
      } // for n
      if (verbose) {
        if (not do_recompose) {
          CCTK_INFO
            ("Centres have not moved sufficiently; skipping regridding");
        }
      }
    }
    
    if (do_recompose) {
      * last_iteration = cctk_iteration;
    }
    
    if (do_recompose) {
      
      vector <gh::mregs> & regssss =
        * static_cast <vector <gh::mregs>  *> (regssss_);
      
      // Make multigrid unaware
      vector< vector <vector <region_t> > > regsss (Carpet::maps);
      for (int m = 0; m < maps; ++ m) {
        regsss.at(m) = regssss.at(m).at(0);
      }
      
      BEGIN_MAP_LOOP (cctkGH, CCTK_GF) {
        Regrid (cctkGH, regsss.at(Carpet::map));
      } END_MAP_LOOP;
      
      // Count levels
      vector <int> rls (maps);
      for (int m = 0; m < maps; ++ m) {
        rls.at(m) = regsss.at(m).size();
      }
      int maxrl = 0;
      for (int m = 0; m < maps; ++ m) {
        maxrl = max (maxrl, rls.at(m));
      }
      
      // Make multiprocessor aware
      for (int rl = 0; rl < maxrl; ++ rl) {
        vector <vector <region_t> > regss (maps);
        for (int m = 0; m < maps; ++ m) {
          if (rl < rls.at(m)) {
            regss.at(m) = regsss.at(m).at(rl);
          }
        }
        Carpet::SplitRegionsMaps (cctkGH, regss);
        for (int m = 0; m < maps; ++ m) {
          if (rl < rls.at(m)) {
            regsss.at(m).at(rl) = regss.at(m);
          }
        }
      } // for rl
      
      // Make multigrid aware
      Carpet::MakeMultigridBoxesMaps (cctkGH, regsss, regssss);
      
      // Remember current positions
      for (int n = 0; n < num_centres; ++ n) {
        old_position_x[n] = position_x[n];
        old_position_y[n] = position_y[n];
        old_position_z[n] = position_z[n];
      }
      
    } // if do_recompose
    
    return do_recompose;
  }
  
  
  
} // namespace CarpetRegrid2
