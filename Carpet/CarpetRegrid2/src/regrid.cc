#include <algorithm>
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
#include <CarpetTimers.hh>

#include "indexing.hh"



namespace CarpetRegrid2 {
  
  using namespace std;
  using namespace Carpet;
  
  
  
  struct centre_description {
    int            _num_levels;
    int            _active;
    rvect          _position;
    vector <rvect> _radius;
    
    centre_description (cGH const * cctkGH, int n);
  };
  
  
  
  centre_description::
  centre_description (cGH const * const cctkGH, int const n)
  {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;
    
    assert (n >= 0 and n < num_centres);
    
    bool found_error = false;
    
    int lsh[2];
    getvectorindex2 (cctkGH, "CarpetRegrid2::radii", lsh);
    
    this->_num_levels = num_levels[n];
    this->_active = active[n];
    this->_position = rvect (position_x[n], position_y[n], position_z[n]);
    if (any (not isfinite(this->_position))) {
      CCTK_VWarn (CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "The position of region %d is [%g,%g,%g], which is not finite",
                  n + 1,
                  double(this->_position[0]),
                  double(this->_position[1]),
                  double(this->_position[2]));
      found_error = true;
    }
    this->_radius.resize (this->_num_levels);
    this->_radius.at(0) = rvect(-1.0, -1.0, -1.0); // unused
    for (int rl = 1; rl < this->_num_levels; ++ rl) {
      int const ind = index2 (lsh, rl, n);
      CCTK_REAL const rx = radius_x[ind] < 0.0 ? radius[ind] : radius_x[ind];
      CCTK_REAL const ry = radius_y[ind] < 0.0 ? radius[ind] : radius_y[ind];
      CCTK_REAL const rz = radius_z[ind] < 0.0 ? radius[ind] : radius_z[ind];
      rvect const rad (rx, ry, rz);
      if (any (not isfinite(rad))) {
        CCTK_VWarn (CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "The radius of refinement level %d of region %d is [%g,%g,%g], which is not finite",
                    rl, n + 1,
                    double(this->_radius.at(rl)[0]),
                    double(this->_radius.at(rl)[1]),
                    double(this->_radius.at(rl)[2]));
        found_error = true;
      }
      if (any (rad < 0.0)) {
        CCTK_VWarn (CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "The radius of refinement level %d of region %d is [%g,%g,%g], which is non-negative",
                    rl, n + 1,
                    double(this->_radius.at(rl)[0]),
                    double(this->_radius.at(rl)[1]),
                    double(this->_radius.at(rl)[2]));
        found_error = true;
      }
      this->_radius.at(rl) = rad;
    }
    
    if (this->_num_levels > maxreflevels) {
      CCTK_VWarn (CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Region %d has %d levels active, which is larger than the maximum number of refinement levels %d",
                  n + 1, this->_num_levels, maxreflevels);
      found_error = true;
    }
    
    if (found_error) {
      CCTK_WARN (CCTK_WARN_ABORT,
                 "Errors found in grid structure specification");
    }
  }
  
  
  
  // Convert a coordinate location to an index location.  For cell
  // centring, shift upwards.
  ivect
  rpos2ipos (rvect const & rpos,
             rvect const & origin, rvect const & scale,
             gh const & hh, int const rl)
  {
    ivect const istride  = hh.baseextents.at(0).at(rl).stride();
    ivect const bistride = hh.baseextents.at(0).at(0).stride();
    
    if (hh.refcent == cell_centered) {
      assert (all (istride % 2 == 0));
    }
    
#if 1
    ivect const ipos =
      hh.refcent == vertex_centered
      ? ivect (floor (((rpos - origin) * scale                    ) / rvect(istride) + rvect(0.5))) * istride
      : ivect (floor (((rpos - origin) * scale - rvect(bistride/2)) / rvect(istride)             )) * istride + istride/2 + bistride/2;
#else
    ivect const ipos =
      hh.refcent == vertex_centered
      ? ivect (floor (((rpos - origin) * scale                    ) / rvect(istride) + rvect(0.5))) * istride
      : ivect (floor (((rpos - origin) * scale - rvect(bistride/2)) / rvect(istride) + rvect(0.5))) * istride + istride/2 + bistride/2;
#endif
    
    return ipos;
  }
  
  // Convert a coordinate location to an index location, rounding in
  // the opposite manner as rpos2ipos.  For cell centring, shift
  // downwards instead of upwards.
  ivect
  rpos2ipos1 (rvect const & rpos,
              rvect const & origin, rvect const & scale,
              gh const & hh, int const rl)
  {
    ivect const istride  = hh.baseextents.at(0).at(rl).stride();
    ivect const bistride = hh.baseextents.at(0).at(0).stride();
    
    if (hh.refcent == cell_centered) {
      assert (all (istride % 2 == 0));
    }
    
#if 1
    ivect const ipos =
      hh.refcent == vertex_centered
      ? ivect (ceil (((rpos - origin) * scale                    ) / rvect(istride) - rvect(0.5))) * istride
      : ivect (ceil (((rpos - origin) * scale - rvect(bistride/2)) / rvect(istride)             )) * istride - istride/2 + bistride/2;
#else
    ivect const ipos =
      hh.refcent == vertex_centered
      ? ivect (ceil (((rpos - origin) * scale                    ) / rvect(istride) - rvect(0.5))) * istride
      : ivect (ceil (((rpos - origin) * scale - rvect(bistride/2)) / rvect(istride) - rvect(0.5))) * istride - istride/2 + bistride/2;
#endif
    
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
  
  // Snap (enlarge) a bbox to the next coarser level, if desired
  ibbox
  snap_ibbox (ibbox const & ib,
              gh const & hh, int const rl)
  {
    DECLARE_CCTK_PARAMETERS;
    
    assert (not ib.empty());
    assert (rl > 0);
    
    if (not snap_to_coarse) return ib;
    
    ibbox const & base  = hh.baseextents.at(0).at(rl);
    ibbox const & cbase = hh.baseextents.at(0).at(rl-1);
    assert (all (cbase.stride() % base.stride() == 0));
    ivect const reffact = cbase.stride() / base.stride();
    
    ivect const lo  = ib.lower();
    ivect const up  = ib.upper();
    ivect const str = ib.stride();
    assert (all (str == base.stride()));
    
    if (veryverbose) {
      cout << "Snapping: coarse is " << cbase << ", current is " << base << "\n";
    }
    
#warning "TODO: shift/unshift boxes, because we are looking at grid points, not cell boundaries"

    return ib.expand(reffact-1).contracted_for(cbase).expanded_for(base);
  }
  
  
  
  void
  get_boundary_specification (jjvect & nboundaryzones,
                              jjvect & is_internal,
                              jjvect & is_staggered,
                              jjvect & shiftout)
  {
    if (CCTK_IsFunctionAliased ("MultiPatch_GetBoundarySpecification")) {
      assert (Carpet::map >= 0);
      CCTK_INT const ierr = MultiPatch_GetBoundarySpecification
        (Carpet::map, 2*dim,
         & nboundaryzones[0][0],
         & is_internal[0][0],
         & is_staggered[0][0],
         & shiftout[0][0]);
      assert (not ierr);
    } else if (CCTK_IsFunctionAliased ("GetBoundarySpecification")) {
      CCTK_INT const ierr = GetBoundarySpecification
        (2*dim,
         & nboundaryzones[0][0],
         & is_internal[0][0],
         & is_staggered[0][0],
         & shiftout[0][0]);
      assert (not ierr);
    } else {
      assert (0);
    }
  }
  
  void
  get_physical_boundary (rvect & physical_lower,
                         rvect & physical_upper,
                         rvect & spacing)
  {
    rvect interior_lower, interior_upper;
    rvect exterior_lower, exterior_upper;
    if (CCTK_IsFunctionAliased ("MultiPatch_GetDomainSpecification")) {
      assert (Carpet::map >= 0);
      CCTK_INT const ierr = MultiPatch_GetDomainSpecification
        (Carpet::map, dim,
         & physical_lower[0], & physical_upper[0],
         & interior_lower[0], & interior_upper[0],
         & exterior_lower[0], & exterior_upper[0],
         & spacing[0]);
      assert (not ierr);
    } else if (CCTK_IsFunctionAliased ("GetDomainSpecification")) {
      CCTK_INT const ierr = GetDomainSpecification
        (dim,
         & physical_lower[0], & physical_upper[0],
         & interior_lower[0], & interior_upper[0],
         & exterior_lower[0], & exterior_upper[0],
         & spacing[0]);
      assert (not ierr);
    } else {
      assert (0);
    }
  }
  
  void
  calculate_exterior_boundary (rvect const & physical_lower,
                               rvect const & physical_upper,
                               rvect       & exterior_lower,
                               rvect       & exterior_upper,
                               rvect const & spacing)
  {
    rvect interior_lower, interior_upper;
    if (CCTK_IsFunctionAliased ("MultiPatch_ConvertFromPhysicalBoundary")) {
      assert (Carpet::map >= 0);
      CCTK_INT const ierr = MultiPatch_ConvertFromPhysicalBoundary
        (Carpet::map, dim,
         & physical_lower[0], & physical_upper[0],
         & interior_lower[0], & interior_upper[0],
         & exterior_lower[0], & exterior_upper[0],
         & spacing[0]);
      assert (not ierr);
    } else if (CCTK_IsFunctionAliased ("ConvertFromPhysicalBoundary")) {
      CCTK_INT const ierr =
        ConvertFromPhysicalBoundary
        (dim,
         & physical_lower[0], & physical_upper[0],
         & interior_lower[0], & interior_upper[0],
         & exterior_lower[0], & exterior_upper[0],
         & spacing[0]);
      assert (not ierr);
    } else {
      assert (0);
    }
  }
  
  
  
  extern "C" {
    CCTK_INT
    CarpetRegrid2_Regrid (CCTK_POINTER_TO_CONST const cctkGH_,
                          CCTK_POINTER          const superregss_,
                          CCTK_POINTER          const regsss_,
                          CCTK_INT              const force);
    
    CCTK_INT
    CarpetRegrid2_RegridMaps (CCTK_POINTER_TO_CONST const cctkGH_,
                              CCTK_POINTER          const superregsss_,
                              CCTK_POINTER          const regssss_,
                              CCTK_INT              const force);
  }
  
  
  
  void
  Regrid (cGH const * const cctkGH,
          gh::rregs & regss)
  {
    DECLARE_CCTK_PARAMETERS;
    
    if (verbose or veryverbose) CCTK_INFO ("Regridding");
    
    assert (is_singlemap_mode());
    gh const & hh = * vhh.at (Carpet::map);
    dh const & dd = * vdd.at (Carpet::map);
    
    //
    // Find extent of domain
    //
    
    // This requires that CoordBase is used
#warning "TODO: check this (for Carpet, and maybe also for CartGrid3D)"
#warning "TODO: (the check that these two are consistent should be in Carpet)"
    
    jjvect nboundaryzones, is_internal;
    jjvect is_staggered, shiftout;
    get_boundary_specification (nboundaryzones, is_internal,
                                is_staggered, shiftout);
    
    rvect physical_lower, physical_upper;
    rvect spacing;
    get_physical_boundary (physical_lower, physical_upper, spacing);
    
    // Adapt spacing for convergence level
    spacing *= ipow ((CCTK_REAL) mgfact, basemglevel);
    
    rvect exterior_lower, exterior_upper;
    calculate_exterior_boundary (physical_lower, physical_upper,
                                 exterior_lower, exterior_upper,
                                 spacing);
    
    //
    // Determine which refinement levels may be changed
    //
    int min_rl = 1;             // we cannot change the coarsest level
    if (freeze_unaligned_levels or freeze_unaligned_parent_levels) {
      while (min_rl < int(regss.size())) {
        // Increase min_rl until we find a level that can be changed
#if 0
#warning "think about this a bit more"
#warning "use this taper-checking also in Comm.cc"
        bool in_sync = true;
        if (freeze_unaligned_parent_levels) {
          int const parent_do_every =
            ipow(mgfact, mglevel) *
            (maxtimereflevelfact / timereffacts.at(min_rl-1));
          in_sync =
            cctkGH->cctk_iteration == 0 or
            (cctkGH->cctk_iteration-1) % parent_do_every == 0;
        }
#else
        bool in_sync = true;
        if (freeze_unaligned_parent_levels) {
          // Assume that non-existing levels are always aligned
          if (min_rl < reflevels) {
            CCTK_REAL const mytime     = tt->get_time (mglevel, min_rl, 0);
            CCTK_REAL const parenttime = tt->get_time (mglevel, min_rl - 1, 0);
            CCTK_REAL const eps = 1.0e-12;
            in_sync =
              abs (mytime - parenttime) <= eps * abs (delta_time);
          }
        }
#endif
        int const do_every =
          ipow(mgfact, mglevel) *
          (maxtimereflevelfact / timereffacts.at(min_rl));
        bool const is_active =
          cctkGH->cctk_iteration == 0 or
          (cctkGH->cctk_iteration-1) % do_every == 0;
        if (is_active and in_sync) break;
        ++ min_rl;
      }
      if (verbose or veryverbose) {
        CCTK_VInfo (CCTK_THORNSTRING, "Regridding levels %d and up", min_rl);
      }
    }
    
    //
    // Calculate the union of the bounding boxes for all levels
    //
    
    // Find out whether the grid staggering corresponds to the
    // boundary staggering.  If there is a mismatch, then there cannot
    // be refinement near the corresponding boundary.
    b2vect const boundary_staggering_mismatch =
      xpose ((hh.refcent == vertex_centered) != (is_staggered == 0));
#warning "TODO: This is too strict"
    assert (all (all (not boundary_staggering_mismatch)));
    
    // This is the physical boundary
    rvect const origin (exterior_lower);
#warning "TODO: use hh.baseextents[0] [rl] instead of [0]"
    rvect const scale (rvect (hh.baseextents.at(0).at(0).stride()) / spacing);
    
    // This is the location of the outermost grid points.  For cell
    // centring, these are 1/2 grid spacing inside of the boundary.
    ivect const physical_ilower =
      rpos2ipos (physical_lower, origin, scale, hh, 0);
    ivect const physical_iupper =
      rpos2ipos1 (physical_upper, origin, scale, hh, 0);
    
    // The set of refined regions
    vector <ibset> regions (min_rl);
    
    // Set up coarse levels
    for (int rl = 0; rl < min_rl; ++ rl) {
      if (veryverbose) {
        cout << "Refinement level " << rl << ": will not be changed" << endl;
      }
      for (size_t c = 0; c < regss.at(rl).size(); ++ c) {
        regions.at(rl) += regss.at(rl).at(c).extent;
      }
      if (veryverbose) {
        cout << "Refinement level " << rl << ": regions are " << regions.at(rl) << endl;
      }
    }

    // Refine only patch 0
    if (Carpet::map == 0) {
      // Loop over all centres
      for (int n = 0; n < num_centres; ++ n) {
        centre_description centre (cctkGH, n);
        if (centre._active) {
          
          // Loop over all levels for this centre
          for (int rl = min_rl; rl < centre._num_levels; ++ rl) {
            
            if (veryverbose) {
              cout << "Refinement level " << rl << ": importing refined region settings..." << endl;
            }
            
            // Calculate a bbox for this region
            rvect const rmin = centre._position - centre._radius.at(rl);
            rvect const rmax = centre._position + centre._radius.at(rl);
            
            if (veryverbose) {
              cout << "Centre " << n+1 << " refinement level " << rl << ": coordinate region is (" << rmin << ":" << rmax << ")\n";
            }
            
            // Convert to an integer bbox
            ivect const istride = hh.baseextents.at(0).at(rl).stride();
            ivect const imin =
              rpos2ipos (rmin, origin, scale, hh, rl)
              - boundary_shiftout * istride;
            ivect const imax =
              rpos2ipos1 (rmax, origin, scale, hh, rl)
              + boundary_shiftout * istride;
            
            if (veryverbose) {
              cout << "Centre " << n+1 << " refinement level " << rl << ": integer region is (" << imin << ":" << imax << ")\n";
            }
            
            ibbox const region =
              snap_ibbox (ibbox (imin, imax, istride), hh, rl);
            
            if (veryverbose) {
              cout << "Centre " << n+1 << " refinement level " << rl << ": snapped integer region is " << region << "\n";
            }
            
            // Add this region to the list of regions
            if (static_cast <int> (regions.size()) < rl+1) {
              regions.resize (rl+1);
            }
            regions.at(rl) |= region;
            
            if (veryverbose) {
              cout << "Refinement level " << rl << ": preliminary regions are " << regions.at(rl) << endl;
            }
            
          } // for rl
          
        } // if centre is active
      }   // for n
    }     // if map==0
    
    
    
    regss.resize (regions.size());
    
    assert (regions.size() > 0);
    for (int rl = regions.size() - 1; rl >= min_rl; -- rl) {
      
      // Sanity check
      assert (not regions.at(rl).empty());
      
      
      
      //
      // Ensure that this grid contains the next finer grid
      //
      if (ensure_proper_nesting) {
        if (rl < int(regions.size()) - 1) {
          
          if (veryverbose) {
            cout << "Refinement level " << rl << ": ensuring proper nesting..." << endl;
          }
          
          assert (not regions.at(rl).empty());
          ibbox const coarse0 = * regions.at(rl).begin();
          
          // This is the location of the outermost grid points.  For
          // cell centring, these are 1/2 grid spacing inside of the
          // boundary.
          ivect const level_physical_ilower =
            rpos2ipos (physical_lower, origin, scale, hh, rl);
          ivect const level_physical_iupper =
            rpos2ipos1 (physical_upper, origin, scale, hh, rl);
          
          i2vect const fdistance = dd.ghost_widths.at(rl);
          i2vect const cdistance =
            i2vect (min_distance + dd.prolongation_stencil_size(rl));
          
          for (ibset::const_iterator ibb = regions.at(rl+1).begin();
               ibb != regions.at(rl+1).end();
               ++ ibb)
          {
            ibbox const & fbb = * ibb;
            
            bvect const lower_is_outer = fbb.lower() <= level_physical_ilower;
            bvect const upper_is_outer = fbb.upper() >= level_physical_iupper;
            b2vect const ob (lower_is_outer, upper_is_outer);
            
            ibbox const ebb = fbb.expand (i2vect (not ob) * fdistance);
            ibbox const cbb = ebb.expanded_for (coarse0);
            ibbox const ecbb = cbb.expand (i2vect (not ob) * cdistance);
            
            // Enlarge this level
            regions.at(rl) |= snap_ibbox (ecbb, hh, rl);
          }
          
          if (veryverbose) {
            cout << "Refinement level " << rl << ": enlarged regions to " << regions.at(rl) << endl;
          }
          
        }
      } // if ensure proper nesting
      
      
      
      //
      // Add buffer zones
      //
      {
        if (veryverbose) {
          cout << "Refinement level " << rl << ": adding buffer zones..." << endl;
        }

        ibset buffered;
        for (ibset::const_iterator ibb = regions.at(rl).begin();
             ibb != regions.at(rl).end();
             ++ ibb)
        {
          ibbox const & bb = * ibb;
          
          ibbox const bbb = bb.expand (dd.buffer_widths.at(rl));

          buffered |= bbb;
        }
        regions.at(rl) = buffered;
        
        if (veryverbose) {
          cout << "Refinement level " << rl << ": enlarged regions to " << regions.at(rl) << endl;
        }
      }
      
      
      
      //
      // Check whether it would be worthwhile to combine all regions
      // into a single region
      //
      // TODO: Check this also for pairs of regions
      //
      // TODO: Check after clipping
      {
        if (veryverbose) {
          cout << "Refinement level " << rl << ": checking whether regions should be combined..." << endl;
        }

        ibbox single;
        for (ibset::const_iterator ibb = regions.at(rl).begin();
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
          regions.at(rl) = ibset (single);
          if (veryverbose) {
            cout << "Refinement level " << rl << ": combining regions to " << regions.at(rl) << endl;
          }
        } else {
          if (veryverbose) {
            cout << "Refinement level " << rl << ": not combining" << endl;
          }
        }
      }
      
      
      
      // Find the location of the outer boundary
      if (veryverbose) {
        cout << "Refinement level " << rl << ": determining outer boundary..." << endl;
      }
      
      rvect const level_physical_lower = physical_lower;
      rvect const level_physical_upper = physical_upper;
      rvect const level_spacing = spacing / rvect (hh.reffacts.at(rl));
      if (veryverbose) {
        cout << "Refinement level " << rl << ": physical coordinate boundary is at " << r2vect(level_physical_lower, level_physical_upper) << endl;
        cout << "Refinement level " << rl << ": spacing is " << level_spacing << endl;
      }
#if 0
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
#endif
      rvect level_exterior_lower, level_exterior_upper;
      calculate_exterior_boundary (level_physical_lower, level_physical_upper,
                                   level_exterior_lower, level_exterior_upper,
                                   level_spacing);
      if (veryverbose) {
        cout << "Refinement level " << rl << ": exterior coordinate boundary is at " << r2vect(level_exterior_lower, level_exterior_upper) << endl;
      }
      
      ibbox const & baseextent = hh.baseextents.at(0).at(rl);
      ivect const istride = baseextent.stride();
      if (veryverbose) {
        cout << "Refinement level " << rl << ": stride is " << istride << endl;
      }
      
      // This is the location of the outermost grid points.  For cell
      // centring, these are 1/2 grid spacing inside of the boundary.
      ivect const level_physical_ilower =
        rpos2ipos (physical_lower, origin, scale, hh, rl);
      ivect const level_physical_iupper =
        rpos2ipos1 (physical_upper, origin, scale, hh, rl);
      if (veryverbose) {
        cout << "Refinement level " << rl << ": physical boundary is at " << i2vect(level_physical_ilower, level_physical_iupper) << endl;
        cout << "Refinement level " << rl << ": reconstructed physical coordinate boundary is at "
             << r2vect(ipos2rpos(level_physical_ilower - (hh.refcent == cell_centered ? istride/2 : 0), origin, scale, hh, rl),
                       ipos2rpos(level_physical_iupper + (hh.refcent == cell_centered ? istride/2 : 0), origin, scale, hh, rl)) << endl;
      }
      
      ivect const level_exterior_ilower =
        rpos2ipos (level_exterior_lower, origin, scale, hh, rl);
      ivect const level_exterior_iupper =
        rpos2ipos1 (level_exterior_upper, origin, scale, hh, rl);
      assert (all (level_exterior_ilower >= baseextent.lower()));
      assert (all (level_exterior_iupper <= baseextent.upper()));
      if (veryverbose) {
        cout << "Refinement level " << rl << ": exterior boundary is at " << i2vect(level_exterior_ilower, level_exterior_iupper) << endl;
        cout << "Refinement level " << rl << ": reconstructed exterior coordinate boundary is at "
             << r2vect(ipos2rpos(level_exterior_ilower, origin, scale, hh, rl),
                       ipos2rpos(level_exterior_iupper, origin, scale, hh, rl)) << endl;
      }
      
      // Find the minimum necessary distance away from the outer
      // boundary due to buffer and ghost zones.  This is e.g. the
      // distance that the lower boundary of a bbox has to have from
      // the lower boundary.  This is in terms of grid points.
      i2vect const min_bnd_dist_away = dd.ghost_widths.at(rl);
      // Find the minimum necessary distance from the outer boundary
      // due to buffer and ghost zones.  This is e.g. the distance
      // that the upper boundary of a bbox has to have from the lower
      // boundary.  This is in terms of grid points.
      i2vect const min_bnd_dist_incl = dd.ghost_widths.at(rl);
      // TODO: The above is required only near symmetry boundaries.
      
      
      
      //
      // Make the boxes rotating-90 symmetric
      //
      if (symmetry_rotating90) {
        if (veryverbose) {
          cout << "Refinement level " << rl << ": making regions rotating-90 symmetric" << endl;
        }
        
        ibset added;
        for (ibset::const_iterator ibb = regions.at(rl).begin();
             ibb != regions.at(rl).end();
             ++ ibb)
        {
          ibbox const & bb = * ibb;
          if (veryverbose) {
            cout << "   considering box " << bb << "..." << endl;
          }
          
          bvect const lower_is_outside_lower =
            bb.lower() - min_bnd_dist_away[0] * bb.stride() <= level_physical_ilower;
          
          // Treat both x and y directions
          for (int dir=0; dir<=1; ++dir) {
            if (veryverbose) {
              cout << "   symmetrising in " << "xy"[dir] << " direction..." << endl;
            }
            if (lower_is_outside_lower[dir]) {
              ivect const ilo = bb.lower();
              ivect const iup = bb.upper();
              ivect const istr = bb.stride();
              
              // Origin
              rvect const axis (physical_lower[0],
                                physical_lower[1],
                                physical_lower[2]); // z component is unused
              ivect const iaxis0 =
                rpos2ipos (axis, origin, scale, hh, rl);
              assert (all (iaxis0 % istr == 0));
              ivect const iaxis1 =
                rpos2ipos1 (axis, origin, scale, hh, rl);
              assert (all (iaxis1 % istr == 0));
              ivect const offset = iaxis1 - iaxis0;
              assert (all (offset % istr == 0));
              assert (all (offset >= 0 and offset < 2*istr));
              assert (all ((iaxis0 + iaxis1 - offset) % (2*istr) == 0));
              ivect const iaxis = (iaxis0 + iaxis1 - offset) / 2;
              // negated (reflected) domain boundaries
              ivect const neg_ilo = (2*iaxis+offset) - ilo;
              ivect const neg_iup = (2*iaxis+offset) - iup;
              // offset to add when permuting directions
              ivect const permute01 (-iaxis[0]+iaxis[1], -iaxis[1]+iaxis[0], 0);
              
              // Rotate 90 degrees about z axis
              ivect new_ilo, new_iup;
              if (dir==0) {
                // rotate clockwise
                new_ilo = ivect (ilo[1], neg_iup[0], ilo[2]) + permute01;
                new_iup = ivect (iup[1], neg_ilo[0], iup[2]) + permute01;
              }
              if (dir==1) {
                // rotate counterclockwise
                new_ilo = ivect (neg_iup[1], ilo[0],  ilo[2]) + permute01;
                new_iup = ivect (neg_ilo[1],  iup[0],  iup[2]) + permute01;
              }
              ivect const new_istr (istr);
              
              ibbox const new_bb (new_ilo, new_iup, new_istr);
              // Will be clipped later
              //assert (new_bb.is_contained_in (baseextent));
              
              added |= new_bb;
              if (veryverbose) {
                cout << "      adding box " << new_bb << endl;
              }
            }
          }
        }
          
        regions.at(rl) |= added;
        
        if (veryverbose) {
          cout << "Refinement level " << rl << ": symmetrised regions are " << regions.at(rl) << endl;
        }
      }
      
      
      
      //
      // Make the boxes rotating-180 symmetric
      //
      if (symmetry_rotating180) {
        if (veryverbose) {
          cout << "Refinement level " << rl << ": making regions rotating-180 symmetric" << endl;
        }
        
        ibset added;
        for (ibset::const_iterator ibb = regions.at(rl).begin();
             ibb != regions.at(rl).end();
             ++ ibb)
        {
          ibbox const & bb = * ibb;
          if (veryverbose) {
            cout << "   considering box " << bb << "..." << endl;
          }
          
          bvect const lower_is_outside_lower =
            bb.lower() - min_bnd_dist_away[0] * bb.stride() <= level_physical_ilower;
          
          // Treat x direction
          int const dir = 0;
          if (veryverbose) {
            cout << "   symmetrising in " << "x"[dir] << " direction..." << endl;
          }
          if (lower_is_outside_lower[dir]) {
            ivect const ilo = bb.lower();
            ivect const iup = bb.upper();
            ivect const istr = bb.stride();
            assert (istr[0] == istr[1]);
            
            // Origin
            assert (hh.refcent == vertex_centered or all (istr % 2 == 0));
            rvect const axis (physical_lower[0],
                              (physical_lower[1] + physical_upper[1]) / 2,
                              physical_lower[2]); // z component is unused
            ivect const iaxis0 =
              rpos2ipos (axis, origin, scale, hh, rl);
            assert (all ((iaxis0 - baseextent.lower()) % istr == 0));
            ivect const iaxis1 =
              rpos2ipos1 (axis, origin, scale, hh, rl);
            assert (all ((iaxis1 - baseextent.lower()) % istr == 0));
            ivect const offset = iaxis1 - iaxis0;
            assert (all (offset % istr == 0));
            if (hh.refcent == vertex_centered) {
              assert (all (offset >= 0 and offset < 2*istr));
              assert (all ((iaxis0 + iaxis1 - offset) % (2*istr) == 0));
            } else {
              // The offset may be negative because both boundaries
              // are shifted inwards by 1/2 grid spacing, and
              // therefore iaxis0 < iaxis1 + istr
              assert (all (offset >= -istr and offset < istr));
              assert (all ((iaxis0 + iaxis1 - offset) % (2*istr) == istr));
              assert (all (istr % 2 == 0));
            }
            ivect const iaxis = (iaxis0 + iaxis1 - offset) / 2;
            ivect const neg_ilo = (2*iaxis+offset) - ilo;
            ivect const neg_iup = (2*iaxis+offset) - iup;
            
            // Rotate 180 degrees about z axis
            ivect const new_ilo (neg_iup[0], neg_iup[1], ilo[2]);
            ivect const new_iup (neg_ilo[0], neg_ilo[1], iup[2]);
            ivect const new_istr (istr);
            
            ibbox const new_bb (new_ilo, new_iup, new_istr);
            // Will be clipped later
            //assert (new_bb.is_contained_in (baseextent));
            
            added |= new_bb;
            if (veryverbose) {
              cout << "      adding box " << new_bb << endl;
            }
          }
        }
        
        regions.at(rl) |= added;
        
        if (veryverbose) {
          cout << "Refinement level " << rl << ": symmetrised regions are " << regions.at(rl) << endl;
        }
      } // if symmetry_rotating180
      
      
      
      //
      // Clip at the outer boundary
      //
      {
        if (veryverbose) {
          cout << "Refinement level " << rl << ": clipping at outer boundary..." << endl;
        }
        
        ibset clipped;
        for (ibset::const_iterator ibb = regions.at(rl).begin();
             ibb != regions.at(rl).end();
             ++ ibb)
        {
          ibbox const & bb = * ibb;
          if (veryverbose) {
            cout << "   considering box " << bb << "..." << endl;
          }
          
          // Clip boxes that extend outside the boundary.  Enlarge
          // boxes that are inside but too close to the outer
          // boundary.
          bvect const lower_is_outside_lower =
            bb.lower() - min_bnd_dist_away[0] * bb.stride() <= level_physical_ilower;
          // Remove bboxes that are completely outside.
          bvect const upper_is_outside_lower =
            bb.upper() < level_physical_ilower;
          // Enlarge bboxes that extend not far enough inwards.
          bvect const upper_is_almost_outside_lower =
            bb.upper() < level_physical_ilower + min_bnd_dist_incl[0] * bb.stride();
          
          // Ditto for the upper boundary.
          bvect const upper_is_outside_upper =
            bb.upper() + min_bnd_dist_away[1] * bb.stride() >= level_physical_iupper;
          bvect const lower_is_outside_upper =
            bb.lower() > level_physical_iupper;
          bvect const lower_is_almost_outside_upper =
            bb.lower() > level_physical_iupper - min_bnd_dist_incl[1] * bb.stride();
          
          assert (not any (lower_is_almost_outside_upper and
                           lower_is_outside_lower));
          assert (not any (upper_is_almost_outside_lower and
                           upper_is_outside_upper));
          
          if (any (upper_is_outside_lower or lower_is_outside_upper)) {
            // The box is completely outside.  Ignore it.
            if (veryverbose) {
              cout << "      box is completely outside; dropping it" << endl;
            }
            continue;
          }
          
          if (any ((lower_is_outside_lower and
                    boundary_staggering_mismatch[0]) or
                   (upper_is_outside_upper and
                    boundary_staggering_mismatch[1])))
          {
            ostringstream msg;
            msg << "Level " << rl << " of the refinement hierarchy has inconsistent bountary staggering."
                << "  The refined region extends up to the boundary, but the staggering of the boundary is different from the staggering of the mesh refinement."
                << "  lower_is_outside_lower=" << lower_is_outside_lower
                << "  upper_is_outside_upper=" << upper_is_outside_upper
                << "  boundary_staggering_mismatch=" << boundary_staggering_mismatch
                << "  level_physical_ilower=" << level_physical_ilower
                << "  level_physical_iupper=" << level_physical_iupper
                << "  baseextent=" << baseextent;
            CCTK_WARN (CCTK_WARN_ABORT, msg.str().c_str());
          }
          
          ibbox const clipped_bb
            (either (lower_is_outside_lower,
                     level_exterior_ilower,
                     either (lower_is_almost_outside_upper,
                             (level_physical_iupper -
                              min_bnd_dist_incl[1] * bb.stride()),
                             bb.lower())),
             either (upper_is_outside_upper,
                     level_exterior_iupper,
                     either (upper_is_almost_outside_lower,
                             (level_physical_ilower +
                              min_bnd_dist_incl[0] * bb.stride()),
                             bb.upper())),
             bb.stride());
          if (not clipped_bb.is_contained_in (baseextent)) {
            ostringstream msg;
            msg << "Level " << rl << " of the refinement hierarchy is not contained in the simulation domain."
                << "  (There may be too many ghost or buffer zones.)"
                << "  One bbox is " << clipped_bb << "."
                << "  lower_is_outside_lower=" << lower_is_outside_lower
                << "  upper_is_outside_upper=" << upper_is_outside_upper
                << "  lower_is_almost_outside_upper=" << lower_is_almost_outside_upper
                << "  upper_is_almost_outside_lower=" << upper_is_almost_outside_lower
                << "  level_exterior_ilower=" << level_exterior_ilower
                << "  level_exterior_iupper=" << level_exterior_iupper
                << "  level_physical_ilower=" << level_physical_ilower
                << "  level_physical_iupper=" << level_physical_iupper
                << "  baseextent=" << baseextent;
            CCTK_WARN (CCTK_WARN_ABORT, msg.str().c_str());
          }
          assert (clipped_bb.is_contained_in (baseextent));
          
          clipped |= clipped_bb;
          if (veryverbose) {
            cout << "      Clipped box is " << clipped_bb << endl;
          }
          
        } // for ibb
        
        regions.at(rl) = clipped;
        if (veryverbose) {
          cout << "Refinement level " << rl << ": clipped regions are " << regions.at(rl) << endl;
        }
      }
      
      
      
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
        for (ibset::const_iterator ibb = regions.at(rl).begin();
             ibb != regions.at(rl).end();
             ++ ibb)
        {
          ibbox const & bb = * ibb;
          assert (bb.is_contained_in (hh.baseextents.at(0).at(rl)));
          
          bvect const lower_is_outer = bb.lower() <= level_physical_ilower;
          bvect const upper_is_outer = bb.upper() >= level_physical_iupper;
          
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
    
    //
    // Check whether all grids are contained in the next coarser grid
    //
    for (int rl = regions.size() - 1; rl >= min_rl; -- rl) {
      
      assert (not regions.at(rl-1).empty());
      ibbox const coarse0 = * regions.at(rl-1).begin();
      
      // This is the location of the outermost grid points.  For cell
      // centring, these are 1/2 grid spacing inside of the boundary.
      ivect const level_physical_ilower =
        rpos2ipos (physical_lower, origin, scale, hh, rl);
      ivect const level_physical_iupper =
        rpos2ipos1 (physical_upper, origin, scale, hh, rl);
      
      i2vect const fdistance = dd.ghost_widths.at(rl);
      i2vect const cdistance =
        i2vect (min_distance + dd.prolongation_stencil_size(rl));
      
      bool is_properly_nested = true;
      
      for (ibset::const_iterator ibb = regions.at(rl).begin();
           ibb != regions.at(rl).end();
           ++ ibb)
      {
        ibbox const & fbb = * ibb;
        
        bvect const lower_is_outer = fbb.lower() <= level_physical_ilower;
        bvect const upper_is_outer = fbb.upper() >= level_physical_iupper;
        b2vect const ob (lower_is_outer, upper_is_outer);
        
        ibbox const cbb = fbb.expanded_for (coarse0);
        
        is_properly_nested = is_properly_nested and cbb <= regions.at(rl-1);
      }
      
      if (not is_properly_nested) {
        ostringstream msg;
        msg << "Level " << rl << " of the refinement hierarchy is not properly nested.  It is not contained in level " << (rl-1) << ".";
        CCTK_WARN (CCTK_WARN_ALERT, msg.str().c_str());
      }
    } // for rl
    
  }
  
  
  
  CCTK_INT
  CarpetRegrid2_Regrid (CCTK_POINTER_TO_CONST const cctkGH_,
                        CCTK_POINTER          const superregss_,
                        CCTK_POINTER          const regsss_,
                        CCTK_INT              const force)
  {
    cGH const * const cctkGH = static_cast <cGH const *> (cctkGH_);
    
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;
    
    char const * const where = "CarpetRegrid2::Regrid";
    static Carpet::Timer timer (where);
    timer.start();
    
    assert (is_singlemap_mode());
    
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
    
    if (verbose or veryverbose) {
      if (do_recompose) {
        for (int n = 0; n < num_centres; ++ n) {
          if (active[n]) {
            CCTK_VInfo (CCTK_THORNSTRING,
                        "Centre %d is at position [%g,%g,%g] with %d levels",
                        n + 1,
                        static_cast <double> (position_x[n]),
                        static_cast <double> (position_y[n]),
                        static_cast <double> (position_z[n]),
                        static_cast <int> (num_levels[n]));
          } else {
            CCTK_VInfo (CCTK_THORNSTRING,
                        "Centre %d is not active", n + 1);
          }
        }
      }
    }
    
    if (not force and do_recompose and * last_iteration != -1) {
      
      int lsh[2];
      getvectorindex2 (cctkGH, "CarpetRegrid2::radii", lsh);
      
      // Regrid only if the regions have changed sufficiently
      do_recompose = false;
      for (int n = 0; n < num_centres; ++ n) {
        
        // Regrid if a region became active or inactive
        do_recompose = active[n] != old_active[n];
        if (do_recompose) break;
        
        // Check only active regions
        if (not active[n]) continue;
        
        // Regrid if the number of levels changed
        do_recompose = num_levels[n] != old_num_levels[n];
        if (do_recompose) break;
        
        // Regrid if the positions have changed sufficiently
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
        do_recompose = dist2 > pow (mindist, 2);
        if (do_recompose) break;
        
        // Regrid if the radii have changed sufficiently
        for (int rl = 1; rl < num_levels[n]; ++ rl) {
          int const ind = index2 (lsh, rl, n);
          CCTK_REAL const rx = radius_x[ind] < 0 ? radius[ind] : radius_x[ind];
          CCTK_REAL const ry = radius_y[ind] < 0 ? radius[ind] : radius_y[ind];
          CCTK_REAL const rz = radius_z[ind] < 0 ? radius[ind] : radius_z[ind];
          rvect const rad (rx, ry, rz);
          rvect const oldrad
            (old_radius_x[ind], old_radius_y[ind], old_radius_z[ind]);
          CCTK_REAL const drfac = 
	    (sqrt (sum (ipow (rad - oldrad, 2))))/(sqrt (sum (ipow (oldrad, 2))));
          CCTK_REAL mindrfac;
          switch (n) {
          case 0: mindrfac = radius_rel_change_threshold_1; break;
          case 1: mindrfac = radius_rel_change_threshold_2; break;
          case 2: mindrfac = radius_rel_change_threshold_3; break;
          case 3: mindrfac = radius_rel_change_threshold_4; break;
          case 4: mindrfac = radius_rel_change_threshold_5; break;
          case 5: mindrfac = radius_rel_change_threshold_6; break;
          case 6: mindrfac = radius_rel_change_threshold_7; break;
          case 7: mindrfac = radius_rel_change_threshold_8; break;
          case 8: mindrfac = radius_rel_change_threshold_9; break;
          case 9: mindrfac = radius_rel_change_threshold_10; break;
          default: assert (0);
          }
          do_recompose = drfac > mindrfac;
          if (do_recompose) break;
        } // for rl
        if (do_recompose) break;
        
      } // for n
      if (verbose or veryverbose) {
        if (not do_recompose) {
          CCTK_INFO
            ("Refined regions have not changed sufficiently; skipping regridding");
        }
      }
    }
    
    if (do_recompose) {
      * last_iteration = cctk_iteration;
      * last_map = Carpet::map;
    }
    
    if (do_recompose) {
      
      gh::rregs & superregss = * static_cast <gh::rregs *> (superregss_);
      gh::mregs & regsss     = * static_cast <gh::mregs *> (regsss_);
      
      Regrid (cctkGH, superregss);
      
      // Make multiprocessor aware
      vector <vector <region_t> > regss (superregss.size());
      for (size_t rl = 0; rl < regss.size(); ++ rl) {
        SplitRegions (cctkGH, superregss.at(rl), regss.at(rl));
      } // for rl
      
      // Make multigrid aware
      MakeMultigridBoxes (cctkGH, Carpet::map, regss, regsss);
      
      // Remember current positions
      for (int n = 0; n < num_centres; ++ n) {
        old_active[n] = active[n];
        
        old_num_levels[n] = num_levels[n];
        
        old_position_x[n] = position_x[n];
        old_position_y[n] = position_y[n];
        old_position_z[n] = position_z[n];
        
        old_radius_x[n] = radius_x[n] < 0 ? radius[n] : radius_x[n];
        old_radius_y[n] = radius_y[n] < 0 ? radius[n] : radius_y[n];
        old_radius_z[n] = radius_z[n] < 0 ? radius[n] : radius_z[n];
      }
      
    } // if do_recompose
    
    timer.stop();
    
    return do_recompose;
  }
  
  
  
  CCTK_INT
  CarpetRegrid2_RegridMaps (CCTK_POINTER_TO_CONST const cctkGH_,
                            CCTK_POINTER          const superregsss_,
                            CCTK_POINTER          const regssss_,
                            CCTK_INT              const force)
  {
    cGH const * const cctkGH = static_cast <cGH const *> (cctkGH_);
    
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;
    
    char const * const where = "CarpetRegrid2::RegridMaps";
    static Carpet::Timer timer (where);
    timer.start();
    
    assert (is_level_mode());
    
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
    
    if (verbose or veryverbose) {
      if (do_recompose) {
        for (int n = 0; n < num_centres; ++ n) {
          if (active[n]) {
            CCTK_VInfo (CCTK_THORNSTRING,
                        "Centre %d is at position [%g,%g,%g] with %d levels",
                        n + 1,
                        static_cast <double> (position_x[n]),
                        static_cast <double> (position_y[n]),
                        static_cast <double> (position_z[n]),
                        static_cast <int> (num_levels[n]));
          } else {
            CCTK_VInfo (CCTK_THORNSTRING,
                        "Centre %d is not active", n + 1);
          }
        }
      }
    }
    
    if (not force and do_recompose and * last_iteration != -1) {
      
      int lsh[2];
      getvectorindex2 (cctkGH, "CarpetRegrid2::radii", lsh);
      
      // Regrid only if the regions have changed sufficiently
      do_recompose = false;
      for (int n = 0; n < num_centres; ++ n) {
        
        // When debugging, sneakily add a new level, but skip the
        // initial regrid, and the regrid before the first time step
        if (add_levels_automatically and cctk_iteration > 1) {
          num_levels[n] = min (num_levels[n] + 1, maxreflevels);
          CCTK_VInfo (CCTK_THORNSTRING,
                      "Increasing number of levels of centre %d to %d (it=%d)",
                      n + 1, 
                      static_cast <int> (num_levels[n]),
                      cctk_iteration);
        }
        
        // Regrid if a region became active or inactive
        do_recompose = active[n] != old_active[n];
        if (do_recompose) break;
        
        // Check only active regions
        if (not active[n]) continue;
        
        // Regrid if the number of levels changed
        do_recompose = num_levels[n] != old_num_levels[n];
        if (do_recompose) break;
        
        // Regrid if the positions have changed sufficiently
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
        do_recompose = dist2 > pow (mindist, 2);
        if (do_recompose) break;
        
        // Regrid if the radii have changed sufficiently
        for (int rl = 1; rl < num_levels[n]; ++ rl) {
          int const ind = index2 (lsh, rl, n);
          CCTK_REAL const rx = radius_x[ind] < 0 ? radius[ind] : radius_x[ind];
          CCTK_REAL const ry = radius_y[ind] < 0 ? radius[ind] : radius_y[ind];
          CCTK_REAL const rz = radius_z[ind] < 0 ? radius[ind] : radius_z[ind];
          rvect const rad (rx, ry, rz);
          rvect const oldrad
            (old_radius_x[ind], old_radius_y[ind], old_radius_z[ind]);
          CCTK_REAL const drfac = 
	    (sqrt (sum (ipow (rad - oldrad, 2))))/(sqrt (sum (ipow (oldrad, 2))));
          CCTK_REAL mindrfac;
          switch (n) {
          case 0: mindrfac = radius_rel_change_threshold_1; break;
          case 1: mindrfac = radius_rel_change_threshold_2; break;
          case 2: mindrfac = radius_rel_change_threshold_3; break;
          case 3: mindrfac = radius_rel_change_threshold_4; break;
          case 4: mindrfac = radius_rel_change_threshold_5; break;
          case 5: mindrfac = radius_rel_change_threshold_6; break;
          case 6: mindrfac = radius_rel_change_threshold_7; break;
          case 7: mindrfac = radius_rel_change_threshold_8; break;
          case 8: mindrfac = radius_rel_change_threshold_9; break;
          case 9: mindrfac = radius_rel_change_threshold_10; break;
          default: assert (0);
          }
          do_recompose = drfac > mindrfac;
          if (do_recompose) break;
        } // for rl
        if (do_recompose) break;
        
      } // for n
      if (verbose or veryverbose) {
        if (not do_recompose) {
          CCTK_INFO
            ("Refined regions have not changed sufficiently; skipping regridding");
        }
      }
    }
    
    if (do_recompose) {
      * last_iteration = cctk_iteration;
    }
    
    if (do_recompose) {
      
      vector <gh::rregs> & superregsss =
        * static_cast <vector <gh::rregs>  *> (superregsss_);
      vector <gh::mregs> & regssss =
        * static_cast <vector <gh::mregs>  *> (regssss_);
      
      BEGIN_MAP_LOOP (cctkGH, CCTK_GF) {
        Regrid (cctkGH, superregsss.at(Carpet::map));
      } END_MAP_LOOP;
      
      vector< vector <vector <region_t> > > regsss (maps);
      
      // Count levels
      vector <int> rls (maps);
      for (int m = 0; m < maps; ++ m) {
        rls.at(m) = superregsss.at(m).size();
      }
      int maxrl = 0;
      for (int m = 0; m < maps; ++ m) {
        maxrl = max (maxrl, rls.at(m));
      }
      // All maps must have the same number of levels
      for (int m = 0; m < maps; ++ m) {
        superregsss.at(m).resize (maxrl);
        regsss.at(m).resize (maxrl);
      }
      
      // Make multiprocessor aware
      for (int rl = 0; rl < maxrl; ++ rl) {
        vector <vector <region_t> > superregss (maps);
        for (int m = 0; m < maps; ++ m) {
          superregss.at(m) = superregsss.at(m).at(rl);
        }
        vector <vector <region_t> > regss (maps);
        SplitRegionsMaps (cctkGH, superregss, regss);
        
        for (int m = 0; m < maps; ++ m) {
          superregsss.at(m).at(rl) = superregss.at(m);
          regsss.at(m).at(rl) = regss.at(m);
        }
      } // for rl
      
      // Make multigrid aware
      MakeMultigridBoxesMaps (cctkGH, regsss, regssss);
      
      int lsh[2];
      getvectorindex2 (cctkGH, "CarpetRegrid2::radii", lsh);
      
      // Remember current positions
      for (int n = 0; n < num_centres; ++ n) {
        old_active[n] = active[n];
        
        old_num_levels[n] = num_levels[n];
        
        old_position_x[n] = position_x[n];
        old_position_y[n] = position_y[n];
        old_position_z[n] = position_z[n];
        
	for (int rl = 1; rl < num_levels[n]; ++ rl) {
	  int const ind = index2 (lsh, rl, n);
	  old_radius_x[ind] = radius_x[ind] < 0 ? radius[ind] : radius_x[ind];
	  old_radius_y[ind] = radius_y[ind] < 0 ? radius[ind] : radius_y[ind];
	  old_radius_z[ind] = radius_z[ind] < 0 ? radius[ind] : radius_z[ind];
	}
      }
      
    } // if do_recompose
    
    timer.stop();
    
    return do_recompose;
  }
  
  
  
} // namespace CarpetRegrid2
