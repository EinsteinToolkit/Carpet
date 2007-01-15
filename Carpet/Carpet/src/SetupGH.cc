#include <cassert>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <sstream>
#include <stack>
#include <vector>

#include <cctk.h>
#include <cctk_Parameters.h>

#include <util_ErrorCodes.h>
#include <util_Table.h>

#include <bbox.hh>
#include <defs.hh>
#include <dist.hh>
#include <ggf.hh>
#include <gh.hh>
#include <region.hh>
#include <vect.hh>

#include "carpet.hh"



namespace Carpet {
  
  using namespace std;
  
  
  
  static void setup_model_information (cGH * cctkGH);
  static void setup_multigrid_information (cGH * cctkGH);
  static void setup_refinement_information ();
  static void setup_map_information ();
  static void setup_domain_extents (cGH const * cctkGH);
  static void allocate_grid_hierarchy (int m, ivect const & npoints);
  static void
  allocate_data_hierarchy (int m,
                           ivect const & lghosts, ivect const & ughosts);
  static void allocate_time_hierarchy (int m);
  static void setup_grid_hierarchy (cGH const * cctkGH);
  static void
  set_base_extent (int m,
                   vector<vector<region_t> > & regss);
  
  static void allocate_group_data (cGH const * cctkGH);
  static void
  allocate_group_hierarchies (int group,
                              ivect const & sizes,
                              ivect const & lghosts,
                              ivect const & ughosts);
  static void
  setup_group_grid_hierarchy (cGH const * cctkGH,
                              int group, cGroup const & gdata,
                              ivect const & convpowers,
                              ivect const & convoffsets);
  static void
  initialise_group_info (cGH const * cctkGH,
                         int group, cGroup const & gdata);
  static void set_state (cGH * cctkGH);
  static void enable_storage_for_all_groups (cGH const * cctkGH);
  
  
  
  static ivect get_ghostzones ();
  static ivect get_npoints ();
  static void  get_boundary_specification (int m);
  static void
  get_domain_specification (int m,
                            ivect const & npoints,
                            rvect & physical_min,
                            rvect & physical_max,
                            rvect & base_spacing);
  static void
  adapt_domain_specification (int m,
                              rvect const & physical_min,
                              rvect const & physical_max,
                              rvect const & base_spacing,
                              rvect & exterior_min,
                              rvect & exterior_max,
                              rvect & spacing);
  static void
  calculate_grid_points (int m,
                         ivect const & lghosts,
                         ivect const & ughosts,
                         rvect const & exterior_min,
                         rvect const & exterior_max,
                         rvect const & spacing,
                         ivect & npoints);
  static void
  find_processor_decomposition
  (cGH const * cctkGH,
   vector<vector<vector<region_t> > >          & regsss,
   vector<vector<vector<vector<region_t> > > > & regssss);
  
  static void
  get_group_size (int group, cGroup const & gdata,
                  ivect & sizes,
                  ivect & lghosts, ivect & ughosts);
  static void
  adapt_group_size_mglevel (int group, cGroup const & gdata,
                            ivect & sizes,
                            ivect & convpowers,
                            ivect & convoffsets);
  static void
  get_convergence_options (int group, cGroup const & gdata,
                           ivect & convpowers,
                           ivect & convoffsets);
  static void
  adapt_group_size_disttype (cGH const * cctkGH,
                             int group, cGroup const & gdata,
                             ivect & sizes,
                             ivect const & lghosts, ivect const & ughosts);
  static void output_group_statistics (cGH const * cctkGH);
  static operator_type
  get_transport_operator (cGH const * cctkGH,
                          int group, cGroup const & gdata);
  static bool
  can_transfer_variable_type (cGH const * cctkGH,
                              int group, cGroup const & gdata);
  
  
  
  static void ensure_CartGrid3D_type ();
  static void ensure_CartGrid3D_avoid_origin ();
  static void ensure_ReflectionSymmetry_avoid_origin ();
  static void
  ensure_ghostzones (int m,
                     ivect const & lghosts, ivect const & ughosts);
  static void ensure_group_options (int group, cGroup const & gdata);
  


  
  
  
  void *
  SetupGH (tFleshConfig * const fc,
           int const convLevel,
           cGH * const cctkGH)
  {
    DECLARE_CCTK_PARAMETERS;
    
    // Initialise current position (must be the very first thing,
    // before the first output)
    mglevel      = -1;
    reflevel     = -1;
    mc_grouptype = -1;
    map          = -1;
    component    = -1;
    
    // Say hello
    Waypoint ("Setting up the grid hierarchy");
    Output ("Carpet is running on %d processors", CCTK_nProcs(cctkGH));
    
    // Check arguments:
    // Only a specific number of dimensions is supported
    assert (cctkGH->cctk_dim == dim);
    // Not sure what to do with that
    assert (convLevel == 0);
    
    // Set up descriptions from user parameters
    setup_model_information (cctkGH);
    setup_multigrid_information (cctkGH);
    setup_refinement_information ();
    setup_map_information ();
    
    // Calculate domain extents for each map
    setup_domain_extents (cctkGH);
    
    // Set up grid hierarchy
    setup_grid_hierarchy (cctkGH);
    
    // Allocate space for group descriptions
    allocate_group_data (cctkGH);
    
    // Set times, origin, and spacings, and go to meta mode
    set_state (cctkGH);
    
    // Enable prolongating
    do_prolongate = true;
    do_warn_about_storage = false; // This is enabled later
    
    if (enable_all_storage) {
      enable_storage_for_all_groups (cctkGH);
    }
    
    Waypoint ("Done with setting up the grid hierarchy");
    
    return & carpetGH;
  }
  
  
  
  //////////////////////////////////////////////////////////////////////////////
  // Routines which perform actions ////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  
  
  
  void
  setup_model_information (cGH * const cctkGH)
  {
    DECLARE_CCTK_PARAMETERS;
    
    cctkGH->identity = strdup (model);
  }
  
  
  
  void
  setup_multigrid_information (cGH * const cctkGH)
  {
    DECLARE_CCTK_PARAMETERS;
    
    basemglevel = convergence_level;
    mglevels = num_convergence_levels;
    mgfact = convergence_factor;
    maxmglevelfact = ipow(mgfact, mglevels-1);
    cctkGH->cctk_convfac = mgfact;
  }
  
  
  
  void
  setup_refinement_information ()
  {
    DECLARE_CCTK_PARAMETERS;
    
    // Set maximum number of refinement levels
    maxreflevels = max_refinement_levels;
    
    // Set the per-level spatial refinement factors
    if (CCTK_EQUALS (space_refinement_factors, "")) {
      // Calculate them from the default refinement factor
      spacereffacts.resize (maxreflevels);
      for (int n=0; n<maxreflevels; ++n) {
        spacereffacts.at(n) = ivect (ipow (int(refinement_factor), n));
      }
    } else {
      // Read them from the parameter
      try {
        istringstream srf (space_refinement_factors);
        srf >> spacereffacts;
      } catch (input_error) {
        CCTK_WARN
          (0, "Could not parse parameter \"space_refinement_factors\"");
      }
    }
    // TODO: turn these into real error messages
    assert (int(spacereffacts.size()) >= maxreflevels);
    assert (all (spacereffacts.front() == 1));
    for (int n=1; n<maxreflevels; ++n) {
      assert (all (spacereffacts.at(n) >= spacereffacts.at(n-1)));
      assert (all (spacereffacts.at(n) % spacereffacts.at(n-1) == 0));
    }
    
    // Set the per-level temporal refinement factors
    if (CCTK_EQUALS (time_refinement_factors, "")) {
      // Calculate them from the default refinement factor
      timereffacts.resize (maxreflevels);
      for (int n=0; n<maxreflevels; ++n) {
        timereffacts.at(n) = ipow (int(refinement_factor), n);
      }
    } else {
      // Read them from the parameter
      try {
        istringstream trf (time_refinement_factors);
        trf >> timereffacts;
      } catch (input_error) {
        CCTK_WARN (0, "Could not parse parameter \"time_refinement_factors\"");
      }
    }
    // TODO: turn these into real error messages
    assert (int(timereffacts.size()) >= maxreflevels);
    assert (timereffacts.front() == 1);
    for (int n=1; n<maxreflevels; ++n) {
      assert (timereffacts.at(n) >= timereffacts.at(n-1));
      assert (timereffacts.at(n) % timereffacts.at(n-1) == 0);
    }
    
    // Calculate the maximum refinement factors
    maxtimereflevelfact = timereffacts.at (maxreflevels-1);
    maxspacereflevelfact = spacereffacts.at (maxreflevels-1);
  }
  
  
  
  void
  setup_map_information ()
  {
    DECLARE_CCTK_PARAMETERS;
    
    if (domain_from_multipatch) {
      assert (num_maps == 1);   // must be the default to avoid confusion
      assert (CCTK_IsFunctionAliased ("MultiPatch_GetSystemSpecification"));
      CCTK_INT maps1;
      int const ierr = MultiPatch_GetSystemSpecification (& maps1);
      assert (not ierr);
      maps = maps1;
    } else {
      maps = num_maps;
    }
    carpetGH.maps = maps;
  }
  
  
  
  void
  setup_domain_extents (cGH const * const cctkGH)
  {
    DECLARE_CCTK_PARAMETERS;
    
    for (int m=0; m<maps; ++m) {
      
      // Number of ghost zones
      ivect const lghosts = get_ghostzones();
      ivect const ughosts = lghosts;
      
      // Number of grid points
      ivect npoints = get_npoints ();
      
      // Boundary description
      get_boundary_specification (m);
      
      // Grid size
      rvect physical_min, physical_max;
      rvect base_spacing;
      get_domain_specification
        (m,
         npoints,
         physical_min, physical_max,
         base_spacing);
      
      if (max_refinement_levels > 1) {
        // Ensure that ReflectionSymmetry::avoid_origin = no
        ensure_ReflectionSymmetry_avoid_origin ();
      }
      
      // Adapt domain specification for convergence level
      rvect exterior_min, exterior_max;
      rvect spacing;
      adapt_domain_specification
        (m,
         physical_min, physical_max,
         base_spacing,
         exterior_min, exterior_max,
         spacing);
      
      // Calculate global number of grid points
      calculate_grid_points
        (m,
         lghosts, ughosts,
         exterior_min, exterior_max,
         spacing,
         npoints);
      
      // Allocate hierarchies
      allocate_grid_hierarchy (m, npoints);
      allocate_data_hierarchy (m, lghosts, ughosts);
      allocate_time_hierarchy (m);
      
    } // for m
  }
  
  
  
  void
  allocate_grid_hierarchy (int const m,
                           ivect const & npoints)
  {
    DECLARE_CCTK_PARAMETERS;
    
    // Centering
    centering refcentering;
    int reffactdenom;
    if (CCTK_EQUALS(refinement_centering, "vertex")) {
      refcentering = vertex_centered;
      reffactdenom = 1;
    } else if (CCTK_EQUALS(refinement_centering, "cell")) {
      refcentering = cell_centered;
      reffactdenom = 2;
    } else {
      assert (0);
    }
    
    // Base grid extent
    ivect const str (maxspacereflevelfact * reffactdenom);
    ivect const lb (0);
    ivect const ub ((npoints - 1) * str);
    ibbox const baseext (lb, ub, str);
    
    // Allocate grid hierarchy
    vhh.resize(maps);
    vhh.at(m) = new gh (spacereffacts, refcentering,
                        convergence_factor, vertex_centered,
                        baseext);
  }
  
  
  
  static ivect outer_buffer_width (ivect const & ghosts);
  
  void
  allocate_data_hierarchy (int const m,
                           ivect const & lghosts,
                           ivect const & ughosts)
  {
    DECLARE_CCTK_PARAMETERS;
    
    int const taper_factor = use_tapered_grids ? refinement_factor : 1;
    int const inner_buffer_width = use_outer_buffer_zones ? 0 : buffer_width;
    
    ivect const lbuffers =
      taper_factor * outer_buffer_width (lghosts) - lghosts;
    ivect const ubuffers =
      taper_factor * outer_buffer_width (ughosts) - ughosts;
    
    vdd.resize(maps);
    vdd.at(m) = new dh (* vhh.at(m),
                        lghosts, ughosts,
                        prolongation_order_space,
                        inner_buffer_width, lbuffers, ubuffers);
    
    if (max_refinement_levels > 1) {
      ensure_ghostzones (m, lghosts, ughosts);
    }
  }
  
  ivect
  outer_buffer_width (ivect const & ghosts)
  {
    DECLARE_CCTK_PARAMETERS;
    return (use_outer_buffer_zones
            ? ghosts * int(num_integrator_substeps) + int(buffer_width)
            : ghosts);
  }
  
  
  
  void
  allocate_time_hierarchy (int const m)
  {
    DECLARE_CCTK_PARAMETERS;
    
    vtt.resize (maps);
    vtt.at(m) = new th (* vhh.at(m),
                        timereffacts,
                        1.0);
  }
  
  
  
  void
  setup_grid_hierarchy (cGH const * const cctkGH)
  {
    DECLARE_CCTK_PARAMETERS;
    
    vector<vector<vector<region_t> > > regsss (maps);
    for (int m=0; m<maps; ++m) {
      set_base_extent (m, regsss.at(m));
    }
    
    vector<vector<vector<vector<region_t> > > > regssss;
    find_processor_decomposition (cctkGH, regsss, regssss);
    
    for (int m=0; m<maps; ++m) {
      
      // Check the regions
      CheckRegions (regssss.at(m));
      
#if 0
      // Do this later because CactusBase/IO might not yet be
      // initialised
      // Write grid structure to file
      OutputGridStructure (cctkGH, m, regssss.at(m));
#endif
      
      // Recompose grid hierarchy
      vhh.at(m)->regrid (regssss.at(m));
      int const rl = 0;
      vhh.at(m)->recompose (rl, false);
      
    } // for m
    
    CCTK_INFO ("Grid structure (grid points):");
    for (int ml=0; ml<mglevels; ++ml) {
      int const rl = 0;
      for (int m=0; m<maps; ++m) {
        for (int c=0; c<vhh.at(m)->components(rl); ++c) {
          ibbox const ext = vhh.at(m)->extent(ml,rl,c);
          ivect const lower = ext.lower();
          ivect const upper = ext.upper();
          int const convfact = ipow(mgfact, ml);
          assert (all(lower % maxspacereflevelfact == 0));
          assert (all(upper % maxspacereflevelfact == 0));
          assert (all(((upper - lower) / maxspacereflevelfact) % convfact == 0));
          cout << "   [" << ml << "][" << rl << "][" << m << "][" << c << "]"
               << "   exterior: "
               << "proc "
               << vhh.at(m)->processor(rl,c)
               << "   "
               << lower / maxspacereflevelfact
               << " : "
               << upper / maxspacereflevelfact
               << "   ("
               << (upper - lower) / maxspacereflevelfact / convfact + 1
               << ")" << endl;
        }
      }
    }
    
    // Assert that all maps have one refinement level
    reflevels = 1;
    for (int m=0; m<maps; ++m) {
      assert (vhh.at(m)->reflevels() == reflevels);
    }
  }
  
  
  
  void
  set_base_extent (int const m,
                   vector<vector<region_t> > & regss)
  {
    DECLARE_CCTK_PARAMETERS;
    
    // Create one refinement level
    int const rl = 0;
    regss.resize(1);
    vector<region_t> & regs = regss.at(rl);
    
    if (CCTK_EQUALS (base_extents, "")) {
      
      // Default: one grid component covering everything
      region_t reg;
      reg.extent = vhh.at(m)->baseextent;
      reg.outer_boundaries = b2vect (bvect (true));
      reg.refinement_boundaries = b2vect (bvect (false));
      reg.map = m;
      regs.push_back (reg);
      
    } else {
      
      // Read explicit grid components
      // TODO: invent something for the other convergence levels
      vector<ibbox> exts;
      istringstream ext_str (base_extents);
      try {
        ext_str >> exts;
      } catch (input_error) {
        CCTK_WARN (0, "Could not parse parameter \"base_extents\"");
      }
      CCTK_VInfo (CCTK_THORNSTRING,
                  "Using %d grid components", int(exts.size()));
      if (exts.size() == 0) {
        CCTK_WARN (0, "Cannot evolve with zero grid components");
      }
      
      vector<bbvect> obs;
      istringstream ob_str (base_outerbounds);
      try {
        ob_str >> obs;
      } catch (input_error) {
        CCTK_WARN (0, "Could not parse parameter \"base_outerbounds\"");
      }
      assert (obs.size() == exts.size());
      
      for (size_t n=0; n<exts.size(); ++n) {
        region_t reg;
        reg.extent = exts.at(n);
        reg.outer_boundaries = xpose(obs.at(n));
        reg.refinement_boundaries = b2vect(bvect(false));
        reg.map = m;
        regs.push_back (reg);
      }
      
    }
  }
  
  
  
  void
  allocate_group_data (cGH const * const cctkGH)
  {
    groupdata.resize (CCTK_NumGroups());
    arrdata.resize (CCTK_NumGroups());
    
    for (int group=0; group<CCTK_NumGroups(); ++group) {
      
      cGroup gdata;
      int const ierr = CCTK_GroupData (group, &gdata);
      assert (not ierr);
      
      // Check for compact, contiguous, and staggered groups
      ensure_group_options (group, gdata);
      
      switch (gdata.grouptype) {
        
        // Grid functions
      case CCTK_GF: {
        
        // All grid function groups must have the standard rank
        assert (gdata.dim == dim);
        
        // Set up one refinement level
        groupdata.at(group).activetimelevels.resize(mglevels);
        for (int ml=0; ml<mglevels; ++ml) {
          groupdata.at(group).activetimelevels.at(ml).resize(1);
        }
        
        // Grid function groups use the global grid descriptors
        arrdata.at(group).resize(maps);
        for (int m=0; m<maps; ++m) {
          arrdata.at(group).at(m).hh = vhh.at(m);
          arrdata.at(group).at(m).dd = vdd.at(m);
          arrdata.at(group).at(m).tt = vtt.at(m);
        }
        
        break;
      }
        
        // Grid scalars and grid arrays are treated in the same way
      case CCTK_SCALAR:
      case CCTK_ARRAY: {
        
        // Use only one refinement level for grid arrays
        groupdata.at(group).activetimelevels.resize(mglevels);
        for (int ml=0; ml<mglevels; ++ml) {
          groupdata.at(group).activetimelevels.at(ml).resize(1);
        }
        
        // Use only one map for grid arrays
        arrdata.at(group).resize(1);
        
        ivect sizes, lghosts, ughosts;
        get_group_size (group, gdata, sizes, lghosts, ughosts);
        
        // Adapt group sizes for convergence level
        ivect convpowers, convoffsets;
        adapt_group_size_mglevel (group, gdata, sizes, convpowers, convoffsets);
        // Adapt group sizes for disttype
        adapt_group_size_disttype
          (cctkGH, group, gdata, sizes, lghosts, ughosts);
        
        allocate_group_hierarchies (group, sizes, lghosts, ughosts);
        
        setup_group_grid_hierarchy
          (cctkGH, group, gdata, convpowers, convoffsets);
        
        break;
      } // case scalar or array
        
      default:
        assert (0);
      } // switch grouptype
      
      initialise_group_info (cctkGH, group, gdata);
      
    } // for groups
    
    output_group_statistics (cctkGH);
  }
  
  
  
  void
  allocate_group_hierarchies (int const group,
                              ivect const & sizes,
                              ivect const & lghosts,
                              ivect const & ughosts)
  {
    DECLARE_CCTK_PARAMETERS;
    
    // Calculate base extent
    ivect const lb(0);
    ivect const ub(sizes-1);
    ivect const str(1);
    ibbox const baseext(lb, ub, str);
    
    // One refinement level
    vector<int> grouptimereffacts(1);
    grouptimereffacts.at(0) = 1;
    vector<ivect> groupspacereffacts(1);
    groupspacereffacts.at(0) = ivect(1);
    
    // There is only one map
    int const m=0;
    
    arrdata.at(group).at(m).hh =
      new gh (groupspacereffacts, vertex_centered,
              convergence_factor, vertex_centered,
              baseext);
    
    arrdata.at (group).at(m).dd =
      new dh (*arrdata.at (group).at(m).hh, lghosts, ughosts, 0, 0, 0, 0);
    
    arrdata.at (group).at(m).tt =
      new th (*arrdata.at (group).at(m).hh, grouptimereffacts, 1.0);
  }
  
  
  
  void
  setup_group_grid_hierarchy (cGH const * const cctkGH,
                              int const group,
                              cGroup const & gdata,
                              ivect const & convpowers,
                              ivect const & convoffsets)
  {
    // Set refinement structure for scalars and arrays
    vector<region_t> regs(1);
    int const m=0;
    regs.at(0).extent = arrdata.at(group).at(m).hh->baseextent;
    regs.at(0).outer_boundaries = b2vect(true);
    regs.at(0).refinement_boundaries = b2vect(false);
    regs.at(m).map = m;
    
    // Split it into components, one for each processor
    switch (gdata.disttype) {
    case CCTK_DISTRIB_DEFAULT: {
      SplitRegions_Automatic (cctkGH, regs);
      break;
    }
    case CCTK_DISTRIB_CONSTANT: {
      int const d = gdata.dim==0 ? 0 : gdata.dim-1;
      SplitRegions_AlongDir (cctkGH, regs, d);
      break;
    }
    default:
      assert (0);
    }
    
    // Only one refinement level
    vector<vector<region_t> > regss(1);
    regss.at(0) = regs;
    
    // Create all multigrid levels
    vector<vector<vector<region_t> > > regsss (mglevels);
    ivect mgfact1;
    i2vect offset;
    for (int d=0; d<dim; ++d) {
      mgfact1[d] = ipow (mgfact, convpowers[d]);
      offset[0][d] = 0;
      offset[1][d] = convoffsets[d];
    }
    regsss.at(0) = regss;
    
    for (int ml=1; ml<mglevels; ++ml) {
      int const rl = 0;
      for (int c=0; c<int(regss.at(rl).size()); ++c) {
        // this base
        ivect const baselo = ivect(0);
        ivect const baselo1 = baselo;
        // next finer grid
        ivect const flo = regsss.at(ml-1).at(rl).at(c).extent.lower();
        ivect const fhi = regsss.at(ml-1).at(rl).at(c).extent.upper();
        ivect const fstr = regsss.at(ml-1).at(rl).at(c).extent.stride();
        // this grid
        ivect const str = fstr * mgfact1;
        ivect const lo = flo + regsss.at(ml-1).at(rl).at(c).outer_boundaries[0].ifthen
          (+ (offset[0] - mgfact1 * offset[0]) * fstr,
           ivect(0));
        ivect const hi = fhi + regsss.at(ml-1).at(rl).at(c).outer_boundaries[1].ifthen
          (- (offset[1] - mgfact1 * offset[1]) * fstr,
           ivect(0));
        ivect const lo1 = baselo1 + (lo - baselo1 + str - 1) / str * str;
        ivect const hi1 = lo1 + (hi - lo1) / str * str;
        regsss.at(ml).at(rl).at(c) = regsss.at(ml-1).at(rl).at(c);
        regsss.at(ml).at(rl).at(c).extent = ibbox(lo1, hi1, str);
      }
    }
    
    // Recompose for this map
    {
      char * const groupname = CCTK_GroupName (group);
      assert (groupname);
      Checkpoint ("Recomposing grid array group \"%s\"...", groupname);
      arrdata.at(group).at(0).hh->regrid (regsss);
      arrdata.at(group).at(0).hh->recompose (0, false);
      Checkpoint ("Done recomposing grid array group \"%s\".", groupname);
      free (groupname);
    }
  }
  
  
  
  void
  initialise_group_info (cGH const * const cctkGH,
                         int const group,
                         cGroup const & gdata)
  {
    DECLARE_CCTK_PARAMETERS;
    
    // Initialise group information
    groupdata.at(group).info.dim         = gdata.dim;
    groupdata.at(group).info.gsh         = new int [dim];
    groupdata.at(group).info.lsh         = new int [dim];
    groupdata.at(group).info.lbnd        = new int [dim];
    groupdata.at(group).info.ubnd        = new int [dim];
    groupdata.at(group).info.bbox        = new int [2*dim];
    groupdata.at(group).info.nghostzones = new int [dim];
    
    groupdata.at(group).transport_operator =
      get_transport_operator (cctkGH, group, gdata);
    
    groupdata.at(group).info.activetimelevels = 0;
    
    // Initialise group variables
    for (size_t m=0; m<arrdata.at(group).size(); ++m) {
      
      arrdata.at(group).at(m).data.resize(CCTK_NumVarsInGroupI(group));
      for (size_t var=0; var<arrdata.at(group).at(m).data.size(); ++var) {
        arrdata.at(group).at(m).data.at(var) = 0;
      }
      
    }
  }
  
  
  
  void
  set_state (cGH * const cctkGH)
  {
    // Allocate level times
    leveltimes.resize (mglevels);
    for (int ml=0; ml<mglevels; ++ml) {
      leveltimes.at(ml).resize (1);
    }
    
    // Allocate orgin and spacings
    origin_space.resize (maps);
    delta_space.resize (maps);
    for (int m=0; m<maps; ++m) {
      origin_space.at(m).resize (mglevels);
    }
    
    // Current state
    mglevelfact = 1;
    cctkGH->cctk_time = 0.0;
    cctkGH->cctk_delta_time = 1.0;
    for (int d=0; d<dim; ++d) {
      cctkGH->cctk_origin_space[d] = 0.0;
      cctkGH->cctk_delta_space[d] = 1.0;
    }
    
    // Set up things as if in local mode
    mglevel      = 0;
    reflevel     = 0;
    mc_grouptype = CCTK_GF;
    map          = 0;
    component    = 0;
    
    // Leave everything, so that everything is set up correctly
    leave_local_mode     (cctkGH);
    leave_singlemap_mode (cctkGH);
    leave_level_mode     (cctkGH);
    leave_global_mode    (cctkGH);
  }
  
  
  
  void
  enable_storage_for_all_groups (cGH const * const cctkGH)
  {
    BEGIN_MGLEVEL_LOOP(cctkGH) {
      BEGIN_REFLEVEL_LOOP(cctkGH) {
        
        for (int group=0; group<CCTK_NumGroups(); ++group) {
          char * const groupname = CCTK_GroupName(group);
          EnableGroupStorage (cctkGH, groupname);
          free (groupname);
        }
        
      } END_REFLEVEL_LOOP;
    } END_MGLEVEL_LOOP;
  }
  
  
  
  //////////////////////////////////////////////////////////////////////////////
  // Routines which do not change state ////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  
  
  
  ivect
  get_ghostzones ()
  {
    DECLARE_CCTK_PARAMETERS;
    
    // Decide which parameters to use
    ivect ghostzones;
    if (ghost_size == -1) {
      ghostzones = ivect (ghost_size_x, ghost_size_y, ghost_size_z);
    } else {
      ghostzones = ivect (ghost_size, ghost_size, ghost_size);
    }
    return ghostzones;
  }
  
  
  
  ivect
  get_npoints ()
  {
    DECLARE_CCTK_PARAMETERS;
    
    // Decide which parameters to use
    ivect npoints;
    if (global_nsize == -1) {
      npoints = ivect (global_nx, global_ny, global_nz);
    } else {
      npoints = ivect (global_nsize, global_nsize, global_nsize);
    }
    
    // Modify npoints for benchmarks
    if (constant_load_per_processor) {
      if (CCTK_EQUALS (processor_topology, "manual")) {
        // Enlarge the domain so that each processor has the specified
        // number of grid points, using the specified processor
        // topology
        assert (processor_topology_3d_x >= 1);
        assert (processor_topology_3d_y >= 1);
        assert (processor_topology_3d_z >= 1);
        ivect const topo =
          ivect (processor_topology_3d_x,
                 processor_topology_3d_y,
                 processor_topology_3d_z);
        npoints *= topo;
      } else if (CCTK_EQUALS (processor_topology, "automatic")) {
        // Enlarge the domain in a smart way so that each processor
        // has the specified number of grid points
        int const nprocs = dist::size();
        // Factorise the number of processors, placing the smallest
        // factors in the tail
        stack<int> factors;
        for (int procsleft = nprocs; procsleft > 1;) {
          for (int divisor = 2; divisor <= procsleft; ++divisor) {
            while (procsleft % divisor == 0) {
              factors.push (divisor);
              procsleft /= divisor;
            }
          }
        }
        // Distribute the factors greedily onto the directions,
        // starting with the largest factor, and preferring to enlarge
        // the x direction
        while (not factors.empty()) {
          int const mindir = minloc (npoints);
          assert (mindir>=0 and mindir<dim);
          int const factor = factors.top();
          factors.pop();
          npoints[mindir] *= factor;
        }
      } else {
        // TODO: handle processor_topology values "along-z" and "along-dir"
        CCTK_WARN (0, "Unsupported value of parameter processor_topology");
      }
    } // if constant_load_per_processor
    return npoints;
  }
  
  
  
  void
  get_boundary_specification (int const m)
  {
    DECLARE_CCTK_PARAMETERS;
    
    jjvect nboundaryzones, is_internal, is_staggered, shiftout;
    
    int const ierr =
      GetBoundarySpecification (2*dim,
                                &nboundaryzones[0][0],
                                &is_internal[0][0],
                                &is_staggered[0][0],
                                &shiftout[0][0]);
    assert (not ierr);
    
    ostringstream buf;
    buf << "CoordBase boundary specification for map " << m << ":" << endl
        << "   nboundaryzones: " << iivect(nboundaryzones) << endl
        << "   is_internal   : " << iivect(is_internal)    << endl
        << "   is_staggered  : " << iivect(is_staggered)   << endl
        << "   shiftout      : " << iivect(shiftout);
    Output (buf.str().c_str());
    
    if (max_refinement_levels > 1) {
      // Ensure that the boundary is not staggered
      for (int d=0; d<dim; ++d) {
        for (int f=0; f<2; ++f) {
          if (CCTK_EQUALS (refinement_centering, "vertex")) {
            if (is_staggered[d][f]) {
              CCTK_WARN (0, "The parameters CoordBase::boundary_staggered specify a staggered boundary.  Carpet does not support staggered boundaries when Carpet::max_refinement_levels > 1 with Carpet::centering = \"vertex\"");
            }
          } else if (CCTK_EQUALS (refinement_centering, "cell")) {
            if (not is_staggered[d][f]) {
              CCTK_WARN (0, "The parameters CoordBase::boundary_staggered specify a non-staggered boundary.  Carpet does not support non-staggered boundaries when Carpet::max_refinement_levels > 1 with Carpet::centering = \"cell\"");
            }
          } else {
            assert (0);
          }
        }
      }
    }
  }
  
  
  
  void
  get_domain_specification (int const m,
                            ivect const & npoints,
                            rvect & physical_min,
                            rvect & physical_max,
                            rvect & base_spacing)
  {
    DECLARE_CCTK_PARAMETERS;
    
    rvect interior_min, interior_max;
    rvect exterior_min, exterior_max;
    
    if (domain_from_multipatch) {
      assert (not domain_from_coordbase);
      
      // TODO: handle CartGrid3D: either add parameter
      // type=multipatch, and make it handle map numbers, or ignore it
      // altogether, maybe creating a new thorn
      
      assert (CCTK_IsFunctionAliased ("MultiPatch_GetDomainSpecification"));
      int const ierr = MultiPatch_GetDomainSpecification
        (m, dim,
         &physical_min[0], &physical_max[0],
         &interior_min[0], &interior_max[0],
         &exterior_min[0], &exterior_max[0], &base_spacing[0]);
      assert (not ierr);
      
    } else if (domain_from_coordbase) {
      
      // Ensure that CartGrid3D::type = "coordbase"
      ensure_CartGrid3D_type ();
      
      int const ierr = GetDomainSpecification
        (dim,
         &physical_min[0], &physical_max[0],
         &interior_min[0], &interior_max[0],
         &exterior_min[0], &exterior_max[0],
         &base_spacing[0]);
      assert (not ierr);
      
    } else {
      // Legacy code
      
      if (max_refinement_levels > 1) {
        // Ensure that CartGrid3D::avoid_origin = no
        ensure_CartGrid3D_avoid_origin ();
      } // if max_refinement_levels > 1
      
      ostringstream buf;
      buf << "Standard grid specification for map " << m << ":" << endl
          << "   number of grid points: " << npoints;
      Output (buf.str().c_str());
      
      // Reduce to physical domain
      // TODO: This is not the true domain specification.  However, it
      // is later written to the domainspec, and it is used by Carpet
      // for screen output.
      exterior_min = 0.0;
      exterior_max = rvect (npoints - 1);
      base_spacing = 1.0;
      int const ierr = ConvertFromExteriorBoundary
        (dim,
         &physical_min[0], &physical_max[0],
         &interior_min[0], &interior_max[0],
         &exterior_min[0], &exterior_max[0],
         &base_spacing[0]);
      assert (not ierr);
      
    } // if legacy domain specification
    
    ostringstream buf;
    buf << "CoordBase domain specification for map " << m << ":" << endl
        << "   physical extent: " << physical_min << " : " << physical_max 
        << "   (" << physical_max - physical_min << ")" << endl
        << "   interior extent: " << interior_min << " : " << interior_max 
        << "   (" << interior_max - interior_min << ")" << endl
        << "   exterior extent: " << exterior_min << " : " << exterior_max 
        << "   (" << exterior_max - exterior_min << ")" << endl
        << "   base_spacing   : " << base_spacing;
    Output (buf.str().c_str());
  }
  
  
  
  void
  adapt_domain_specification (int const m,
                              rvect const & physical_min,
                              rvect const & physical_max,
                              rvect const & base_spacing,
                              rvect & exterior_min,
                              rvect & exterior_max,
                              rvect & spacing)
  {
    DECLARE_CCTK_PARAMETERS;
    
    CCTK_REAL const baseconvfact =
      ipow (static_cast<CCTK_REAL> (convergence_factor), basemglevel);
    spacing = base_spacing * baseconvfact;
    
    rvect interior_min, interior_max;
    int const ierr =
      ConvertFromPhysicalBoundary
      (dim,
       &physical_min[0], &physical_max[0],
       &interior_min[0], &interior_max[0],
       &exterior_min[0], &exterior_max[0],
       &spacing[0]);
    assert (not ierr);
    
    ostringstream buf;
    buf << "Adapted domain specification for map " << m << ":" << endl
        << "   convergence factor: " << convergence_factor << endl
        << "   convergence level : " << basemglevel << endl
        << "   physical extent   : " << physical_min << " : " << physical_max
             << "   (" << physical_max - physical_min << ")" << endl
        << "   interior extent   : " << interior_min << " : " << interior_max
             << "   (" << interior_max - interior_min << ")" << endl
        << "   exterior extent   : " << exterior_min << " : " << exterior_max
             << "   (" << exterior_max - exterior_min << ")" << endl
        << "   spacing           : " << spacing;
    Output (buf.str().c_str());
  }
  
  
  
  void
  calculate_grid_points (int const m,
                         ivect const & lghosts,
                         ivect const & ughosts,
                         rvect const & exterior_min,
                         rvect const & exterior_max,
                         rvect const & spacing,
                         ivect & npoints)
  {
    DECLARE_CCTK_PARAMETERS;
    
    rvect const real_npoints
      = either (spacing,
                (exterior_max - exterior_min) / spacing + rvect(1),
                rvect(1));
    
    ostringstream buf;
    buf << "Base grid specification for map " << m << ":" << endl
        << "   number of grid points : " << real_npoints << endl
        << "   number of ghost points: " << lghosts;
    Output (buf.str().c_str());
    
    npoints = floor (real_npoints + static_cast<CCTK_REAL> (0.5));
    
    // Check domain size
    if (any (abs (rvect (npoints) - real_npoints) >
             static_cast<CCTK_REAL> (0.001)))
    {
      CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "The domain size for map %d scaled for convergence level %d with convergence factor %d is not integer",
                  m, int(basemglevel), int(convergence_factor));
    }
    
    // Sanity check
    assert (all (npoints <= INT_MAX));
    int max = INT_MAX;
    for (int d=0; d<dim; ++d) {
      assert (npoints[d] <= max);
      max /= npoints[d];
    }
    
    // Save domain specification
    domainspecs.resize(maps);
    domainspecs.at(m).exterior_min = exterior_min;
    domainspecs.at(m).exterior_max = exterior_max;
    domainspecs.at(m).npoints = npoints;
  }
  
  
  
  void
  find_processor_decomposition
  (cGH const * const cctkGH,
   vector<vector<vector<region_t> > >          & regsss,
   vector<vector<vector<vector<region_t> > > > & regssss)
  {
    DECLARE_CCTK_PARAMETERS;
    
    assert (regssss.empty());
    regssss.resize (maps);
    
    if (not regrid_in_level_mode) {
      // Distribute each map independently
      
      for (int m=0; m<maps; ++m) {
        int const rl=0;
        
        // Distribute onto the processors
        SplitRegions (cctkGH, regsss.at(m).at(rl));
        
        // Create all multigrid levels
        MakeMultigridBoxes (cctkGH, regsss.at(m), regssss.at(m));
      } // for m
      
    } else {
      // Distribute all maps at the same time
      
      int const rl=0;
      
      vector<vector<region_t> >  regss(maps);
      for (int m=0; m<maps; ++m) {
        regss.at(m) = regsss.at(m).at(rl);
      }
      
      SplitRegionsMaps (cctkGH, regss);
      
      for (int m=0; m<maps; ++m) {
        regsss.at(m).at(rl) = regss.at(m);
      }
      
      // Create all multigrid levels
      MakeMultigridBoxesMaps (cctkGH, regsss, regssss);
      
    } // if
  }
  
  
  
  void
  get_group_size (int const group,
                  cGroup const & gdata,
                  ivect & sizes,
                  ivect & lghosts,
                  ivect & ughosts)
  {
    // Default values
    sizes = 1;
    lghosts = 0;
    ughosts = 0;
    
    switch (gdata.grouptype) {
      
    case CCTK_SCALAR:
      // treat scalars as DIM=0, DISTRIB=const arrays
      assert (gdata.dim==0);
      assert (gdata.disttype == CCTK_DISTRIB_CONSTANT);
      break;
      
    case CCTK_ARRAY: {
      assert (gdata.dim>=1 or gdata.dim<=dim);
      CCTK_INT const * const * const sz  = CCTK_GroupSizesI (group);
      CCTK_INT const * const * const gsz = CCTK_GroupGhostsizesI (group);
      // Decode group sizes
      for (int d=0; d<gdata.dim; ++d) {
        if (sz) {
          sizes[d] = *sz[d];
        }
        if (gsz) {
          lghosts[d] = *gsz[d];
          ughosts[d] = *gsz[d];
        }
      }
      break;
    }
      
    default:
      assert (0);
    } // switch grouptype
  }
  
  
  
  void
  adapt_group_size_mglevel (int const group,
                            cGroup const & gdata,
                            ivect & sizes,
                            ivect & convpowers,
                            ivect & convoffsets)
  {
    DECLARE_CCTK_PARAMETERS;
    
    // Adapt array sizes for convergence level
    get_convergence_options (group, gdata, convpowers, convoffsets);
    
    ivect baseconvpowers = convpowers * int(basemglevel);
    rvect real_sizes =
      (((rvect (sizes)
         - rvect (convoffsets))
        / ipow (rvect(convergence_factor), baseconvpowers))
       + rvect (convoffsets));
    // Do not modify extra dimensions
    for (int d=gdata.dim; d<dim; ++d) {
      real_sizes[d] = sizes[d];
    }
    
    // Round group sizes
    sizes = floor (real_sizes + static_cast<CCTK_REAL> (0.5));
    
    if (any(sizes < 0)) {
      char * const groupname = CCTK_GroupName (group);
      CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "The shape of group \"%s\" scaled for convergence level %d with convergence factor %d is negative",
                  groupname, int(basemglevel), int(convergence_factor));
      free (groupname);
    }
    
    if (any(abs(rvect(sizes) - real_sizes) > static_cast<CCTK_REAL> (1.0e-8))) {
      char * const groupname = CCTK_GroupName(group);
      CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "The shape of group \"%s\" scaled for convergence level %d with convergence factor %d is not integer",
                  groupname, int(basemglevel), int(convergence_factor));
      free (groupname);
    }
  }
  
  
  
  void
  get_convergence_options (int const group,
                           cGroup const & gdata,
                           ivect & convpowers,
                           ivect & convoffsets)
  {
    if (gdata.tagstable >= 0) {
      {
        jvect convpowers1;
        int const status = Util_TableGetIntArray
          (gdata.tagstable,
           gdata.dim, &convpowers1[0], "convergence_power");
        if (status == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
          // use default: independent of convergence level
          convpowers = 0;
        } else if (status == 1) {
          // a scalar was given
          convpowers = convpowers1[0];
        } else if (status == gdata.dim) {
          convpowers = convpowers1;
        } else {
          char * const groupname = CCTK_GroupName (group);
          CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                      "The key \"convergence_power\" in the tags table of group \"%s\" is wrong",
                      groupname);
          free (groupname);
        }
        assert (all (convpowers >= 0));
      }
      
      {
        jvect convoffsets1;
        int const status = Util_TableGetIntArray
          (gdata.tagstable,
           gdata.dim, &convoffsets1[0], "convergence_offset");
        if (status == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
          // use default: offset is 0
          convoffsets = 0;
        } else if (status == 1) {
          // a scalar was given
          
        } else if (status == gdata.dim) {
          convoffsets = convoffsets1;
        } else {
          char * const groupname = CCTK_GroupName (group);
          CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                      "The key \"convergence_offset\" in the tags table of group \"%s\" is wrong",
                      groupname);
          free (groupname);
        }
      }
    }
  }
  
  
  
  void
  adapt_group_size_disttype (cGH const * const cctkGH,
                             int const group,
                             cGroup const & gdata,
                             ivect & sizes,
                             ivect const & lghosts,
                             ivect const & ughosts)
  {
    switch (gdata.disttype) {
      
    case CCTK_DISTRIB_DEFAULT: {
      // do nothing
      break;
    }
      
    case CCTK_DISTRIB_CONSTANT: {
      if (not all (lghosts == 0 and ughosts == 0)) {
        char * const groupname = CCTK_GroupName (group);
        CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "The group \"%s\" has DISTRIB=constant, but its "
                    "ghostsize is not 0",
                    groupname);
        free (groupname);
      }
      assert (all (lghosts == 0 and ughosts == 0));
      
      int const nprocs = CCTK_nProcs (cctkGH);
      
      // Find dimension which should be extended
      int const d = gdata.dim==0 ? 0 : gdata.dim-1;
      // Extend group sizes
      sizes[d] = ((sizes[d] - (lghosts[d] + ughosts[d])) * nprocs
                  + (lghosts[d] + ughosts[d]));
      assert (sizes[d] >= 0);
      break;
    }
      
    default:
      assert (0);
    } // switch disttype
  }
    
    
    
  void
  output_group_statistics (cGH const * const cctkGH)
  {
    int num_gf_groups = 0;
    int num_gf_vars = 0;
    vector<int> num_array_groups(dim+1), num_array_vars(dim+1);
    for (int d=0; d<=dim; ++d) {
      num_array_groups.at(d) = 0;
      num_array_vars.at(d) = 0;
    }
    
    for (int group=0; group<CCTK_NumGroups(); ++group) {
      
      cGroup gdata;
      int const ierr = CCTK_GroupData (group, & gdata);
      assert (not ierr);
      
      switch (gdata.grouptype) {
      case CCTK_GF:
        num_gf_groups += 1;
        num_gf_vars += gdata.numvars * gdata.numtimelevels;
        break;
      case CCTK_SCALAR:
      case CCTK_ARRAY:
        assert (gdata.dim<=dim);
        num_array_groups.at(gdata.dim) += 1;
        num_array_vars.at(gdata.dim) += gdata.numvars * gdata.numtimelevels;
        break;
      default:
        assert (0);
      }
    } // for group
    
    CCTK_INFO ("Group and variable statistics:");
    CCTK_VInfo (CCTK_THORNSTRING,
                "   There are %d grid functions in %d groups",
                num_gf_vars, num_gf_groups);
    CCTK_VInfo (CCTK_THORNSTRING,
                "   There are %d grid scalars in %d groups",
                num_array_vars.at(0), num_array_groups.at(0));
    for (int d=1; d<=dim; ++d) {
      CCTK_VInfo (CCTK_THORNSTRING,
                  "   There are %d %d-dimensional grid arrays in %d groups",
                  num_array_vars.at(d), d, num_array_groups.at(d));
    }
    CCTK_VInfo (CCTK_THORNSTRING,
                "   (The number of variables counts all time levels)");
  }
  
  
  
  operator_type
  get_transport_operator (cGH const * const cctkGH,
                          int const group,
                          cGroup const & gdata)
  {
    assert (group>=0 and group<CCTK_NumGroups());
    
    if (CCTK_GroupTypeI(group) != CCTK_GF) {
      // Ignore everything but true grid functions
      return op_error;
    }
    
    bool const can_transfer = can_transfer_variable_type (cctkGH, group, gdata);
    
    // Get prolongation method
    char prolong_string[1000];
    bool have_prolong_string = false;
    {
      int const prolong_length = Util_TableGetString
        (gdata.tagstable,
         sizeof prolong_string, prolong_string, "Prolongation");
      if (prolong_length >= 0) {
        have_prolong_string = true;
      } else if (prolong_length == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
        // do nothing
      } else {
        assert (0);
      }
    }
    
    // Get prolongation parameter name
    char prolong_param_string[1000];
    bool have_prolong_param_string = false;
    {
      int const prolong_param_length = Util_TableGetString
        (gdata.tagstable,
         sizeof prolong_param_string, prolong_param_string,
         "ProlongationParameter");
      if (prolong_param_length >= 0) {
        have_prolong_param_string = true;
      } else if (prolong_param_length == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
        // do nothing
      } else {
        assert (0);
      }
    }
    
    // Complain if both are given
    if (have_prolong_string and have_prolong_param_string) {
      char * const groupname = CCTK_GroupName (group);
      CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Group \"%s\" has both the tags \"Prolongation\" and \"ProlongationParameter\".  This is not possible.",
                  groupname);
      free (groupname);
    }
    
    // Map the parameter name
    if (have_prolong_param_string) {
      char * thorn;
      char * name;
      int const ierr2
        = CCTK_DecomposeName (prolong_param_string, &thorn, &name);
      if (ierr2 < 0) {
        char * const groupname = CCTK_GroupName (group);
        CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Group \"%s\" has the \"ProlongationParameter\" tag \"%s\".  This is not a valid parameter name.",
                    groupname, prolong_param_string);
        free (groupname);
      }
      int type;
      char const * const * const value
        = (static_cast<char const * const *>
           (CCTK_ParameterGet (name, thorn, &type)));
      if (not value or not *value) {
        char * const groupname = CCTK_GroupName (group);
        CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Group \"%s\" has the \"ProlongationParameter\" tag \"%s\".  This parameter does not exist.",
                    groupname, prolong_param_string);
        free (groupname);
      }
      if (type != PARAMETER_KEYWORD and type != PARAMETER_STRING) {
        char * const groupname = CCTK_GroupName (group);
        CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Group \"%s\" has the \"ProlongationParameter\" tag \"%s\".  This parameter has the wrong type; it must be either KEYWORD or STRING.",
                    groupname, prolong_param_string);
        free (groupname);
      }
      free (thorn);
      free (name);
      assert (strlen(*value) < sizeof prolong_string);
      strcpy (prolong_string, *value);
      have_prolong_string = true;
    }
    
    // Select a default, if necessary
    if (not have_prolong_string) {
      if (can_transfer) {
        // Use the default
        if (gdata.numtimelevels == 1) {
          // Only one time level:
          char * const groupname = CCTK_GroupName (group);
          CCTK_VWarn (2, __LINE__, __FILE__, CCTK_THORNSTRING,
                      "Group \"%s\" has only one time level; therefore it will not be prolongated or restricted.",
                      groupname);
          free (groupname);
          return op_none;
        } else {
          // Several time levels: use the default
          return op_Lagrange;
        }
      } else {
        if (gdata.grouptype == CCTK_GF) {
          char * const groupname = CCTK_GroupName (group);
          CCTK_VWarn (2, __LINE__, __FILE__, CCTK_THORNSTRING,
                      "Group \"%s\" has the variable type \"%s\" which cannot be prolongated or restricted.",
                      groupname, CCTK_VarTypeName(gdata.vartype));
          free (groupname);
          return op_none;
        } else {
          return op_error;
        }
      }
    }
    
    // Select the prolongation method
    assert (have_prolong_string);
    if (CCTK_Equals(prolong_string, "none")) {
      return op_none;
    } else if (CCTK_Equals(prolong_string, "copy")) {
      return op_copy;
    } else if (CCTK_Equals(prolong_string, "Lagrange")) {
      return op_Lagrange;
    } else if (CCTK_Equals(prolong_string, "ENO")) {
      return op_ENO;
    } else if (CCTK_Equals(prolong_string, "WENO")) {
      return op_WENO;
    } else {
      char * const groupname = CCTK_GroupName (group);
      CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Group \"%s\" has the unknown prolongation method \"%s\".",
                  groupname, prolong_string);
      free (groupname);
      return op_error;
    }
    return op_error;
  }
  
  
  
  bool
  can_transfer_variable_type (cGH const * const cctkGH,
                              int const group,
                              cGroup const & gdata)
  {
    // Find out which types correspond to the default types
#if CCTK_INTEGER_PRECISION_1
#  define CCTK_DEFAULT_INTEGER_TYPE CCTK_VARIABLE_INT1
#elif CCTK_INTEGER_PRECISION_2
#  define CCTK_DEFAULT_INTEGER_TYPE CCTK_VARIABLE_INT2
#elif CCTK_INTEGER_PRECISION_4
#  define CCTK_DEFAULT_INTEGER_TYPE CCTK_VARIABLE_INT4
#elif CCTK_INTEGER_PRECISION_8
#  define CCTK_DEFAULT_INTEGER_TYPE CCTK_VARIABLE_INT8
#else
#  error "Unsupported default integer type"
#endif
    
#if CCTK_REAL_PRECISION_4
#  define CCTK_DEFAULT_REAL_TYPE CCTK_VARIABLE_REAL4
#  define CCTK_DEFAULT_COMPLEX_TYPE CCTK_VARIABLE_COMPLEX8
#elif CCTK_REAL_PRECISION_8
#  define CCTK_DEFAULT_REAL_TYPE CCTK_VARIABLE_REAL8
#  define CCTK_DEFAULT_COMPLEX_TYPE CCTK_VARIABLE_COMPLEX16
#elif CCTK_REAL_PRECISION_16
#  define CCTK_DEFAULT_REAL_TYPE CCTK_VARIABLE_REAL16
#  define CCTK_DEFAULT_COMPLEX_TYPE CCTK_VARIABLE_COMPLEX32
#else
#  error "Unsupported default real type"
#endif
    
    int const type0 = gdata.vartype;
    int type1;
    
    switch (type0) {
    case CCTK_VARIABLE_INT:
      type1 = CCTK_DEFAULT_INTEGER_TYPE;
      break;
    case CCTK_VARIABLE_REAL:
      type1 = CCTK_DEFAULT_REAL_TYPE;
      break;
    case CCTK_VARIABLE_COMPLEX:
      type1 = CCTK_DEFAULT_COMPLEX_TYPE;
      break;
    default:
      type1 = type0;
    }
    switch (type1) {
      
#ifdef HAVE_CCTK_REAL8
    case CCTK_VARIABLE_REAL8:
      // This type is supported.
      return true;
#endif
      
#ifdef HAVE_CCTK_REAL4
    case CCTK_VARIABLE_REAL4:
#endif
#ifdef HAVE_CCTK_REAL16
    case CCTK_VARIABLE_REAL16:
#endif
#ifdef HAVE_CCTK_COMPLEX8
    case CCTK_VARIABLE_COMPLEX8:
#endif
#ifdef HAVE_CCTK_COMPLEX16
    case CCTK_VARIABLE_COMPLEX16:
#endif
#ifdef HAVE_CCTK_COMPLEX32
    case CCTK_VARIABLE_COMPLEX32:
#endif
      // This type is not supported, but could be.
      return false;
      
    case CCTK_VARIABLE_BYTE:
#ifdef HAVE_CCTK_INT1
    case CCTK_VARIABLE_INT1:
#endif
#ifdef HAVE_CCTK_INT2
    case CCTK_VARIABLE_INT2:
#endif
#ifdef HAVE_CCTK_INT4
    case CCTK_VARIABLE_INT4:
#endif
#ifdef HAVE_CCTK_INT8
    case CCTK_VARIABLE_INT8:
#endif
      // This type is not supported, and cannot be.
      return false;
      
    default:
      {
        CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Internal error: encountered variable type %d (%s) for group %d (%s)",
                    type1, CCTK_VarTypeName(type1),
                    group, CCTK_GroupName(group));
      }
    }
    
    // not reached
    abort ();
    return false;
  }



  //////////////////////////////////////////////////////////////////////////////
  // Parameter and consistency checking ////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  
  
  
  // Ensure that CartGrid3D::type = "coordbase"
  void
  ensure_CartGrid3D_type ()
  {
    if (CCTK_IsThornActive ("CartGrid3D")) {
      int type;
      void const * const ptr = CCTK_ParameterGet ("type", "CartGrid3D", & type);
      assert (ptr != 0);
      assert (type == PARAMETER_KEYWORD);
      char const * const coordtype
        = * static_cast<char const * const *> (ptr);
      if (not CCTK_EQUALS (coordtype, "coordbase")) {
        CCTK_WARN (0, "When Carpet::domain_from_coordbase = yes, and when thorn CartGrid3D is active, then you also have to set CartGrid3D::type = \"coordbase\"");
      }
    }
  }
  
  
  
  // Ensure that CartGrid3D::avoid_origin = no
  void
  ensure_CartGrid3D_avoid_origin ()
  {
    if (CCTK_IsThornActive ("CartGrid3D")) {
      int type;
      void const * ptr;
      
      ptr = CCTK_ParameterGet ("no_origin", "CartGrid3D", & type);
      assert (ptr != 0);
      assert (type == PARAMETER_BOOLEAN);
      CCTK_INT const no_origin = * static_cast<CCTK_INT const *> (ptr);
      
      ptr = CCTK_ParameterGet ("no_originx", "CartGrid3D", & type);
      assert (ptr != 0);
      assert (type == PARAMETER_BOOLEAN);
      CCTK_INT const no_originx = * static_cast<CCTK_INT const *> (ptr);
      
      ptr = CCTK_ParameterGet ("no_originy", "CartGrid3D", & type);
      assert (ptr != 0);
      assert (type == PARAMETER_BOOLEAN);
      CCTK_INT const no_originy = * static_cast<CCTK_INT const *> (ptr);
      
      ptr = CCTK_ParameterGet ("no_originz", "CartGrid3D", & type);
      assert (ptr != 0);
      assert (type == PARAMETER_BOOLEAN);
      CCTK_INT const no_originz = * static_cast<CCTK_INT const *> (ptr);
      
      ptr = CCTK_ParameterGet ("avoid_origin", "CartGrid3D", & type);
      assert (ptr != 0);
      assert (type == PARAMETER_BOOLEAN);
      CCTK_INT const avoid_origin = * static_cast<CCTK_INT const *> (ptr);
      
      ptr = CCTK_ParameterGet ("avoid_originx", "CartGrid3D", & type);
      assert (ptr != 0);
      assert (type == PARAMETER_BOOLEAN);
      CCTK_INT const avoid_originx = * static_cast<CCTK_INT const *> (ptr);
      
      ptr = CCTK_ParameterGet ("avoid_originy", "CartGrid3D", & type);
      assert (ptr != 0);
      assert (type == PARAMETER_BOOLEAN);
      CCTK_INT const avoid_originy = * static_cast<CCTK_INT const *> (ptr);
      
      ptr = CCTK_ParameterGet ("avoid_originz", "CartGrid3D", & type);
      assert (ptr != 0);
      assert (type == PARAMETER_BOOLEAN);
      CCTK_INT const avoid_originz = * static_cast<CCTK_INT const *> (ptr);
      
      ptr = CCTK_ParameterGet ("domain", "CartGrid3D", & type);
      assert (ptr != 0);
      assert (type == PARAMETER_KEYWORD);
      char const * const domain = * static_cast<char const * const *> (ptr);
      
      bvect stag;
      stag[0] = no_origin and no_originx and avoid_origin and avoid_originx;
      stag[1] = no_origin and no_originy and avoid_origin and avoid_originy;
      stag[2] = no_origin and no_originz and avoid_origin and avoid_originz;
      
      // TODO: Check only if there is actually a symmetry boundary
      if (not CCTK_EQUALS (domain, "full") and any(stag)) {
        CCTK_WARN (0, "When Carpet::domain_from_coordbase = no, when Carpet::max_refinement_levels > 1, and when thorn CartGrid3D provides symmetry boundaries, then you have to set CartGrid3D::avoid_origin = no");
      }
    }
  }
  
  
  
  // Ensure that ReflectionSymmetry::avoid_origin = no
  void
  ensure_ReflectionSymmetry_avoid_origin ()
  {
    if (CCTK_IsThornActive ("ReflectionSymmetry")) {
      int type;
      void const * ptr;
      
      ptr = CCTK_ParameterGet ("avoid_origin_x", "ReflectionSymmetry", & type);
      assert (ptr != 0);
      assert (type == PARAMETER_BOOLEAN);
      CCTK_INT const avoid_origin_x = * static_cast<CCTK_INT const *> (ptr);
      
      ptr = CCTK_ParameterGet ("avoid_origin_y", "ReflectionSymmetry", & type);
      assert (ptr != 0);
      assert (type == PARAMETER_BOOLEAN);
      CCTK_INT const avoid_origin_y = * static_cast<CCTK_INT const *> (ptr);
      
      ptr = CCTK_ParameterGet ("avoid_origin_z", "ReflectionSymmetry", & type);
      assert (ptr != 0);
      assert (type == PARAMETER_BOOLEAN);
      CCTK_INT const avoid_origin_z = * static_cast<CCTK_INT const *> (ptr);
      
      ptr = CCTK_ParameterGet ("reflection_x", "ReflectionSymmetry", & type);
      assert (ptr != 0);
      assert (type == PARAMETER_BOOLEAN);
      CCTK_INT const reflection_x = * static_cast<CCTK_INT const *> (ptr);
      
      ptr = CCTK_ParameterGet ("reflection_y", "ReflectionSymmetry", & type);
      assert (ptr != 0);
      assert (type == PARAMETER_BOOLEAN);
      CCTK_INT const reflection_y = * static_cast<CCTK_INT const *> (ptr);
      
      ptr = CCTK_ParameterGet ("reflection_z", "ReflectionSymmetry", & type);
      assert (ptr != 0);
      assert (type == PARAMETER_BOOLEAN);
      CCTK_INT const reflection_z = * static_cast<CCTK_INT const *> (ptr);
      
      if ((reflection_x and avoid_origin_x) or
          (reflection_y and avoid_origin_y) or
          (reflection_z and avoid_origin_z))
      {
        CCTK_WARN (0, "When Carpet::max_refinement_levels > 1, and when ReflectionSymmetry::symmetry_[xyz] = yes, then you also have to set ReflectionSymmetry::avoid_origin_[xyz] = no");
      }
    }
  }
  
  
  
  void
  ensure_ghostzones (int const m,
                     ivect const & lghosts, ivect const & ughosts)
  {
    DECLARE_CCTK_PARAMETERS;
    
    int const prolongation_stencil_size
      = vdd.at(m)->prolongation_stencil_size();
    int const min_nghosts
      = ((prolongation_stencil_size + refinement_factor - 1)
         / (refinement_factor - 1));
    if (any (min (lghosts, ughosts) < min_nghosts)) {
      CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "There are not enough ghost zones for the desired spatial prolongation order on map %d.  With Carpet::prolongation_order_space=%d, you need at least %d ghost zones.",
                  m, int(prolongation_order_space), min_nghosts);
    }
  }
  
  
  
  void
  ensure_group_options (int const group,
                        cGroup const & gdata)
  {
#ifdef CCTK_HAVE_COMPACT_GROUPS
    if (gdata.compact) {
      char * const groupname = CCTK_GroupName (group);
      CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "The group \"%s\" has COMPACT=1.  Compact groups are not yet supported",
                  groupname);
      free (groupname);
    }
#endif
    
#ifdef CCTK_HAVE_CONTIGUOUS_GROUPS
    if (gdata.contiguous) {
      char * const groupname = CCTK_GroupName (group);
      CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "The group \"%s\" has CONTIGUOUS=1.  Contiguous groups are not yet supported",
                  groupname);
      free (groupname);
    }
#endif
    
    // Staggered groups are not supported
    if (gdata.stagtype != 0) {
      char * const groupname = CCTK_GroupName (group);
      CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "The group \"%s\" is staggered.  Staggered groups are not yet supported",
                  groupname);
      free (groupname);
    }
  }
  
  
  
} // namespace Carpet
