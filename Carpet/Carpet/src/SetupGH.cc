#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <iostream>
#include <sstream>
#include <vector>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "util_ErrorCodes.h"
#include "util_Table.h"

#include "bbox.hh"
#include "defs.hh"
#include "dist.hh"
#include "ggf.hh"
#include "gh.hh"
#include "vect.hh"

#include "carpet.hh"

extern "C" {
  static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/SetupGH.cc,v 1.81 2004/06/26 11:34:02 schnetter Exp $";
  CCTK_FILEVERSION(Carpet_Carpet_SetupGH_cc);
}



namespace Carpet {
  
  using namespace std;
  
  
  
  static bool CanTransferVariableType (const cGH * const cgh, const int group)
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
    
    if (CCTK_NumVarsInGroupI(group) == 0) return true;
    
    const int var0 = CCTK_FirstVarIndexI(group);
    const int type0 = CCTK_VarTypeI(var0);
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
      
#ifdef CCTK_REAL8
    case CCTK_VARIABLE_REAL8:
      // This type is supported.
      return true;
#endif
      
#ifdef CCTK_REAL4
    case CCTK_VARIABLE_REAL4:
#endif
#ifdef CCTK_REAL16
    case CCTK_VARIABLE_REAL16:
#endif
#ifdef CCTK_REAL4 /* CCTK_COMPLEX8 */
    case CCTK_VARIABLE_COMPLEX8:
#endif
#ifdef CCTK_REAL8 /* CCTK_COMPLEX16 */
    case CCTK_VARIABLE_COMPLEX16:
#endif
#ifdef CCTK_REAL16 /* CCTK_COMPLEX32 */
    case CCTK_VARIABLE_COMPLEX32:
#endif
      // This type is not supported, but could be.
      return false;
      
    case CCTK_VARIABLE_BYTE:
#ifdef CCTK_INT1
    case CCTK_VARIABLE_INT1:
#endif
#ifdef CCTK_INT2
    case CCTK_VARIABLE_INT2:
#endif
#ifdef CCTK_INT4
    case CCTK_VARIABLE_INT4:
#endif
#ifdef CCTK_INT8
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
    return false;
  }
  
  
  
  static operator_type GetTransportOperator (const cGH * const cgh,
                                             const int group)
  {
    assert (group>=0 && group<CCTK_NumGroups());
    
    if (CCTK_GroupTypeI(group) != CCTK_GF) {
      // Ignore everything but true grid functions
      return op_error;
    }
    
    const bool can_transfer = CanTransferVariableType (cgh, group);
    
    cGroup gp;
    const int ierr = CCTK_GroupData (group, &gp);
    assert (!ierr);
    
    // Get prolongation method
    char prolong_string[1000];
    bool have_prolong_string = false;
    {
      const int prolong_length = Util_TableGetString
        (gp.tagstable, sizeof prolong_string, prolong_string, "Prolongation");
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
      const int prolong_param_length = Util_TableGetString
        (gp.tagstable, sizeof prolong_param_string, prolong_param_string,
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
    if (have_prolong_string && have_prolong_param_string) {
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
      int const ierr
        = CCTK_DecomposeName (prolong_param_string, &thorn, &name);
      if (ierr < 0) {
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
      if (! value || ! *value) {
        char * const groupname = CCTK_GroupName (group);
        CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Group \"%s\" has the \"ProlongationParameter\" tag \"%s\".  This parameter does not exist.",
                    groupname, prolong_param_string);
        free (groupname);
      }
      if (type != PARAMETER_KEYWORD && type != PARAMETER_STRING) {
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
    if (! have_prolong_string) {
      if (can_transfer) {
        // Use the default
        if (gp.numtimelevels == 1) {
          // Only one time level: do not prolongate
          char * const groupname = CCTK_GroupName (group);
          CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                      "Group \"%s\" has only one time level; therefore it will not be prolongated or restricted.",
                      groupname);
          free (groupname);
          return op_none;
        } else {
          // Several time levels: use the default
          return op_Lagrange;
        }
      } else {
        if (gp.grouptype == CCTK_GF) {
          char * const groupname = CCTK_GroupName (group);
          CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                      "Group \"%s\" has the variable type \"%s\" which cannot be prolongated or restricted.",
                      groupname, CCTK_VarTypeName(gp.vartype));
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
    } else if (CCTK_Equals(prolong_string, "Lagrange")) {
      return op_Lagrange;
    } else if (CCTK_Equals(prolong_string, "TVD")) {
      return op_TVD;
    } else if (CCTK_Equals(prolong_string, "ENO")) {
      return op_ENO;
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
  
  
  
  void* SetupGH (tFleshConfig* fc, int convLevel, cGH* cgh)
  {
    DECLARE_CCTK_PARAMETERS;
    
    int ierr;
    
    assert (cgh->cctk_dim == dim);
    
    // Not sure what to do with that
    assert (convLevel==0);
    
    dist::pseudoinit();
    
    // Initialise current position
    mglevel   = -1;
    reflevel  = -1;
    map       = -1;
    component = -1;
    
    
    
    Waypoint ("Setting up the grid hierarchy");
    
    // Processor information
    Output ("Carpet is running on %d processors", CCTK_nProcs(cgh));
    
    // Multigrid information
    basemglevel = convergence_level;
    mglevels = num_convergence_levels;
    mgfact = convergence_factor;
    maxmglevelfact = ipow(mgfact, mglevels-1);
    cgh->cctk_convfac = mgfact;
    
    // Refinement information
    maxreflevels = max_refinement_levels;
    reffact = refinement_factor;
    maxreflevelfact = ipow(reffact, maxreflevels-1);
    
    // Map information
    carpetGH.maps = maps = num_maps;;
    
    
    
    // Allocate space for groups
    groupdata.resize(CCTK_NumGroups());
    arrdata.resize(CCTK_NumGroups());
    
    vhh.resize(maps);
    vdd.resize(maps);
    vtt.resize(maps);
    
    // Loop over maps
    for (int m=0; m<maps; ++m) {
      
      // Get boundary description
      jjvect nboundaryzones, is_internal, is_staggered, shiftout;
      ierr = GetBoundarySpecification
        (2*dim, &nboundaryzones[0][0], &is_internal[0][0],
         &is_staggered[0][0], &shiftout[0][0]);
      assert (!ierr);
      
      {
        ostringstream buf;
        buf << "CoordBase boundary specification for map " << m << ":" << endl
            << "   nboundaryzones: " << nboundaryzones << endl
            << "   is_internal   : " << is_internal    << endl
            << "   is_staggered  : " << is_staggered   << endl
            << "   shiftout      : " << shiftout;
        Output (buf.str().c_str());
      }
      
      // Ghost zones
      ivect lghosts, ughosts;
      if (ghost_size == -1) {
        lghosts = ivect(ghost_size_x, ghost_size_y, ghost_size_z);
        ughosts = ivect(ghost_size_x, ghost_size_y, ghost_size_z);
      } else {
        lghosts = ivect(ghost_size, ghost_size, ghost_size);
        ughosts = ivect(ghost_size, ghost_size, ghost_size);
      }
      
      // Grid size
      rvect physical_min, physical_max;
      rvect interior_min, interior_max;
      rvect exterior_min, exterior_max;
      rvect base_spacing;
      
      if (domain_from_coordbase) {
        
        ierr = GetDomainSpecification
          (dim, &physical_min[0], &physical_max[0],
           &interior_min[0], &interior_max[0],
           &exterior_min[0], &exterior_max[0], &base_spacing[0]);
        assert (!ierr);
        
      } else {
        // Legacy code
        
        // specify global number of grid points
        ivect npoints;
        if (global_nsize == -1) {
          npoints = ivect(global_nx, global_ny, global_nz);
        } else {
          npoints = ivect(global_nsize, global_nsize, global_nsize);
        }
        ostringstream buf;
        buf << "Standard grid specification for map " << m << ":" << endl
            << "   number of grid points: " << npoints;
        Output (buf.str().c_str());
        
        // reduce to physical domain
        exterior_min = 0.0;
        exterior_max = rvect(npoints - 1);
        base_spacing = 1.0;
        ierr = ConvertFromExteriorBoundary
          (dim, &physical_min[0], &physical_max[0],
           &interior_min[0], &interior_max[0],
           &exterior_min[0], &exterior_max[0], &base_spacing[0]);
        assert (!ierr);
        
      }
      
      {
        ostringstream buf;
        buf << "CoordBase domain specification for map " << m << ":" << endl
            << "   physical extent: " << physical_min << " : " << physical_max << "   (" << physical_max - physical_min << ")" << endl
            << "   interior extent: " << interior_min << " : " << interior_max << "   (" << interior_max - interior_min << ")" << endl
            << "   exterior extent: " << exterior_min << " : " << exterior_max << "   (" << exterior_max - exterior_min << ")" << endl
            << "   base_spacing   : " << base_spacing;
        Output (buf.str().c_str());
      }
      
      // Adapt for convergence level
      rvect const spacing
        = base_spacing * pow (CCTK_REAL(convergence_factor), basemglevel);
      
      // Calculate global number of grid points
      ierr = ConvertFromPhysicalBoundary
        (dim, &physical_min[0], &physical_max[0],
         &interior_min[0], &interior_max[0],
         &exterior_min[0], &exterior_max[0], &spacing[0]);
      assert (!ierr);
      
      {
        ostringstream buf;
        buf << "Adapted domain specification for map " << m << ":" << endl
            << "   convergence factor: " << convergence_factor << endl
            << "   convergence level : " << basemglevel << endl
            << "   physical extent   : " << physical_min << " : " << physical_max << "   (" << physical_max - physical_min << ")" << endl
            << "   interior extent   : " << interior_min << " : " << interior_max << "   (" << interior_max - interior_min << ")" << endl
            << "   exterior extent   : " << exterior_min << " : " << exterior_max << "   (" << exterior_max - exterior_min << ")" << endl
            << "   spacing           : " << spacing;
        Output (buf.str().c_str());
      }
      
      rvect const real_npoints = (exterior_max - exterior_min) / spacing + 1;
      
      {
        ostringstream buf;
        buf << "Base grid specification for map " << m << ":" << endl
            << "   number of grid points : " << real_npoints << endl
            << "   number of ghost points: " << lghosts;
        Output (buf.str().c_str());
      }
      
      const ivect npoints = floor(real_npoints + 0.5);
      if (any(abs(rvect(npoints) - real_npoints) > 0.001)) {
        CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "The domain size for map %d scaled for convergence level %d with convergence factor %d is not integer",
                    m, basemglevel, convergence_factor);
      }
      
      // Sanity check
      // (if this fails, someone requested an insane amount of memory)
      assert (all(npoints <= INT_MAX));
      {
        int max = INT_MAX;
        for (int d=0; d<dim; ++d) {
          assert (npoints[d] <= max);
          max /= npoints[d];
        }
      }
      
      
      
      // Base grid extent
      const int stride = maxreflevelfact;
      const ivect str(stride);
      const ivect lb(0);
      const ivect ub((npoints - 1) * str);
      const ibbox baseext(lb, ub, str);
      
      // Allocate grid hierarchy
      vhh.at(m) = new gh<dim>(refinement_factor, vertex_centered,
                              convergence_factor, vertex_centered, baseext);
      
      // Allocate data hierarchy
      vdd.at(m) = new dh<dim>(*vhh.at(m), lghosts, ughosts,
                              prolongation_order_space, buffer_width);
      
      // Allocate time hierarchy
      vtt.at(m) = new th<dim>(*vhh.at(m), 1.0);
      
      if (max_refinement_levels > 1) {
        const int prolongation_stencil_size
          = vdd.at(m)->prolongation_stencil_size();
        const int min_nghosts
          = ((prolongation_stencil_size + refinement_factor - 1)
             / (refinement_factor-1));
        if (any(min(lghosts,ughosts) < min_nghosts)) {
          CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                      "There are not enough ghost zones for the desired spatial prolongation order on map %d.  With Carpet::prolongation_order_space=%d, you need at least %d ghost zones.",
                      m, prolongation_order_space, min_nghosts);
        }
        
      }
      
      
      
      // Set initial refinement structure
      
      vector<ibbox> bbs;
      vector<bbvect> obs;
      if (strcmp(base_extents, "") == 0) {
        
        // default: one grid component covering everything
        bbs.push_back (vhh.at(m)->baseextent);
        obs.push_back (bbvect(true));
        
      } else {
        
        // explicit grid components
        // TODO: invent something for the other convergence levels
        istringstream ext_str(base_extents);
        try {
          ext_str >> bbs;
        } catch (input_error) {
          CCTK_WARN (0, "Could not parse parameter \"base_extents\"");
        }
        CCTK_VInfo (CCTK_THORNSTRING, "Using %d grid patches", bbs.size());
        cout << "grid-patches-are " << bbs << endl;
        if (bbs.size()<=0) {
          CCTK_WARN (0, "Cannot evolve with 0 grid patches");
        }
        istringstream ob_str (base_outerbounds);
        try {
          ob_str >> obs;
        } catch (input_error) {
          CCTK_WARN (0, "Could not parse parameter \"base_outerbounds\"");
        }
        assert (obs.size() == bbs.size());
        
      }
      
      // Distribute onto the processors
      // (TODO: this should be done globally for all maps)
      vector<int> ps;
      SplitRegions (cgh, bbs, obs, ps);
      
      // Create all multigrid levels
      vector<vector<ibbox> > bbss;
      MakeMultigridBoxes (cgh, bbs, obs, bbss);
      
      // Only one refinement level
      vector<vector<vector<ibbox> > > bbsss(1);
      vector<vector<bbvect> > obss(1);
      vector<vector<int> > pss(1);
      bbsss.at(0) = bbss;
      obss.at(0) = obs;
      pss.at(0) = ps;
      
      // Check the regions
      CheckRegions (bbsss, obss, pss);
      
      // Write grid structure to file
      OutputGridStructure (cgh, m, bbsss, obss, pss);
      
      // Recompose grid hierarchy
      vhh.at(m)->recompose (bbsss, obss, pss);
      
      CCTK_INFO ("Grid structure (grid points):");
      const int rl = 0;
      for (int c=0; c<vhh.at(m)->components(rl); ++c) {
        for (int ml=0; ml<vhh.at(m)->mglevels(rl,c); ++ml) {
          const ivect lower = vhh.at(m)->extents.at(rl).at(c).at(ml).lower();
          const ivect upper = vhh.at(m)->extents.at(rl).at(c).at(ml).upper();
          const int convfact = ipow(mgfact, ml);
          assert (all(lower % maxreflevelfact == 0));
          assert (all(upper % maxreflevelfact == 0));
          assert (all(((upper - lower) / maxreflevelfact) % convfact == 0));
          cout << "   [" << ml << "][" << rl << "][" << m << "][" << c << "]"
               << "   exterior extent: " << lower / maxreflevelfact
               << " : " << upper / maxreflevelfact
               << "   (" << (upper - lower) / maxreflevelfact / convfact + 1
               << ")" << endl;
        }
      }
      
    } // loop over maps
    
    reflevels = 1;
    for (int m=0; m<maps; ++m) {
      assert (vhh.at(m)->reflevels() == reflevels);
    }
    
    
    
    // Allocate space for variables in group (but don't enable storage
    // yet)
    for (int group=0; group<CCTK_NumGroups(); ++group) {
      
      cGroup gp;
      ierr = CCTK_GroupData (group, &gp);
      assert (!ierr);
      
      switch (gp.grouptype) {
	
      case CCTK_GF: {
	assert (gp.dim == dim);
        arrdata.at(group).resize(maps);
        for (int m=0; m<maps; ++m) {
          arrdata.at(group).at(m).hh = vhh.at(m);
          arrdata.at(group).at(m).dd = vdd.at(m);
          arrdata.at(group).at(m).tt = vtt.at(m);
        }
	break;
      }
        
      case CCTK_SCALAR:
      case CCTK_ARRAY: {
        
        arrdata.at(group).resize(1);
        
	ivect sizes(1), ghostsizes(0);
	
	switch (gp.grouptype) {
	  
	case CCTK_SCALAR:
	  // treat scalars as DIM=0, DISTRIB=const arrays
	  assert (gp.dim==0);
	  assert (gp.disttype == CCTK_DISTRIB_CONSTANT);
	  break;
	  
	case CCTK_ARRAY: {
	  assert (gp.dim>=1 || gp.dim<=dim);
	  const CCTK_INT * const * const sz  = CCTK_GroupSizesI(group);
	  const CCTK_INT * const * const gsz = CCTK_GroupGhostsizesI(group);
	  for (int d=0; d<gp.dim; ++d) {
	    if (sz) sizes[d] = *sz[d];
	    if (gsz) ghostsizes[d] = *gsz[d];
	  }
	  break;
	}
          
	default:
	  assert (0);
	}
        
        ivect alghosts(0), aughosts(0);
        for (int d=0; d<gp.dim; ++d) {
          alghosts[d] = ghostsizes[d];
          aughosts[d] = ghostsizes[d];
        }
        
        
        
        // Adapt array sizes for convergence level
        jvect convpowers (0);
        jvect convoffsets (0);
        
        if (gp.tagstable >= 0) {
          int status;
          
          status = Util_TableGetIntArray
            (gp.tagstable, gp.dim, &convpowers[0], "convergence_power");
          if (status == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
            // keep default: independent of convergence level
          } else if (status == 1) {
            // a scalar was given
            convpowers = convpowers[0];
          } else if (status == gp.dim) {
            // do nothing
          } else {
            char * const groupname = CCTK_GroupName(group);
            CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                        "The key \"convergence_power\" in the tags table of group \"%s\" is wrong",
                        groupname);
            free (groupname);
          }
          assert (all (convpowers >= 0));
          
          status = Util_TableGetIntArray
            (gp.tagstable, gp.dim, &convoffsets[0], "convergence_offset");
          if (status == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
            // keep default: offset is 0
          } else if (status == 1) {
            // a scalar was given
            convoffsets = convoffsets[0];
          } else if (status == gp.dim) {
            // do nothing
          } else {
            char * const groupname = CCTK_GroupName(group);
            CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                        "The key \"convergence_offset\" in the tags table of group \"%s\" is wrong",
                        groupname);
            free (groupname);
          }
          
        } // if there is a group tags table
        
        rvect real_sizes
          = ((sizes - convoffsets)
             / pow(rvect(convergence_factor), convpowers * basemglevel)
             + convoffsets);
        for (int d=gp.dim; d<dim; ++d) {
          real_sizes[d] = sizes[d];
        }
        sizes = floor(real_sizes + 0.5);
        if (any(sizes < 0)) {
          char * const groupname = CCTK_GroupName(group);
          CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                      "The shape of group \"%s\" scaled for convergence level %d with convergence factor %d is negative",
                      groupname, basemglevel, convergence_factor);
          free (groupname);
        }
        if (any(abs(sizes - real_sizes) > 0.001)) {
          char * const groupname = CCTK_GroupName(group);
          CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                      "The shape of group \"%s\" scaled for convergence level %d with convergence factor %d is not integer",
                      groupname, basemglevel, convergence_factor);
          free (groupname);
        }
        
        
        
        assert (gp.disttype==CCTK_DISTRIB_CONSTANT
                || gp.disttype==CCTK_DISTRIB_DEFAULT);
        
        if (gp.disttype==CCTK_DISTRIB_CONSTANT) {
          if (! all (ghostsizes == 0)) {
            char * const groupname = CCTK_GroupName(group);
            CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                        "The group \"%s\" has DISTRIB=constant, but its ghostsize is not 0",
                        groupname);
            free (groupname);
          }
          assert (all (ghostsizes == 0));
          const int d = gp.dim==0 ? 0 : gp.dim-1;
          sizes[d] = (sizes[d] - 2*ghostsizes[d]) * CCTK_nProcs(cgh) + 2*ghostsizes[d];
          assert (sizes[d] >= 0);
        }
        
        
        
        const ivect alb(0);
        const ivect aub(sizes-1);
        const ivect astr(1);
        const ibbox abaseext(alb, aub, astr);
        
        assert (all(convpowers == convpowers[0]));
        const int amgfact1 = ipow(mgfact, convpowers[0]);
        
        arrdata.at(group).at(0).hh
          = new gh<dim>(refinement_factor, vertex_centered,
                        amgfact1, vertex_centered,
                        abaseext);
        
        arrdata.at(group).at(0).dd
          = new dh<dim>(*arrdata.at(group).at(0).hh, alghosts, aughosts, 0, 0);
        
        arrdata.at(group).at(0).tt
          = new th<dim>(*arrdata.at(group).at(0).hh, 1.0);
	
        
        
        // Set refinement structure for scalars and arrays
        vector<ibbox> bbs;
        vector<bbvect> obs;
        bbs.push_back (abaseext);
        obs.push_back (bbvect(true));
        
        // Split it into components, one for each processor
        vector<int> ps;
        if (gp.disttype==CCTK_DISTRIB_CONSTANT) {
          SplitRegions_AlongDir (cgh, bbs, obs, ps, gp.dim==0 ? 0 : gp.dim-1);
        } else {
          SplitRegions_Automatic (cgh, bbs, obs, ps);
        }
        
        // Create all multigrid levels
        vector<vector<ibbox> > bbss (bbs.size());
        ivect amgfact;
        iivect aoffset;
        for (int d=0; d<dim; ++d) {
          amgfact[d] = ipow(mgfact, convpowers[d]);
          aoffset[d][0] = 0;
          aoffset[d][1] = convoffsets[d];
        }
        for (size_t c=0; c<bbs.size(); ++c) {
          bbss.at(c).resize (mglevels);
          bbss.at(c).at(0) = bbs.at(c);
          for (int ml=1; ml<mglevels; ++ml) {
            // this base
            ivect const baselo = ivect(0);
            ivect const baselo1 = baselo;
            // next finer grid
            ivect const flo = bbss.at(c).at(ml-1).lower();
            ivect const fhi = bbss.at(c).at(ml-1).upper();
            ivect const fstr = bbss.at(c).at(ml-1).stride();
            // this grid
            ivect const str = fstr * amgfact;
            ivect const lo = flo + xpose(obs.at(c))[0].ifthen (   (xpose(aoffset)[0] - amgfact * xpose(aoffset)[0]) * fstr, ivect(0));
            ivect const hi = fhi + xpose(obs.at(c))[1].ifthen ( - (xpose(aoffset)[1] - amgfact * xpose(aoffset)[1]) * fstr, ivect(0));
            ivect const lo1 = baselo1 + (lo - baselo1 + str - 1) / str * str;
            ivect const hi1 = lo1 + (hi - lo1) / str * str;
            bbss.at(c).at(ml) = ibbox(lo1, hi1, str);
          }
        }
        
	// Only one refinement level
	vector<vector<vector<ibbox> > > bbsss(1);
	vector<vector<bbvect> > obss(1);
        vector<vector<int> > pss(1);
	bbsss.at(0) = bbss;
	obss.at(0) = obs;
        pss.at(0) = ps;
	
	// Recompose for this map
        char * const groupname = CCTK_GroupName (group);
        assert (groupname);
        Checkpoint ("Recomposing grid array group \"%s\"", groupname);
        free (groupname);
	arrdata.at(group).at(0).hh->recompose (bbsss, obss, pss);
	
	break;
      } // case of scalar or array
	
      default:
	assert (0);
      } // switch on group type
      
      // Initialise group information
      groupdata.at(group).info.dim         = gp.dim;
      groupdata.at(group).info.gsh         = new int [dim];
      groupdata.at(group).info.lsh         = new int [dim];
      groupdata.at(group).info.lbnd        = new int [dim];
      groupdata.at(group).info.ubnd        = new int [dim];
      groupdata.at(group).info.bbox        = new int [2*dim];
      groupdata.at(group).info.nghostzones = new int [dim];
      
      groupdata.at(group).transport_operator = GetTransportOperator (cgh, group);
      
      // Initialise group variables
      for (int m=0; m<(int)arrdata.at(group).size(); ++m) {
        
        arrdata.at(group).at(m).data.resize(CCTK_NumVarsInGroupI(group));
        for (int var=0; var<(int)arrdata.at(group).at(m).data.size(); ++var) {
          arrdata.at(group).at(m).data.at(var) = 0;
        }
        
      }
      
    } // for group
    
    
    
    // Allocate level times
    leveltimes.resize (mglevels);
    for (int ml=0; ml<mglevels; ++ml) {
      leveltimes.at(ml).resize (maxreflevels);
    }
    origin_space.resize (mglevels);
    
    // Enable prolongating
    do_prolongate = true;
    
    
    
    // Finish initialisation
    mglevelfact = 1;
    cgh->cctk_time = 0;
    cgh->cctk_delta_time = 1.0;
    for (int d=0; d<dim; ++d) {
      cgh->cctk_origin_space[d] = 0.0;
      cgh->cctk_delta_space[d] = 1.0;
    }
    
    mglevel   = 0;
    reflevel  = 0;
    map       = 0;
    component = 0;
    
    leave_local_mode     (cgh);
    leave_singlemap_mode (cgh);
    leave_level_mode     (cgh);
    leave_global_mode    (cgh);
    
    
    
    // Some statistics
    if (verbose || veryverbose) {
      
      int num_gf_groups = 0;
      int num_gf_vars = 0;
      vector<int> num_array_groups(dim+1), num_array_vars(dim+1);
      for (int d=0; d<=dim; ++d) {
        num_array_groups.at(d) = 0;
        num_array_vars.at(d) = 0;
      }
      
      for (int group=0; group<CCTK_NumGroups(); ++group) {
        
        cGroup data;
        const int ierr = CCTK_GroupData (group, &data);
        assert (!ierr);
        
        switch (data.grouptype) {
        case CCTK_GF:
          num_gf_groups += 1;
          num_gf_vars += data.numvars * data.numtimelevels;
          break;
        case CCTK_SCALAR:
        case CCTK_ARRAY:
          assert (data.dim<=dim);
          num_array_groups.at(data.dim) += 1;
          num_array_vars.at(data.dim) += data.numvars * data.numtimelevels;
          break;
        default:
          assert (0);
        }
      }
      CCTK_INFO ("Group and variable statistics:");
      CCTK_VInfo (CCTK_THORNSTRING,
                  "   There are %d grid functions in %d groups",
                  num_gf_vars, num_gf_groups);
      CCTK_VInfo (CCTK_THORNSTRING,
                  "   There are %d grid scalars in %d groups",
                  num_array_vars.at(0), num_array_groups.at(0));
      for (int dim=1; dim<=3; ++dim) {
        CCTK_VInfo (CCTK_THORNSTRING,
                    "   There are %d %d-dimensional grid arrays in %d groups",
                    num_array_vars.at(dim), dim, num_array_groups.at(dim));
      }
      CCTK_VInfo (CCTK_THORNSTRING,
                  "   (The number of variables counts all time levels)");
    }
    
    
    
    // Enable storage for all groups if desired
    if (true || enable_all_storage) {
      BEGIN_MGLEVEL_LOOP(cgh) {
        BEGIN_REFLEVEL_LOOP(cgh) {
	  for (int group=0; group<CCTK_NumGroups(); ++group) {
            char * const groupname = CCTK_GroupName(group);
	    EnableGroupStorage (cgh, groupname);
            free (groupname);
	  }
        } END_REFLEVEL_LOOP;
      } END_MGLEVEL_LOOP;
    }
    
    Waypoint ("Done with setting up the grid hierarchy");
    
    return &carpetGH;
  }
  
} // namespace Carpet
