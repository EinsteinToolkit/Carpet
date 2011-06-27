#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <cctk_Version.h>

#include <cassert>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include <hdf5.h>
#include <F5/F5F.h>
#include <F5/F5R.h>
#include <F5/F5X.h>
#include <F5/F5iterate.h>
#include <F5/F5uniform.h>

#include "CactusBase/IOUtil/src/ioutil_CheckpointRecovery.h"

#include <bbox.hh>
#include <vect.hh>

#include <carpet.hh>

#include "iof5.hh"



namespace CarpetIOF5 {
  
  using namespace std;
  using namespace Carpet;
  
  
  
  class output_iterator_t {
    // Can't be cGH const, since the mode loops change change its
    // entries
    cGH* const cctkGH;
    
    vector<bool> const output_var;    // whether to output this variable
    bool const output_everything;
    bool const is_multipatch;
    
    int group_type;             // CCTK_GF or CCTK_ARRAY
    int group_index;            // if group_type != CCTK_GF; else -1
    vector<int> vindices;       // variable indices to output
    
    string gridname;
    string chartname;
    string topologyname;
    string fragmentname;
    
    
    
    struct map_indices_t {
      int dim;
      ivect gsh;
      rvect origin, delta;
      rvect lower, upper;
      
      map_indices_t (cGH const* const cctkGH, int const gindex)
      {
        DECLARE_CCTK_ARGUMENTS;
        
        if (gindex == -1) {
          // grid function
          dim = ::dim;
          for (int d=0; d<(::dim); ++d) {
            gsh[d]    = cctk_gsh[d];
            origin[d] = CCTK_ORIGIN_SPACE(d);
            delta[d]  = CCTK_DELTA_SPACE(d);
          }
        } else {
          // grid array
          cGroupDynamicData dyndata;
          int const ierr = CCTK_GroupDynamicData (cctkGH, gindex, &dyndata);
          assert (not ierr);
          // HDF5 and F5 can't handle dim=0
          dim = max(dyndata.dim, 1);
          for (int d=0; d<dyndata.dim; ++d) {
            gsh[d]    = dyndata.gsh[d];
            origin[d] = 0.0;
            delta[d]  = 1.0;
          }
          for (int d=dyndata.dim; d<(::dim); ++d) {
            gsh[d]    = 1;
            origin[d] = 0.0;
            delta[d]  = 1.0;
          }
        }
        for (int d=0; d<(::dim); ++d) {
          lower[d] = origin[d];
          upper[d] = lower[d] + (gsh[d]-1) * delta[d];
        }
      }
    };
    
    struct component_indices_t: map_indices_t {
      // elements >=dim remain undefined
      ivect lbnd, lsh, lghosts, ughosts;
      ivect imin, imax, ioff, ilen;
      rvect clower, cupper;
      
      component_indices_t (cGH const* const cctkGH, int const gindex)
        : map_indices_t(cctkGH, gindex)
      {
        DECLARE_CCTK_ARGUMENTS;
        DECLARE_CCTK_PARAMETERS;
        
        if (gindex == -1) {
          // grid function
          for (int d=0; d<(::dim); ++d) {
            lbnd[d]    = cctk_lbnd[d];
            lsh[d]     = cctk_lsh[d];
            lghosts[d] = cctk_bbox[2*d  ] ? 0 : cctk_nghostzones[d];
            ughosts[d] = cctk_bbox[2*d+1] ? 0 : cctk_nghostzones[d];
          }
        } else {
          // grid array
          cGroupDynamicData dyndata;
          int const ierr = CCTK_GroupDynamicData (cctkGH, gindex, &dyndata);
          assert (not ierr);
          for (int d=0; d<dyndata.dim; ++d) {
            lbnd[d]    = dyndata.lbnd[d];
            lsh[d]     = dyndata.lsh[d];
            lghosts[d] = dyndata.bbox[2*d  ] ? 0 : dyndata.nghostzones[d];
            ughosts[d] = dyndata.bbox[2*d+1] ? 0 : dyndata.nghostzones[d];
          }
          for (int d=dyndata.dim; d<(::dim); ++d) {
            lbnd[d]    = 0;
            lsh[d]     = 1;
            lghosts[d] = 0;
            ughosts[d] = 0;
          }
        }
        for (int d=0; d<(::dim); ++d) {
          imin[d] = 0;
          imax[d] = lsh[d];
          if (not output_ghost_points) {
            int const overlap = min(ughosts[d], minimum_component_overlap);
            imin[d] += lghosts[d];
            imax[d] -= ughosts[d] - overlap;
            lghosts[d] = 0;
            ughosts[d] = overlap;
          }
          ioff[d] = lbnd[d] + imin[d];
          ilen[d] = imax[d] - imin[d];
          clower[d] = lower[d] + ioff[d] * delta[d];
          cupper[d] = clower[d] + (ilen[d]-1) * delta[d];
        }
      }
    };
    
    
    
  public:
    output_iterator_t (cGH* const cctkGH_,
                       vector<bool> const& output_var_,
                       bool const output_everything_)
      : cctkGH(cctkGH_),
        output_var(output_var_),
        output_everything(output_everything_),
        is_multipatch
        (CCTK_IsFunctionAliased("MultiPatch_GetSystemSpecification"))
    {
    }
    
    void iterate (hid_t const file)
    {
      // Iterate over the variables in groups, first all grid
      // functions, then all non-GF groups
      group_type = CCTK_GF;
      group_index = -1;
      vindices.clear();
      for (int vindex=0; vindex<CCTK_NumVars(); ++vindex) {
        if (output_var.at(vindex)) {
          int const gindex = CCTK_GroupIndexFromVarI(vindex);
          if (CCTK_GroupTypeI(gindex) == CCTK_GF) {
            vindices.push_back (vindex);
          }
        }
      }
      if (not vindices.empty()) {
        output_simulation (file);
      }
      
      group_type = CCTK_ARRAY;
      for (group_index=0; group_index<CCTK_NumGroups(); ++group_index) {
        if (CCTK_GroupTypeI(group_index) != CCTK_GF) {
          vindices.clear();
          int const first_vindex = CCTK_FirstVarIndexI (group_index);
          int const num_vars = CCTK_NumVarsInGroupI (group_index);
          for (int vindex=first_vindex; vindex<first_vindex+num_vars; ++vindex)
          {
            if (output_var.at(vindex)) {
              vindices.push_back (vindex);
            }
          }
          if (not vindices.empty()) {
            output_simulation (file);
          }
        }
      }
    }
    
    
    
  private:
    
    void output_simulation (hid_t const file)
    {
      DECLARE_CCTK_PARAMETERS;
      indent_t indent;
      bool error_flag = false;
      
      gridname = generate_gridname(cctkGH);
      cout << indent << "simulation=" << gridname << "\n";
      
      assert (is_global_mode());
      
      chartname = generate_chartname(cctkGH);
      
      int const max_rl = group_type == CCTK_GF ? reflevels : 1;
      for (int rl=0; rl<max_rl; ++rl) {
        // Continue only if this level exists at this iteration
        assert (maxtimereflevelfact % timereffacts.AT(rl) == 0);
        int const do_every =
          group_type == CCTK_GF ? maxtimereflevelfact / timereffacts.AT(rl) : 1;
        if (cctkGH->cctk_iteration % do_every == 0) {
          ENTER_LEVEL_MODE(cctkGH, rl) {
            DECLARE_CCTK_ARGUMENTS;
            
            assert (timelevel == 0);
            int const max_tl =
              output_everything or output_all_timelevels ?
              (group_type == CCTK_GF ?
               timelevels :
               CCTK_NumTimeLevelsI(group_index)) :
              1;
            for (timelevel=0; timelevel<max_tl; ++timelevel) {
              cctkGH->cctk_time = tt->get_time (mglevel, reflevel, timelevel);
              
              // Choose (arbitrarily) the root level as default
              // topology, for readers which don't understand AMR
              if (reflevel == 0) {
                ivect const reffact = spacereffacts.AT(reflevel);
                F5Path *const globalpath = F5Rcreate_vertex_refinement3D
                  (file, cctk_time, gridname.c_str(), &v2h(reffact)[0],
                   chartname.c_str());
                assert (globalpath);
                FAILWARN (F5Rlink_default_vertex_topology (globalpath,
                                                           &v2h(reffact)[0]));
                
                // Define iteration
                FAILWARN (F5Rset_timestep (globalpath, cctk_iteration));
                
                // Attach Cactus/Carpet metadata
                if (timelevel == 0) {
                  // hid_t const metadata_group = globalpath->Grid_hid;
                  ostringstream pathname;
                  pathname << FIBER_CONTENT_GRIDS << "/" << gridname;
                  hid_t group;
                  group = H5Gopen (globalpath->ContentsGroup_hid,
                                   pathname.str().c_str(),
                                   H5P_DEFAULT);
                  if (group < 0) {
                    group = H5Gcreate (globalpath->ContentsGroup_hid,
                                       pathname.str().c_str(),
                                       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                  }
                  assert (group >= 0);
                  write_metadata (cctkGH, group);
                  herr_t const herr = H5Gclose (group);
                  assert (not herr);
                }
                
                // Close topology
                F5close (globalpath);
                
              } // if reflevel==0
              
              output_reflevel (file);
              
            } // for timelevel
            timelevel = 0;
            cctkGH->cctk_time = tt->get_time (mglevel, reflevel, timelevel);
            
          } LEAVE_LEVEL_MODE;
        } // if do_every
      }   // for rl
      
    }
    
    
    
    void output_reflevel (hid_t const file)
    {
      DECLARE_CCTK_ARGUMENTS;
      indent_t indent;
      
      cout << indent << "reflevel=" << reflevel << "\n";
      
      assert (is_level_mode());
      
      ivect const reffact = spacereffacts.AT(reflevel);
      topologyname = generate_topologyname(cctkGH, group_index, reffact);
      
      // Define grid hierarchy
      map_indices_t const mi(cctkGH, group_index);
      int const indexdepth = vhh.at(0)->refcent == vertex_centered ? 0 : 1;
      F5Path *const path =
        F5Rcreate_coordinate_topology (file, &cctk_time,
                                       gridname.c_str(), chartname.c_str(),
                                       topologyname.c_str(),
                                       indexdepth,
                                       mi.dim, mi.dim, &v2h(reffact)[0]);
      assert (path);
      
      BEGIN_LOCAL_MAP_LOOP (cctkGH, CCTK_GF) {
        output_map (path);
      } END_LOCAL_MAP_LOOP;
      
      // Close topology
      F5close (path);
    }
    
    
    
    void output_map (F5Path *const path)
    {
      DECLARE_CCTK_ARGUMENTS;
      indent_t indent;
      bool error_flag = false;
      
      cout << indent << "map=" << Carpet::map << "\n";
      
      assert (is_singlemap_mode());
      
      if (not is_multipatch) {
        // Define level geometry
        map_indices_t const mi(cctkGH, group_index);
        F5_vec3_double_t const vlower = v2d(mi.lower);
        F5_vec3_double_t const vupper = v2d(mi.upper);
        F5_vec3_double_t const vdelta = v2d(mi.delta);
        static vector<ChartDomain_IDs*> charts;
        if (charts.size() < dim+1) {
          charts.resize(dim+1, NULL);
          charts.at(3) = F5B_standard_cartesian_chart3D();
        }
        if (not charts.at(mi.dim)) {
          assert (mi.dim != 0);
          char const* coordnames[] = {"x", "y", "z"};
          ostringstream chartnamebuf;
          chartnamebuf << "Cartesian " << mi.dim << "D";
          charts.at(mi.dim) =
            F5B_new_global_float_chart(coordnames,
                                       mi.dim, chartnamebuf.str().c_str(),
                                       F5_FORTRAN_ORDER);
          assert (charts.at(mi.dim));
        }
        hid_t const type = charts.at(mi.dim)->DoublePrecision.Point_hid_t;
        FAILWARN (F5Fwrite_linear (path, FIBER_HDF5_POSITIONS_STRING,
                                   mi.dim, &v2h(mi.gsh)[0],
                                   type,
                                   &vlower, &vdelta));
#warning "TODO: path and chart don't match"
        FAILWARN (F5Fset_range (path, &vlower, &vupper));
      }
      
      BEGIN_LOCAL_COMPONENT_LOOP (cctkGH, CCTK_GF) {
        output_component (path);
      } END_LOCAL_COMPONENT_LOOP;
    }
    
    
    
    void output_component (F5Path *const path)
    {
      DECLARE_CCTK_ARGUMENTS;
      DECLARE_CCTK_PARAMETERS;
      indent_t indent;
      bool error_flag = false;
      
      cout << indent
           << "component=" << component << " "
           << "(local_component=" << local_component << ")\n";
      
      assert (is_local_mode());
      
      fragmentname = generate_fragmentname(cctkGH, Carpet::map, component);
      
      // Define coordinates
      // TODO: also define and use is_cartesian for each map
      if (not is_multipatch and group_type == CCTK_GF) {
        // (This is redundant, since the level's overall bounding box
        // was already defined above, but it provides the individual
        // components' bounding boxes.)
        component_indices_t const ci(cctkGH, group_index);
        FAILWARN (F5Fwrite_linear_fraction (path, FIBER_HDF5_POSITIONS_STRING,
                                            ci.dim,
                                            &v2h(ci.gsh)[0], &v2h(ci.ilen)[0],
                                            F5T_COORD3_DOUBLE,
                                            &ci.clower, &ci.delta,
                                            &v2h(ci.ioff)[0],
                                            fragmentname.c_str()));
      } else {
        // Output coordinates
        output_variable (path, CCTK_VarIndex("grid::x"), true);
      }
      
      // Output variables
      for (vector<int>::const_iterator
             vi = vindices.begin(); vi != vindices.end(); ++vi)
      {
        output_variable (path, *vi);
      }
    }
    
    
    
    void output_variable (F5Path *const path, int const var,
                          bool const write_positions = false)
    {
      DECLARE_CCTK_ARGUMENTS;
      DECLARE_CCTK_PARAMETERS;
      indent_t indent;
      bool error_flag = false;
      
      int ierr;
      
      assert (var >= 0);
      if (write_positions) {
        cout << indent << "positions\n";
      } else {
        char *const fullname = CCTK_FullName(var);
        char *const groupname = CCTK_GroupNameFromVarI(var);
        cout << indent
             << "variable=" << fullname << " (group=" << groupname << ")\n";
        free (groupname);
        free (fullname);
      }
      
      int const group = CCTK_GroupIndexFromVarI(var);
      int const v0 = CCTK_FirstVarIndexI(group);
      
      cGroup groupdata;
      ierr = CCTK_GroupData (group, &groupdata);
      assert (not ierr);
      
      assert ((groupdata.grouptype == CCTK_GF) == (group_type == CCTK_GF));
      
      // Output distrib=constant variables only on process 0
      switch (groupdata.disttype) {
      case CCTK_DISTRIB_CONSTANT:
        if (CCTK_MyProc(cctkGH) != 0) return;
        break;
      case CCTK_DISTRIB_DEFAULT:
        // do nothing
        break;
      default:
        assert (0);
      }
      
      assert (groupdata.stagtype == 0);
      assert (groupdata.dim == dim);
      
#warning "TODO: Do not output symmetry zones (unless requested by the user)"
#warning "TODO: Do not output buffer zones (is that easily possible?)"
      
      int const will_cover_complete_domain = not is_multipatch and reflevel==0;
      
      cGroupDynamicData dyndata;
      ierr = CCTK_GroupDynamicData (cctkGH, group, &dyndata);
      assert (not ierr);
      
      // Only output active timelevels
      if (timelevel >= dyndata.activetimelevels) return;
      
      
      
      tensortype_t tensortype = tt_error;
      
      int const coordinates_group = CCTK_GroupIndex ("grid::coordinates");
      assert (coordinates_group >= 0);
      
      if (group == coordinates_group) {
        
        // Special case
        if (var >= v0 and var < v0+dim) {
          tensortype = tt_vector;
        } else if (var == v0+dim) {
          tensortype = tt_scalar;
        } else {
          assert (0);
        }
        
      } else {
        
        // TODO: use the tensor type instead
        switch (groupdata.numvars) {
        case 1:
          tensortype = tt_scalar; break;
        case 3:
          tensortype = tt_vector; break;
        case 6:
          tensortype = tt_symtensor; break;
        case 9:
          tensortype = tt_tensor; break;
        default:
          // fallback: output group as scalars
          tensortype = tt_scalar; break;
        }
        
      }
      
      if (tensortype != tt_scalar and var != v0) {
        // TODO: don't do this; instead, keep track of which groups
        // have been output
        indent_t indent2;
        cout << indent2 << "skipping output since it is not the first variable in its group\n";
        return;
      }
      
      
      
      hid_t type = -1;
      int num_comps = 0;
      string name;
      
      switch (tensortype) {
        
      case tt_scalar: {
        // Scalar field
        assert (not write_positions);
        switch (groupdata.vartype) {
        case CCTK_VARIABLE_INT:  type = H5T_NATIVE_INT;    break;
        case CCTK_VARIABLE_REAL: type = H5T_NATIVE_DOUBLE; break;
        default: assert(0);
        }
        num_comps = 1;
        name = generate_fieldname (cctkGH, var, tensortype);
        break;
      }
        
      case tt_vector: {
        // Vector field, or positions
        switch (groupdata.vartype) {
        case CCTK_VARIABLE_REAL:
          type = write_positions ? F5T_COORD3_DOUBLE : F5T_VEC3_DOUBLE;
          break;
        default: assert(0);
        }
        num_comps = dim;
        name =
          write_positions ?
          FIBER_HDF5_POSITIONS_STRING :
          generate_fieldname (cctkGH, var, tensortype);
        if (write_positions) {
          htri_t const exists =
            H5Lexists (path->Representation_hid, name.c_str(), H5P_DEFAULT);
          assert (exists >= 0);
          if (exists) {
            string const fragmentpath = name + "/" + fragmentname;
            htri_t const exists2 =
              H5Lexists (path->Representation_hid, fragmentpath.c_str(),
                         H5P_DEFAULT);
            assert (exists2 >= 0);
            if (exists2) {
              // Write positions only once
              indent_t indent2;
              cout << indent2 << "skipping output since the positions have already been output\n";
              return;
            }
          }
        }
        break;
      }
        
      case tt_symtensor: {
        // Symmetric tensor field
        assert (not write_positions);
        switch (groupdata.vartype) {
        case CCTK_VARIABLE_REAL: type = F5T_METRIC33_DOUBLE; break;
        default: assert(0);
        }
        num_comps = dim*(dim+1)/2;
        name = generate_fieldname (cctkGH, var, tensortype);
        break;
      }
        
      case tt_tensor: {
        // Non-symmetric tensor field
        assert (not write_positions);
        switch (groupdata.vartype) {
        case CCTK_VARIABLE_REAL: type = F5T_BIVEC3_DOUBLE; break;
        default: assert(0);
        }
        num_comps = dim*dim;
        name = generate_fieldname (cctkGH, var, tensortype);
        break;
      }
        
      default:
        assert (0);
      }
      
      // Write data
      assert (type >= 0);
      assert (num_comps > 0);
      assert (name.size() > 0);
      
      component_indices_t const ci(cctkGH, group_index);
      
      CCTK_REAL const* rdata[num_comps];
      vector<vector<CCTK_REAL> > idata(num_comps);
      void const* data[num_comps];
      for (int d=0; d<num_comps; ++d) {
        rdata[d] = (CCTK_REAL const*)CCTK_VarDataPtrI(cctkGH, timelevel, var+d);
        assert (rdata[d]);
        
        int const vartype = CCTK_VarTypeI(var);
        assert (vartype == CCTK_VARIABLE_REAL);
        idata[d].resize (prod(ci.ilen));
        for (int k=0; k<ci.ilen[2]; ++k) {
          for (int j=0; j<ci.ilen[1]; ++j) {
            for (int i=0; i<ci.ilen[0]; ++i) {
              int const isrc =
                ci.imin[0]+i + ci.lsh[0] *
                (ci.imin[1]+j + ci.lsh[1] * (ci.imin[2]+k));
              int const idst = i + ci.ilen[0] * (j + ci.ilen[1] * k);
              idata[d][idst] = rdata[d][isrc];
            }
          }
        }
        data[d] = &idata[d].front();
      }
      
      // Dataset properties
      hid_t const prop = H5Pcreate (H5P_DATASET_CREATE);
      assert (prop >= 0);
      // Enable compression if requested
      if (compression_level >= 0) {
        FAILWARN (H5Pset_chunk (prop, dim, &v2h(ci.ilen).reverse()[0]));
        FAILWARN (H5Pset_deflate (prop, compression_level));
      }
      // Enable checksums if requested
      if (use_checksums) {
        FAILWARN (H5Pset_chunk (prop, dim, &v2h(ci.ilen).reverse()[0]));
        FAILWARN (H5Pset_fletcher32 (prop));
      }
      
      if (num_comps == 1 and not separate_single_component_tensors) {
        // Write single-component tensors into non-separated fractions
        // for convenience (could also use a separated compound
        // instead)
        FAILWARN
          (F5Fwrite_fraction (path, name.c_str(),
                              dim,
                              is_multipatch ? NULL : &v2h(ci.gsh)[0],
                              &v2h(ci.ilen)[0],
                              type, type, data[0],
                              &v2h(ci.ioff)[0],
                              &v2h(ci.lghosts)[0], &v2h(ci.ughosts)[0],
                              fragmentname.c_str(), prop));
      } else {
        FAILWARN
          (F5FSwrite_fraction (path, name.c_str(),
                               dim,
                               is_multipatch ? NULL : &v2h(ci.gsh)[0],
                               &v2h(ci.ilen)[0],
                               type, type, data,
                               &v2h(ci.ioff)[0],
                               &v2h(ci.lghosts)[0], &v2h(ci.ughosts)[0],
                               fragmentname.c_str(), prop,
                               will_cover_complete_domain));
      }
      
      FAILWARN (H5Pclose (prop));
      
    }
    
  };                            // class output_iterator_t
  
  
  
  void write_metadata (cGH const* const cctkGH, hid_t const file)
  {
    DECLARE_CCTK_PARAMETERS;
    
    assert (cctkGH);
    assert (file >= 0);
    
    herr_t herr;
    
    
    CCTK_INFO ("Writing simulation metadata...");
    
    // NOTE: We could write the metadata into only one of the process
    // files (to save space), or write it into a separate metadata
    // file
    
    // Create metadata only once
    // TODO: update metadata at some point
    htri_t const exists = H5Lexists (file, metadata_group, H5P_DEFAULT);
    assert (exists >= 0);
    if (exists) return;
    
    // Create a group to hold all metadata
    hid_t const group =
      H5Gcreate (file, metadata_group, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    assert (group >= 0);
    
    // General information
    WriteAttribute (group, "Cactus version", CCTK_FullVersion());
    
    // Unique identifiers
    if (CCTK_IsFunctionAliased ("UniqueConfigID")) {
      WriteAttribute
        (group, "config id", (char const*) UniqueConfigID(cctkGH));
    }
    if (CCTK_IsFunctionAliased ("UniqueBuildID")) {
      WriteAttribute
        (group, "build id", (char const*) UniqueBuildID(cctkGH));
    }
    if (CCTK_IsFunctionAliased ("UniqueSimulationID")) {
      WriteAttribute
        (group, "simulation id", (char const*) UniqueSimulationID(cctkGH));
    }
    if (CCTK_IsFunctionAliased ("UniqueRunID")) {
      WriteAttribute
        (group, "run id", (char const*) UniqueRunID(cctkGH));
    }
    
    // Don't write this attribute; the number of files may change
    // after recovering, or after recombining files
#if 0
    // Number of I/O processes (i.e. the number of output files)
    int const nprocs = CCTK_nProcs(cctkGH);
    int const nioprocs = nprocs;
    WriteAttribute (group, "nioprocs", nioprocs);
#endif
    
    // Write parameters into a separate dataset (they may be too large
    // for an attribute)
    int const get_all = 1;
    char *const parameters = IOUtil_GetAllParameters (cctkGH, get_all);
    assert (parameters);
    WriteLargeAttribute (group, all_parameters, parameters);
    free (parameters);
    
    // Grid structure
    string const gs = serialise_grid_structure (cctkGH);
    WriteLargeAttribute (group, grid_structure, gs.c_str());
    
    herr = H5Gclose (group);
    assert (not herr);
  }
  
  
  
  void output (cGH const* const cctkGH,
               hid_t const file,
               vector<bool> const& output_var,
               bool const output_everything)
  {
    output_iterator_t
      iterator(const_cast<cGH*>(cctkGH), output_var, output_everything);
    iterator.iterate (file);
  }
  
} // end namespace CarpetIOF5
