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
    cGH *const cctkGH;
    
    bool const is_multipatch;
    
    string gridname;
    string chartname;
    string fragmentname;
    string topologyname;
    
  public:
    output_iterator_t (cGH *const cctkGH_)
      : cctkGH(cctkGH_),
        is_multipatch
        (CCTK_IsFunctionAliased("MultiPatch_GetSystemSpecification"))
    {
    }
    
    void iterate (hid_t const file)
    {
      int const rl = 0;
      ENTER_LEVEL_MODE(cctkGH, rl) {
        output_simulation (file);
      } LEAVE_LEVEL_MODE;
    }
    
    
    
  private:
    
    void output_simulation (hid_t const file)
    {
      assert (is_level_mode() and reflevel==0);
      DECLARE_CCTK_ARGUMENTS;
      indent_t indent;
      
      gridname = generate_gridname(cctkGH);
      chartname = generate_chartname(cctkGH);
      
      cout << indent << "simulation=" << gridname << "\n";
      
      ivect const reffact = spacereffacts.AT(reflevel);
      F5Path *const globalpath = F5Rcreate_vertex_refinement3D
        (file, cctk_time, gridname.c_str(), &v2h(reffact)[0],
         chartname.c_str());
      
      // Choose (arbitrarily) the root level as default topology, for
      // readers which don't understand AMR
      F5Rlink_default_vertex_topology (globalpath, &v2h(reffact)[0]);
      
      // Define iteration
      F5Rset_timestep (globalpath, cctk_iteration);
      
      // Attach Cactus/Carpet metadata
      write_metadata (cctkGH, globalpath->Grid_hid);
      
      // Close topology
      F5close (globalpath);
      
      // Iterate over all maps (on the root level)
      BEGIN_MAP_LOOP(cctkGH, CCTK_GF) {
        output_map (file);
      } END_MAP_LOOP;
    }
    
    
    
    void output_map (hid_t const file)
    {
      DECLARE_CCTK_ARGUMENTS;
      indent_t indent;
      
      assert (is_singlemap_mode() and reflevel==0);
      
      cout << indent << "map=" << Carpet::map << "\n";
      
      fragmentname = generate_fragmentname(cctkGH, Carpet::map);
      
      // Iterate over all refinement levels of this map
      int const m = Carpet::map;
      BEGIN_GLOBAL_MODE(cctkGH) {
        BEGIN_REFLEVEL_LOOP(cctkGH) {
          if (vhh.AT(m)->local_components(reflevel) > 0) {
            // Continue only if this process owns a component of this
            // map (on this level)
            ENTER_SINGLEMAP_MODE(cctkGH, m, CCTK_GF) {
              output_reflevel (file);
            } LEAVE_SINGLEMAP_MODE;
          }
        } END_REFLEVEL_LOOP;
      } END_GLOBAL_MODE;
    }
    
    
    
    void output_reflevel (hid_t const file)
    {
      DECLARE_CCTK_ARGUMENTS;
      indent_t indent;
      
      cout << indent << "reflevel=" << reflevel << "\n";
      
      ivect const reffact = spacereffacts.AT(reflevel);
      topologyname = generate_topologyname(cctkGH, Carpet::map, reffact);
      
      // Define grid hierarchy
      int const indexdepth = 0; // vertices
      F5Path *const path =
        F5Rcreate_coordinate_topology (file, &cctk_time,
                                       gridname.c_str(), chartname.c_str(),
                                       topologyname.c_str(),
                                       indexdepth,
                                       dim, dim, &v2h(reffact)[0]);
      
      // Determine level coordinates
      ivect gsh;
      rvect origin, delta;
      rvect lower, upper;
      for (int d=0; d<dim; ++d) {
        gsh[d]    = cctk_gsh[d];
        origin[d] = CCTK_ORIGIN_SPACE(d);
        delta[d]  = CCTK_DELTA_SPACE(d);
        lower[d]  = origin[d];
        upper[d]  = lower[d] + (gsh[d]-1) * delta[d];
      }
      
      if (not is_multipatch) {
        // Define level geometry
        F5_vec3_double_t const vlower = v2d(lower);
        F5_vec3_double_t const vupper = v2d(upper);
        F5_vec3_double_t const vdelta = v2d(delta);
        F5Fwrite_linear (path, FIBER_HDF5_POSITIONS_STRING,
                         dim, &v2h(gsh)[0],
                         F5T_COORD3_DOUBLE,
                         &vlower, &vdelta);
        F5Fset_range (path, &vlower, &vupper);
      }
      
      BEGIN_LOCAL_COMPONENT_LOOP (cctkGH, CCTK_GF) {
        output_component (path);
      } END_LOCAL_COMPONENT_LOOP;
      
      // Close topology
      F5close (path);
    }
  
  
  
    void output_component (F5Path *const path)
    {
      DECLARE_CCTK_ARGUMENTS;
      indent_t indent;
      
      cout << indent << "component=" << component << "\n";
      
      // Determine component coordinates
      ivect gsh;
      ivect lbnd, lsh, lghosts, ughosts;
      rvect origin, delta;
      rvect clower, cupper;
      for (int d=0; d<dim; ++d) {
        gsh[d]     = cctk_gsh[d];
        lbnd[d]    = cctk_lbnd[d];
        lsh[d]     = cctk_lsh[d];
        // F5 counts the total overlap, which is the sum of the ghost
        // zones on this and the adjacent component
        lghosts[d] = cctk_bbox[2*d  ] ? 0 : 2*cctk_nghostzones[d];
        ughosts[d] = cctk_bbox[2*d+1] ? 0 : 2*cctk_nghostzones[d];
        origin[d]  = CCTK_ORIGIN_SPACE(d);
        delta[d]   = CCTK_DELTA_SPACE(d);
        clower[d]  = origin[d] + cctk_lbnd[d] * delta[d];
        cupper[d]  = clower[d] + (lsh[d]-1) * delta[d];
      }
      
      // Define coordinates
      if (not is_multipatch) {
        // (This is redundant, since the level's overall bounding box
        // was already defined above, but it provides the individual
        // components' bounding boxes.)
        F5Fwrite_linear_fraction (path, FIBER_HDF5_POSITIONS_STRING,
                                  cctk_dim, &v2h(gsh)[0], &v2h(lsh)[0],
                                  F5T_COORD3_DOUBLE,
                                  &clower, &delta, &v2h(lbnd)[0],
                                  fragmentname.c_str());
      } else {
        output_variable (path, CCTK_VarIndex("grid::x"), true);
      }
      
      // Output some variables
      output_variable (path, CCTK_VarIndex("grid::r"));
      
      output_variable (path, CCTK_VarIndex("grid::x"));
      output_variable (path, CCTK_VarIndex("grid::y"));
      output_variable (path, CCTK_VarIndex("grid::z"));
    }
    
    
    
    void output_variable (F5Path *const path, int const var,
                          bool const write_positions = false)
    {
      DECLARE_CCTK_ARGUMENTS;
      DECLARE_CCTK_PARAMETERS;
      indent_t indent;
      
      int ierr;
      
      assert (var >= 0);
      {
        char *const fullname = CCTK_FullName(var);
        cout << indent << "variable=" << fullname << "\n";
        free (fullname);
      }
      
      int const group = CCTK_GroupIndexFromVarI(var);
      int const v0 = CCTK_FirstVarIndexI(group);
      
      cGroup groupdata;
      ierr = CCTK_GroupData (group, &groupdata);
      assert (not ierr);
      
      // TODO
      assert (groupdata.grouptype == CCTK_GF);
      assert (groupdata.vartype == CCTK_VARIABLE_REAL);
      assert (groupdata.disttype == CCTK_DISTRIB_DEFAULT);
      assert (groupdata.stagtype == 0);
      assert (groupdata.dim == cctk_dim);
      
      int const timelevel = 0;
      
      // Determine component coordinates
      int const dim = cctk_dim;
      ivect gsh;
      ivect lbnd, lsh, lghosts, ughosts;
      ivect imin, imax, ioff, ilen;
      ivect const izero(0);
      for (int d=0; d<dim; ++d) {
        gsh[d]  = cctk_gsh[d];
        lbnd[d] = cctk_lbnd[d];
        lsh[d]  = cctk_lsh[d];
        // F5 counts the total overlap, which is the sum of the ghost
        // zones on this and the adjacent component
        lghosts[d] = cctk_bbox[2*d  ] ? 0 : 2*cctk_nghostzones[d];
        ughosts[d] = cctk_bbox[2*d+1] ? 0 : 2*cctk_nghostzones[d];
        imin[d] = 0;
        imax[d] = lsh[d];
        if (not output_ghost_points) {
          imin[d] += lghosts[d] / 2;
          imax[d] -= ughosts[d] / 2;
          lghosts[d] = 0;
          ughosts[d] = 0;
        }
        ioff[d] = lbnd[d] + imin[d];
        ilen[d] = imax[d] - imin[d];
      }
      
#warning "TODO: Do not output symmetry zones (unless requested by the user)"
#warning "TODO: Do not output buffer zones (is that easily possible?)"
      
      
      
      enum tensortype_t {
        tt_scalar, tt_vector, tt_error
      };
      
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
        default:
          // fallback: output group as scalars
          tensortype = tt_scalar; break;
        }
        
      }
      
      if (tensortype != tt_scalar and var != v0) {
        // TODO: don't do this; instead, keep track of which groups
        // have been output
        char *const fullname = CCTK_FullName(var);
        indent_t indent;
        cout << indent << "skipping output since it is not the first variable in its group\n";
        free (fullname);
        return;
      }
      
      
      
      switch (tensortype) {
        
      case tt_scalar: {
        // Scalar field
        CCTK_REAL const *const rdata =
          (CCTK_REAL const*)CCTK_VarDataPtrI (cctkGH, timelevel, var);
        assert (rdata);
        
        int const vartype = CCTK_VarTypeI(var);
        assert (vartype == CCTK_VARIABLE_REAL);
        vector<CCTK_REAL> idata (prod(ilen));
        for (int k=0; k<ilen[2]; ++k) {
          for (int j=0; j<ilen[1]; ++j) {
            for (int i=0; i<ilen[0]; ++i) {
              int const isrc =
                CCTK_GFINDEX3D(cctkGH, imin[0]+i,imin[1]+j,imin[2]+k);
              int const idst = i + ilen[0] * (j + ilen[1] * k);
              idata[idst] = rdata[isrc];
            }
          }
        }
        void const *const data = &idata.front();
        
        char *const groupname = CCTK_GroupName(group);
        char *const fullname = CCTK_FullName(var);
        // Use the variable name instead of the group name if we may
        // output several variables per group
        assert (not write_positions);
        char const *const name = groupdata.numvars == 1 ? groupname : fullname;
#if 0
        F5Fwrite_fraction (path, name,
                           dim,
                           is_multipatch ? NULL : &v2h(gsh)[0], &v2h(lsh)[0],
                           H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE,
                           data,
                           &v2h(lbnd)[0], &v2h(lghosts)[0], &v2h(ughosts)[0],
                           fragmentname.c_str(), H5P_DEFAULT);
#endif
        F5Fwrite_fraction (path, name,
                           dim,
                           is_multipatch ? NULL : &v2h(gsh)[0], &v2h(ilen)[0],
                           H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE,
                           data,
                           &v2h(ioff)[0], &v2h(izero)[0], &v2h(izero)[0],
                           fragmentname.c_str(), H5P_DEFAULT);
        free (groupname);
        free (fullname);
        break;
      }
        
      case tt_vector: {
        // Vector field
        CCTK_REAL const* rdata[cctk_dim];
        vector<vector<CCTK_REAL> > idata(cctk_dim);
        void const* data[cctk_dim];
        for (int d=0; d<cctk_dim; ++d) {
          rdata[d] =
            (CCTK_REAL const*)CCTK_VarDataPtrI (cctkGH, timelevel, v0+d);
          assert (rdata[d]);
          
          int const vartype = CCTK_VarTypeI(var);
          assert (vartype == CCTK_VARIABLE_REAL);
          idata[d].resize (prod(ilen));
          for (int k=0; k<ilen[2]; ++k) {
            for (int j=0; j<ilen[1]; ++j) {
              for (int i=0; i<ilen[0]; ++i) {
                int const isrc =
                  CCTK_GFINDEX3D(cctkGH, imin[0]+i,imin[1]+j,imin[2]+k);
                int const idst = i + ilen[0] * (j + ilen[1] * k);
                idata[d][idst] = rdata[d][isrc];
              }
            }
          }
          data[d] = &idata[d].front();
          
        }
        char *const groupname = CCTK_GroupName(group);
        char const *const name =
          write_positions ? FIBER_HDF5_POSITIONS_STRING : groupname;
        hid_t const type =
          write_positions ? F5T_COORD3_DOUBLE : F5T_VEC3_DOUBLE;
        int const will_cover_complete_domain = 0;
#if 0
        F5FSwrite_fraction (path, name,
                            dim,
                            is_multipatch ? NULL : &v2h(gsh)[0], &v2h(lsh)[0],
                            type, type,
                            data,
                            &v2h(lbnd)[0], &v2h(lghosts)[0], &v2h(ughosts)[0],
                            fragmentname.c_str(), H5P_DEFAULT,
                            will_cover_complete_domain);
#endif
        F5FSwrite_fraction (path, name,
                            dim,
                            is_multipatch ? NULL : &v2h(gsh)[0], &v2h(ilen)[0],
                            type, type,
                            data,
                            &v2h(ioff)[0], &v2h(izero)[0], &v2h(izero)[0],
                            fragmentname.c_str(), H5P_DEFAULT,
                            will_cover_complete_domain);
        free (groupname);
        break;
      }
        
      default:
        assert (0);
      }
    }
    
  };                            // class output_iterator_t
  
  
  
  void write_metadata (cGH *const cctkGH, hid_t const file)
  {
    DECLARE_CCTK_PARAMETERS;
    
    assert (cctkGH);
    assert (file >= 0);
    
    herr_t herr;
    
    
    CCTK_INFO ("Writing simulation metadata...");
    
    // NOTE: We could write the metadata into only one of the process
    // files (to save space), or write it into a separate metadata
    // file
    
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
    
    // Number of I/O processes (i.e. the number of output files)
    int const nprocs = CCTK_nProcs(cctkGH);
    int const nioprocs = nprocs;
    WriteAttribute (group, "nioprocs", nioprocs);
    
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
  
  
  
  extern "C"
  void F5_Output (CCTK_ARGUMENTS)
  {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;
    
    herr_t herr;
    
    
    
    assert (is_global_mode());
    CCTK_VInfo (CCTK_THORNSTRING, "F5_Output: iteration=%d", cctk_iteration);
    
    
    
    // We don't know how to open multiple files yet
    assert (strcmp (file_content, "everything") == 0);
    
    // Open file
    static bool first_time = true;
    
    string const basename =
      generate_basename (cctkGH, CCTK_VarIndex("grid::r"));
    int const myproc = CCTK_MyProc(cctkGH);
    int const proc = myproc;
    string const name =
      create_filename (cctkGH, basename, proc, first_time);
    
    indent_t indent;
    cout << indent << "process=" << proc << "\n";
    
    bool const truncate_file = first_time and IO_TruncateOutputFiles(cctkGH);
    hid_t const file =
      truncate_file ?
      H5Fcreate (name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT) :
      H5Fopen (name.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    assert (file >= 0);
    first_time = false;
    
    // write_metadata (cctkGH, file);
    
    output_iterator_t iterator(cctkGH);
    iterator.iterate (file);
    
    // Close file
    herr = H5Fclose (file);
    assert (not herr);
  }
  
} // end namespace CarpetIOF5
