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
  // Can't be cGH const, since the mode loops change its entries
  cGH *const cctkGH;

  vector<bool> const output_var; // whether to output this variable
  bool const output_past_timelevels;
  bool const output_metadata;
  bool const is_multipatch;

  int group_type;       // CCTK_GF or CCTK_ARRAY
  int group_index;      // if group_type != CCTK_GF; else -1
  vector<int> vindices; // variable indices to output

  string gridname;
  string chartname;
  ivect slice_ipos;
  string topologyname;
  string fragmentname;

  struct map_indices_t {
    int dim;
    ivect gsh;
    rvect origin, delta;
    rvect lower, upper;

    map_indices_t(cGH const *const cctkGH, int const gindex) {
      DECLARE_CCTK_ARGUMENTS;

      if (gindex == -1) {
        // grid function
        dim = CarpetLib::dim;
        for (int d = 0; d < CarpetLib::dim; ++d) {
          gsh[d] = cctk_gsh[d];
          origin[d] = CCTK_ORIGIN_SPACE(d);
          delta[d] = CCTK_DELTA_SPACE(d);
        }
      } else {
        // grid array
        cGroupDynamicData dyndata;
        int const ierr = CCTK_GroupDynamicData(cctkGH, gindex, &dyndata);
        assert(not ierr);
        // HDF5 and F5 can't handle dim=0
        dim = max(dyndata.dim, 1);
        for (int d = 0; d < dyndata.dim; ++d) {
          gsh[d] = dyndata.gsh[d];
          origin[d] = 0.0;
          delta[d] = 1.0;
        }
        for (int d = dyndata.dim; d < CarpetLib::dim; ++d) {
          gsh[d] = 1;
          origin[d] = 0.0;
          delta[d] = 1.0;
        }
      }
      for (int d = 0; d < CarpetLib::dim; ++d) {
        lower[d] = origin[d];
        upper[d] = lower[d] + (gsh[d] - 1) * delta[d];
      }
    }
  };

  struct component_indices_t : map_indices_t {
    // elements >=dim remain undefined
    ivect lbnd, lsh, ash, lghosts, ughosts;
    ivect imin, imax, ioff, ilen;
    rvect clower, cupper;

    component_indices_t(cGH const *const cctkGH, int const gindex,
                        ivect const slice_ipos)
        : map_indices_t(cctkGH, gindex) {
      DECLARE_CCTK_ARGUMENTS;
      DECLARE_CCTK_PARAMETERS;

      if (gindex == -1) {
        // grid function
        for (int d = 0; d < CarpetLib::dim; ++d) {
          ash[d] = cctk_ash[d];
          lbnd[d] = cctk_lbnd[d];
          lsh[d] = cctk_lsh[d];
          lghosts[d] = cctk_bbox[2 * d] ? 0 : cctk_nghostzones[d];
          ughosts[d] = cctk_bbox[2 * d + 1] ? 0 : cctk_nghostzones[d];
        }
      } else {
        // grid array
        cGroupDynamicData dyndata;
        int const ierr = CCTK_GroupDynamicData(cctkGH, gindex, &dyndata);
        assert(not ierr);
        for (int d = 0; d < dyndata.dim; ++d) {
          ash[d] = dyndata.ash[d];
          lbnd[d] = dyndata.lbnd[d];
          lsh[d] = dyndata.lsh[d];
          lghosts[d] = dyndata.bbox[2 * d] ? 0 : dyndata.nghostzones[d];
          ughosts[d] = dyndata.bbox[2 * d + 1] ? 0 : dyndata.nghostzones[d];
        }
        for (int d = dyndata.dim; d < CarpetLib::dim; ++d) {
          ash[d] = 1;
          lbnd[d] = 0;
          lsh[d] = 1;
          lghosts[d] = 0;
          ughosts[d] = 0;
        }
      }
      for (int d = 0; d < CarpetLib::dim; ++d) {
        imin[d] = 0;
        imax[d] = lsh[d];
        if (slice_ipos[d] != -1) {
          imin[d] = max(imin[d], slice_ipos[d] - lbnd[d]);
          imax[d] = min(imax[d], slice_ipos[d] - lbnd[d] + 1);
        }
        if (not output_ghost_points and slice_ipos[d] < 0) {
          // TODO: Don't output ghosts on refinement boundaries;
          // only output ghosts for inter-process boundaries
          int const overlap = min(ughosts[d], int(minimum_component_overlap));
          imin[d] += lghosts[d];
          imax[d] -= ughosts[d] - overlap;
          lghosts[d] = 0;
          ughosts[d] = overlap;
        }
        ioff[d] = lbnd[d] + imin[d];
        ilen[d] = imax[d] - imin[d];
        clower[d] = lower[d] + ioff[d] * delta[d];
        cupper[d] = clower[d] + (ilen[d] - 1) * delta[d];
      }
    }
  };

public:
  output_iterator_t(cGH *const cctkGH_, vector<bool> const &output_var_,
                    bool const output_past_timelevels_,
                    bool const output_metadata_)
      : cctkGH(cctkGH_), output_var(output_var_),
        output_past_timelevels(output_past_timelevels_),
        output_metadata(output_metadata_),
        is_multipatch(
            CCTK_IsFunctionAliased("MultiPatch_GetSystemSpecification")) {}

  void iterate(hid_t const file) {
    // Iterate over the variables in groups, first all grid
    // functions, then all non-GF groups
    group_type = CCTK_GF;
    group_index = -1;
    vindices.clear();
    for (int vindex = 0; vindex < CCTK_NumVars(); ++vindex) {
      if (output_var.at(vindex)) {
        int const gindex = CCTK_GroupIndexFromVarI(vindex);
        if (CCTK_GroupTypeI(gindex) == CCTK_GF) {
          vindices.push_back(vindex);
        }
      }
    }
    if (not vindices.empty()) {
      output_simulation(file);
    }

    group_type = CCTK_ARRAY;
    for (group_index = 0; group_index < CCTK_NumGroups(); ++group_index) {
      if (CCTK_GroupTypeI(group_index) != CCTK_GF) {
        vindices.clear();
        int const first_vindex = CCTK_FirstVarIndexI(group_index);
        int const num_vars = CCTK_NumVarsInGroupI(group_index);
        for (int vindex = first_vindex; vindex < first_vindex + num_vars;
             ++vindex) {
          if (output_var.at(vindex)) {
            vindices.push_back(vindex);
          }
        }
        if (not vindices.empty()) {
          output_simulation(file);
        }
      }
    }
  }

private:
  void output_simulation(hid_t const file) {
    DECLARE_CCTK_PARAMETERS;
    indent_t indent;

    gridname = generate_gridname(cctkGH);
    cout << indent << "simulation=" << gridname << "\n";

    assert(is_global_mode());

    chartname = generate_chartname(cctkGH);
    cout << indent << "chartname=" << chartname << "\n";

    int const max_rl = group_type == CCTK_GF ? reflevels : 1;
    for (int rl = 0; rl < max_rl; ++rl) {
      // Continue only if this level exists at this iteration
      assert(maxtimereflevelfact % timereffacts.AT(rl) == 0);
      int const do_every =
          group_type == CCTK_GF ? maxtimereflevelfact / timereffacts.AT(rl) : 1;
      if (cctkGH->cctk_iteration % do_every == 0) {
        // TODO: don't switch modes
        ENTER_LEVEL_MODE(cctkGH, rl) {
          DECLARE_CCTK_ARGUMENTS;

          assert(timelevel == 0);
          int const max_tl =
              output_past_timelevels or output_all_timelevels
                  ? (group_type == CCTK_GF
                         ? timelevels
                         : CCTK_MaxActiveTimeLevelsVI(cctkGH, group_index))
                  : 1;
          for (timelevel = 0; timelevel < max_tl; ++timelevel) {
            cctkGH->cctk_time = tt->get_time(mglevel, reflevel, timelevel);

#if 0
              // Choose (arbitrarily) the root level as default
              // topology, for readers which don't understand AMR
              if (reflevel == 0) {
                ivect const reffact = spacereffacts.AT(reflevel);
                F5Path *const globalpath = F5Rcreate_vertex_refinement3D
                  (file, cctk_time, gridname.c_str(), &v2h(reffact)[0],
                   chartname.c_str());
                assert(globalpath);
                // TODO: Probably must not call this for cell-centred
                // AMR; this probably makes the call to
                // F5Rcreate_coordinate_topology below fail
                FAILWARN(F5Rlink_default_vertex_topology(globalpath,
                                                         &v2h(reffact)[0]));
                
                // Define iteration
                FAILWARN(F5Rset_timestep(globalpath, cctk_iteration));
                
                // Close topology
                F5close(globalpath);
              }
#endif

            {
              int const my_out_every =
                  out_every == -2 ? IO_out_every : out_every;
              int const my_out0D_every =
                  out0D_every == -2 ? my_out_every : out0D_every;
              int const my_out1D_every =
                  out1D_every == -2 ? my_out_every : out1D_every;
              int const my_out2D_every =
                  out2D_every == -2 ? my_out_every : out2D_every;
              int const my_out3D_every =
                  out3D_every == -2 ? my_out_every : out3D_every;

              rvect const x0(CCTK_ORIGIN_SPACE(0), CCTK_ORIGIN_SPACE(1),
                             CCTK_ORIGIN_SPACE(2));
              rvect const dx(CCTK_DELTA_SPACE(0), CCTK_DELTA_SPACE(1),
                             CCTK_DELTA_SPACE(2));
              rvect slice_pos;

              // 0D output
              if (group_index < 0 and my_out0D_every > 0 and
                  cctk_iteration % my_out0D_every == 0) {
                slice_ipos = rvect(0);
                output_reflevel(file);
              }

              // 1D output
              if (group_index < 0 and my_out1D_every > 0 and
                  cctk_iteration % my_out1D_every == 0) {
                // x
                slice_pos = rvect(0.0, out_xline_y, out_xline_z);
                slice_ipos = lrint((slice_pos - x0) / dx);
                slice_ipos[0] = -1;
                output_reflevel(file);

                // y
                slice_pos = rvect(out_yline_x, 0.0, out_xline_z);
                slice_ipos = lrint((slice_pos - x0) / dx);
                slice_ipos[1] = -1;
                output_reflevel(file);

                // z
                slice_pos = rvect(out_zline_x, out_zline_y, 0.0);
                slice_ipos = lrint((slice_pos - x0) / dx);
                slice_ipos[2] = -1;
                output_reflevel(file);
              }

              // 2D output
              if (group_index < 0 and my_out2D_every > 0 and
                  cctk_iteration % my_out2D_every == 0) {
                // xy
                slice_pos = rvect(0.0, 0.0, out_xyplane_z);
                slice_ipos = lrint((slice_pos - x0) / dx);
                slice_ipos[0] = slice_ipos[1] = -1;
                output_reflevel(file);

                // xz
                slice_pos = rvect(0.0, out_xzplane_y, 0.0);
                slice_ipos = lrint((slice_pos - x0) / dx);
                slice_ipos[0] = slice_ipos[2] = -1;
                output_reflevel(file);

                // yz
                slice_pos = rvect(out_yzplane_x, 0.0, 0.0);
                slice_ipos = lrint((slice_pos - x0) / dx);
                slice_ipos[1] = slice_ipos[2] = -1;
                output_reflevel(file);
              }

              // 3D output
              if (my_out3D_every > 0 and cctk_iteration % my_out3D_every == 0) {
                // xyz
                slice_ipos = rvect(-1);
                output_reflevel(file);
              }
            }

          } // for timelevel
          timelevel = 0;
          cctkGH->cctk_time = tt->get_time(mglevel, reflevel, timelevel);
        }
        LEAVE_LEVEL_MODE;
      } // if do_every
    }   // for rl
  }

  void output_reflevel(hid_t const file) {
    DECLARE_CCTK_ARGUMENTS;
    indent_t indent;
    bool error_flag = false;

    assert(is_level_mode());

    ivect const reffact = spacereffacts.AT(reflevel);
    topologyname =
        generate_topologyname(cctkGH, group_index, reffact, slice_ipos);
    cout << indent << "reflevel=" << reflevel << " "
         << "topologyname=" << topologyname << "\n";

    // Define grid hierarchy
    map_indices_t const mi(cctkGH, group_index);
    F5Path *path = NULL;
    F5Path *coordpath = NULL;
    bool close_coordpath = true;
    if (group_type == CCTK_GF) {
      int const indexdepth = vhh.at(0)->refcent == vertex_centered ? 0 : 1;
      if (indexdepth == 0) {
        path = F5Rcreate_coordinate_topology(
            file, &cctk_time, gridname.c_str(), chartname.c_str(),
            topologyname.c_str(), 0, mi.dim, mi.dim, &v2h(reffact)[0]);
      } else {
        char vertextopologyname[1000];
        TopologyName(vertextopologyname, sizeof vertextopologyname,
                     &v2h(reffact)[0], 0, dim);
        coordpath = F5Rcreate_coordinate_topology(
            file, &cctk_time, gridname.c_str(), chartname.c_str(),
            vertextopologyname, 0, mi.dim, mi.dim, &v2h(reffact)[0]);
        path = F5Rcreate_coordinate_topology(file, &cctk_time, gridname.c_str(),
                                             /*chartname.c_str(),*/
                                             vertextopologyname,
                                             topologyname.c_str(), indexdepth,
                                             mi.dim, mi.dim, &v2h(reffact)[0]);
        // TODO: how should these two topologies be linked?
        // assert(mi.dim == 3);
        // assert(all(reffact == 1));
        // path =
        //   F5Rcreate_hexaedrons_as_vertices_topology(file, &cctk_time,
        //                                             gridname.c_str());
      }
    } else {
      path = F5Rcreate_coordinate_topology(
          file, &cctk_time, gridname.c_str(), chartname.c_str(),
          topologyname.c_str(), 0, mi.dim, mi.dim, &v2h(reffact)[0]);
    }
    if (!coordpath) {
      coordpath = path;
      close_coordpath = false;
    }
    assert(coordpath);
    assert(path);

    // TODO: Attach the number of I/O processes to the topology. WB
    // sent email with suggestions.

    // Define default topology (once per grid)
    if (group_type == CCTK_GF and reflevel == 0 and timelevel == 0 and
        all(slice_ipos < 0)) {
      // TODO: Check that this happens at least once
      FAILWARN(F5Rlink_default_vertex_topology(path, &v2h(reffact)[0]));
    }

    // Define iteration
    if (reflevel == 0 and timelevel == 0) {
      FAILWARN(F5Rset_timestep(path, cctk_iteration));
    }

    // Attach Cactus/Carpet metadata (only once per output)
    if (output_metadata) {
      if (group_type == CCTK_GF and reflevel == 0 and timelevel == 0) {
        // hid_t const metadata_group = path->Grid_hid;
        ostringstream pathname;
        pathname << FIBER_CONTENT_GRIDS << "/" << gridname;
        hid_t group;
        group = H5Gopen(path->ContentsGroup_hid, pathname.str().c_str(),
                        H5P_DEFAULT);
        if (group < 0) {
          group = H5Gcreate(path->ContentsGroup_hid, pathname.str().c_str(),
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        }
        assert(group >= 0);
        write_metadata(cctkGH, group);
        herr_t const herr = H5Gclose(group);
        assert(not herr);
      }
    }

    BEGIN_LOCAL_MAP_LOOP(cctkGH, group_type) { output_map(path); }
    END_LOCAL_MAP_LOOP;

    // Close topologies
    if (close_coordpath)
      F5close(coordpath);
    F5close(path);
  }

  void output_map(F5Path *const path) {
    DECLARE_CCTK_ARGUMENTS;
    indent_t indent;
    bool error_flag = false;

    cout << indent << "map=" << Carpet::map << "\n";

    assert(is_singlemap_mode());

    if (group_type != CCTK_GF or not is_multipatch) {
      // Define level geometry
      map_indices_t const mi(cctkGH, group_index);
      F5_vec3_double_t const vlower = v2d(mi.lower);
      F5_vec3_double_t const vupper = v2d(mi.upper);
      F5_vec3_double_t const vdelta = v2d(mi.delta);
      static vector<ChartDomain_IDs *> charts;
      if (charts.size() < dim + 1) {
        charts.resize(dim + 1, NULL);
        charts.at(3) = F5B_standard_cartesian_chart3D();
      }
      if (not charts.at(mi.dim)) {
        assert(mi.dim != 0);
        char const *coordnames[] = {"x", "y", "z"};
        ostringstream chartnamebuf;
        chartnamebuf << "Cartesian " << mi.dim << "D";
        charts.at(mi.dim) = F5B_new_global_float_chart(
            coordnames, mi.dim, chartnamebuf.str().c_str(), F5_FORTRAN_ORDER);
        assert(charts.at(mi.dim));
      }
      // hid_t const type = charts.at(mi.dim)->DoublePrecision.Point_hid_t;
      hid_t const type = path->myChart->DoublePrecision.Point_hid_t;
      assert(type);
      FAILWARN(F5Fwrite_linear(path, FIBER_HDF5_POSITIONS_STRING, mi.dim,
                               &v2h(mi.gsh)[0], type, &vlower, &vdelta));
      // TODO: path and chart don't match
      FAILWARN(F5Fset_range(path, &vlower, &vupper));
    }

    BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, group_type) { output_component(path); }
    END_LOCAL_COMPONENT_LOOP;
  }

  void output_component(F5Path *const path) {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;
    indent_t indent;
    bool error_flag = false;

    assert(is_local_mode());

    fragmentname = generate_fragmentname(cctkGH, Carpet::map, component);
    cout << indent << "component=" << component << " "
         << "(local_component=" << local_component << ") "
         << "fragmentname=" << fragmentname << "\n";

    if (group_type == CCTK_GF) {
      // Define coordinates
      // TODO: also define and use is_cartesian for each map
      if (not is_multipatch) {
        // (This is redundant, since the level's overall bounding
        // box was already defined above, but it provides the
        // individual components' bounding boxes.)
        component_indices_t const ci(cctkGH, group_index, slice_ipos);
        FAILWARN(F5Fwrite_linear_fraction(
            path, FIBER_HDF5_POSITIONS_STRING, ci.dim, &v2h(ci.gsh)[0],
            &v2h(ci.ilen)[0], F5T_COORD3_DOUBLE, &ci.clower, &ci.delta,
            &v2h(ci.ioff)[0], fragmentname.c_str()));
      } else {
        // Output coordinates
        output_variable(path, CCTK_VarIndex("grid::x"), true);
      }
    }

    // Output variables
    for (vector<int>::const_iterator vi = vindices.begin();
         vi != vindices.end(); ++vi) {
      output_variable(path, *vi);
    }
  }

  void output_variable(F5Path *const path, int const var,
                       bool const write_positions = false) {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;
    indent_t indent;

    int ierr;

    assert(var >= 0);
    if (write_positions) {
      cout << indent << "positions\n";
    } else {
      char *const fullname = CCTK_FullName(var);
      char *const groupname = CCTK_GroupNameFromVarI(var);
      cout << indent << "variable=" << fullname << " (group=" << groupname
           << ")\n";
      free(groupname);
      free(fullname);
    }

    int const group = CCTK_GroupIndexFromVarI(var);
    int const v0 = CCTK_FirstVarIndexI(group);

    // Don't output groups without storage
    if (not CCTK_QueryGroupStorageI(cctkGH, group))
      return;

    cGroup groupdata;
    ierr = CCTK_GroupData(group, &groupdata);
    assert(not ierr);

    assert((groupdata.grouptype == CCTK_GF) == (group_type == CCTK_GF));

    // Output distrib=constant variables only on process 0
    switch (groupdata.disttype) {
    case CCTK_DISTRIB_CONSTANT:
      if (CCTK_MyProc(cctkGH) != 0)
        return;
      break;
    case CCTK_DISTRIB_DEFAULT:
      // do nothing
      break;
    default:
      assert(0);
    }

    // TODO: Do not output symmetry zones (unless requested by the
    // user)
    // TODO: Do not output buffer zones (is that easily possible?)

    cGroupDynamicData dyndata;
    ierr = CCTK_GroupDynamicData(cctkGH, group, &dyndata);
    assert(not ierr);

    // Only output active timelevels
    if (timelevel >= dyndata.activetimelevels)
      return;

    tensortype_t tensortype = tt_error;

    int const coordinates_group = CCTK_GroupIndex("grid::coordinates");
    assert(coordinates_group >= 0);

    if (group == coordinates_group) {

      // Special case
      if (var >= v0 and var < v0 + dim) {
        tensortype = tt_vector;
      } else if (var == v0 + dim) {
        tensortype = tt_scalar;
      } else {
        assert(0);
      }

    } else {

      // TODO: use the tensor type instead
      switch (groupdata.numvars) {
      case 1:
        tensortype = tt_scalar;
        break;
      case 3:
        tensortype = tt_vector;
        break;
      case 6:
        tensortype = tt_symtensor;
        break;
      case 9:
        tensortype = tt_tensor;
        break;
      default:
        // fallback: output group as scalars
        tensortype = tt_scalar;
        break;
      }
    }

    if (tensortype != tt_scalar and var != v0) {
      // TODO: don't do this; instead, keep track of which groups
      // have been output
      indent_t indent2;
      cout << indent2 << "skipping output since it is not the first variable "
                         "in its group\n";
      return;
    }

    hid_t type = -1;
    int num_comps = 0;
    string name;

    switch (tensortype) {

    case tt_scalar: {
      // Scalar field
      assert(not write_positions);
      switch (groupdata.vartype) {
      case CCTK_VARIABLE_INT:
        switch (sizeof(CCTK_INT)) {
        case 1:
          type = H5T_NATIVE_INT8;
          break;
        case 2:
          type = H5T_NATIVE_INT16;
          break;
        case 4:
          type = H5T_NATIVE_INT32;
          break;
        case 8:
          type = H5T_NATIVE_INT64;
          break;
        // case 16: type = H5T_NATIVE_INT128; break;
        default:
          CCTK_ERROR("Unsupported CCTK_INT type");
        }
        break;
      case CCTK_VARIABLE_REAL:
        switch (sizeof(CCTK_REAL)) {
        case 4:
          type = H5T_NATIVE_FLOAT;
          break;
        case 8:
          type = H5T_NATIVE_DOUBLE;
          break;
        case 16:
          type = H5T_NATIVE_LDOUBLE;
          break; // ???
        default:
          CCTK_ERROR("Unsupported CCTK_REAL type");
        }
        break;
      default:
        CCTK_ERROR("Unsupported type");
      }
      num_comps = 1;
      name = generate_fieldname(cctkGH, var, tensortype);
      break;
    }

    case tt_vector: {
      // Vector field, or positions
      switch (groupdata.vartype) {
      case CCTK_VARIABLE_REAL:
        type = write_positions ? F5T_COORD3_DOUBLE : F5T_VEC3_DOUBLE;
        break;
      default:
        assert(0);
      }
      num_comps = dim;
      // TODO: use generate_positionsname
      name = write_positions ? FIBER_HDF5_POSITIONS_STRING
                             : generate_fieldname(cctkGH, var, tensortype);
      if (write_positions) {
        htri_t const exists =
            H5Lexists(path->Representation_hid, name.c_str(), H5P_DEFAULT);
        assert(exists >= 0);
        if (exists) {
          string const fragmentpath = name + "/" + fragmentname;
          htri_t const exists2 = H5Lexists(path->Representation_hid,
                                           fragmentpath.c_str(), H5P_DEFAULT);
          assert(exists2 >= 0);
          if (exists2) {
            // Write positions only once
            indent_t indent2;
            cout << indent2 << "skipping output since the positions have "
                               "already been output\n";
            return;
          }
        }
      }
      break;
    }

    case tt_symtensor: {
      // Symmetric tensor field
      assert(not write_positions);
      switch (groupdata.vartype) {
      case CCTK_VARIABLE_REAL:
        type = F5T_METRIC33_DOUBLE;
        break;
      default:
        assert(0);
      }
      num_comps = dim * (dim + 1) / 2;
      name = generate_fieldname(cctkGH, var, tensortype);
      break;
    }

    case tt_tensor: {
      // Non-symmetric tensor field
      assert(not write_positions);
      switch (groupdata.vartype) {
      case CCTK_VARIABLE_REAL:
        type = F5T_BIVEC3_DOUBLE;
        break;
      default:
        assert(0);
      }
      num_comps = dim * dim;
      name = generate_fieldname(cctkGH, var, tensortype);
      break;
    }

    default:
      assert(0);
    }

    cout << indent << "fieldname=" << name << "\n";

    output_hyperslab(path, var, type, num_comps, name);
  }

  void output_hyperslab(F5Path *const path, int const var, hid_t const type,
                        int const num_comps, string const name) {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;
    indent_t indent;
    bool error_flag = false;

    cout << indent << "hyperslab=" << (slice_ipos < ivect(0)) << "\n";

    // Write data
    assert(type >= 0);       // HDF5 datatype
    assert(num_comps > 0);   // field components
    assert(name.size() > 0); // field name

    // Ensure that the data types match
    int const vartype = CCTK_VarTypeI(var);
    assert(num_comps * CCTK_VarTypeSize(vartype) == (int)H5Tget_size(type));

    component_indices_t const ci(cctkGH, group_index, slice_ipos);
    // Do not output empty datasets, or slices that do not intersect
    // us
    if (any(ci.ilen <= 0))
      return;

    int const will_cover_complete_domain =
        (group_type != CCTK_GF or not is_multipatch) and reflevel == 0;

    void const *data[num_comps];
    switch (vartype) {

    case CCTK_VARIABLE_INT: {
      for (int d = 0; d < num_comps; ++d) {
        CCTK_INT const *const varptr =
            (CCTK_INT const *)CCTK_VarDataPtrI(cctkGH, timelevel, var + d);
        assert(varptr);
        CCTK_INT *const rdata = new CCTK_INT[prod(ci.ilen)];
        for (int k = 0; k < ci.ilen[2]; ++k) {
          for (int j = 0; j < ci.ilen[1]; ++j) {
            for (int i = 0; i < ci.ilen[0]; ++i) {
              int const isrc =
                  (ci.imin[0] + i +
                   ci.ash[0] * (ci.imin[1] + j + ci.ash[1] * (ci.imin[2] + k)));
              int const idst = i + ci.ilen[0] * (j + ci.ilen[1] * k);
              rdata[idst] = varptr[isrc];
            }
          }
        }
        data[d] = rdata;
      }
      break;
    }

    case CCTK_VARIABLE_REAL: {
      for (int d = 0; d < num_comps; ++d) {
        CCTK_REAL const *const varptr =
            (CCTK_REAL const *)CCTK_VarDataPtrI(cctkGH, timelevel, var + d);
        assert(varptr);
        CCTK_REAL *const rdata = new CCTK_REAL[prod(ci.ilen)];
        for (int k = 0; k < ci.ilen[2]; ++k) {
          for (int j = 0; j < ci.ilen[1]; ++j) {
            for (int i = 0; i < ci.ilen[0]; ++i) {
              int const isrc =
                  (ci.imin[0] + i +
                   ci.ash[0] * (ci.imin[1] + j + ci.ash[1] * (ci.imin[2] + k)));
              int const idst = i + ci.ilen[0] * (j + ci.ilen[1] * k);
              rdata[idst] = varptr[isrc];
            }
          }
        }
        data[d] = rdata;
      }
      break;
    }

    default:
      assert(0);
    }

    // Dataset properties
    hid_t const prop = H5Pcreate(H5P_DATASET_CREATE);
    assert(prop >= 0);
    // Use chunked I/O if requested or needed
    bool const set_chunks =
        use_chunks or compression_level >= 0 or use_checksums;
    if (set_chunks) {
      assert(ci.dim > 0);
      int const chunksize = max_chunksize / H5Tget_size(type);
      assert(chunksize > 0);
      int cs = 1, csd = 1;
      while (csd < chunksize) {
        cs <<= 1;
        csd <<= ci.dim;
        assert(ci.dim > 0);
      }
      ivect chunkshape(cs);
      chunkshape[ci.dim - 1] = cs / (csd / chunksize);
      chunkshape = max(chunkshape, 1);
      chunkshape = min(chunkshape, ci.ilen);
      FAILWARN(
          H5Pset_chunk(prop, ci.dim, &v2h(chunkshape).reverse()[dim - ci.dim]));
    }
    // Enable compression if requested
    if (compression_level >= 0) {
      assert(set_chunks);
      indent_t indent2;
      cout << indent2 << "using chunk size " << ci.ilen << " "
           << "(" << prod(ci.ilen) << ")\n";
      FAILWARN(H5Pset_shuffle(prop));
      FAILWARN(H5Pset_deflate(prop, compression_level));
    }
    // Enable checksums if requested
    if (use_checksums) {
      assert(set_chunks);
      FAILWARN(H5Pset_fletcher32(prop));
    }

    if (num_comps == 1 and not separate_single_component_tensors) {
      // Write single-component tensors into non-separated fractions
      // for convenience (could also use a separated compound
      // instead)
      // TODO: Extend F5 API to allow writing non-fragmented datasets
      // TODO: Extend F5 API to accept hyperslab descriptors
      FAILWARN(F5Fwrite_fraction(
          path, name.c_str(), ci.dim, is_multipatch ? NULL : &v2h(ci.gsh)[0],
          &v2h(ci.ilen)[0], type, type, data[0], &v2h(ci.ioff)[0],
          &v2h(ci.lghosts)[0], &v2h(ci.ughosts)[0], fragmentname.c_str(),
          prop));
    } else {
      int const full_coverage = will_cover_complete_domain and
                                not fragment_contiguous_components and
                                all(slice_ipos < 0);
      // F5ls does not seem to support what F5 generates in this
      // case. It seems that F5 does not generate a separated object
      // (it does not generate a group for all datasets), but
      // declares that it does create one (it calls it
      // FragmentedSeparatedCompound). F5ls then gets confused
      // because it cannot open this group.
      if (not full_coverage and num_comps == 1) {
        CCTK_WARN(CCTK_WARN_ALERT,
                  "Outputting scalars in a fragmented, separated way. "
                  "This does not seem to be supported by F5ls "
                  "(or is implemented wrong in the F5 library).");
      }
      // TODO: Extend F5 API to accept hyperslab descriptors
      FAILWARN(F5FSwrite_fraction(
          path, name.c_str(), ci.dim, is_multipatch ? NULL : &v2h(ci.gsh)[0],
          &v2h(ci.ilen)[0], type, type, data, &v2h(ci.ioff)[0],
          &v2h(ci.lghosts)[0], &v2h(ci.ughosts)[0], fragmentname.c_str(), prop,
          full_coverage));
    }

    for (int d = 0; d < num_comps; ++d) {
      if (data[d]) {
        switch (vartype) {
        case CCTK_VARIABLE_INT:
          delete[](CCTK_INT const *) data[d];
          break;
        case CCTK_VARIABLE_REAL:
          delete[](CCTK_REAL const *) data[d];
          break;
        default:
          assert(0);
        }
      }
      data[d] = NULL;
    }

    FAILWARN(H5Pclose(prop));
  }

}; // class output_iterator_t

void write_metadata(cGH const *const cctkGH, hid_t const file) {
  DECLARE_CCTK_PARAMETERS;

  assert(cctkGH);
  assert(file >= 0);

  herr_t herr;

  // Create metadata only once
  // TODO: instead, overwrite the metadata
  htri_t const exists = H5Lexists(file, metadata_group, H5P_DEFAULT);
  assert(exists >= 0);
  if (exists)
    return;

  CCTK_INFO("Writing simulation metadata...");

  // Create a group to hold all metadata
  hid_t const group =
      H5Gcreate(file, metadata_group, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  assert(group >= 0);

  // General information
  WriteAttribute(group, "Cactus version", CCTK_FullVersion());

  // Unique identifiers
  if (CCTK_IsFunctionAliased("UniqueConfigID")) {
    WriteAttribute(group, "config id", (char const *)UniqueConfigID(cctkGH));
  }
  if (CCTK_IsFunctionAliased("UniqueBuildID")) {
    WriteAttribute(group, "build id", (char const *)UniqueBuildID(cctkGH));
  }
  if (CCTK_IsFunctionAliased("UniqueSimulationID")) {
    WriteAttribute(group, "simulation id",
                   (char const *)UniqueSimulationID(cctkGH));
  }
  if (CCTK_IsFunctionAliased("UniqueRunID")) {
    WriteAttribute(group, "run id", (char const *)UniqueRunID(cctkGH));
  }

// Don't write this attribute; the number of files may change
// after recovering, or after recombining files
#if 0
    // Number of I/O processes (i.e. the number of output files)
    int const nprocs = CCTK_nProcs(cctkGH);
    int const nioprocs = nprocs;
    WriteAttribute(group, "nioprocs", nioprocs);
#endif

  // Write parameters into a separate dataset (they may be too large
  // for an attribute)
  int const get_all = 1;
  char *const parameters = IOUtil_GetAllParameters(cctkGH, get_all);
  assert(parameters);
  WriteLargeAttribute(group, all_parameters, parameters);
  free(parameters);

  // Grid structure
  string const gs = serialise_grid_structure(cctkGH);
  WriteLargeAttribute(group, grid_structure, gs.c_str());

  herr = H5Gclose(group);
  assert(not herr);
}

void output(cGH const *const cctkGH, hid_t const file,
            vector<bool> const &output_var, bool const output_past_timelevels,
            bool const output_metadata) {
  output_iterator_t iterator(const_cast<cGH *>(cctkGH), output_var,
                             output_past_timelevels, output_metadata);
  iterator.iterate(file);
}

} // end namespace CarpetIOF5
