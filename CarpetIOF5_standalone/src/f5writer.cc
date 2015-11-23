#include <sstream>
#include <string>

#include <hdf5.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "carpet.hh"

#include "data_region.hh"
#include "f5writer.hh"
#include "file.hh"
#include "meta_data_region.hh"
#include "physical_quantity.hh"
#include "simulation.hh"
#include "tensor_component.hh"
#include "timestep.hh"
#include "topology.hh"

namespace CarpetIOF5 {

f5writer_t::f5writer_t(cGH const *const cctkGH, int const variable)
    : m_cctkGH(cctkGH), m_variable(variable) {}

void f5writer_t::write(F5::file_t &file) const {
  write_meta(file, file.get_have_metafile());
}

void f5writer_t::write_meta(F5::file_t &file, bool const have_metafile) const {
  DECLARE_CCTK_PARAMETERS;

  if (verbose or veryverbose) {
    CCTK_VInfo(CCTK_THORNSTRING, "OutputVarAs/write_meta");
  }

  if (Carpet::is_meta_mode()) {
    for (Carpet::mglevel_iterator mglevel_iter(m_cctkGH);
         not mglevel_iter.done(); mglevel_iter.step()) {
      write_one_mglevel(file, have_metafile);
    }
  } else {
    write_one_mglevel(file, have_metafile);
  }
}

void f5writer_t::write_one_mglevel(F5::file_t &file,
                                   bool const have_metafile) const {
  DECLARE_CCTK_PARAMETERS;

  if (verbose or veryverbose) {
    CCTK_VInfo(CCTK_THORNSTRING, "OutputVarAs/write_one_mglevel mglevel=%d",
               Carpet::mglevel);
  }

  // ostringstream namebuf;
  // namebuf << "convlevel=" << m_cctkGH->cctk_convlevel;
  // string const namestr = namebuf.str();
  // char const * const name = namestr.c_str();

  int const grouptype = CCTK_GroupTypeFromVarI(m_variable);
  assert(grouptype >= 0);
  switch (grouptype) {
  case CCTK_ARRAY:
  case CCTK_SCALAR: {
    if (Carpet::do_global_mode) {
      write_global(file, have_metafile);
    }
  } break;
  case CCTK_GF: {
    if (Carpet::is_global_mode()) {
      for (Carpet::reflevel_iterator reflevel_iter(m_cctkGH);
           not reflevel_iter.done(); reflevel_iter.step()) {
        write_one_reflevel(file, have_metafile);
      }
    } else {
      write_one_reflevel(simulation, have_metafile);
    }
  } break;
  default:
    assert(0);
  }
}

void f5writer_t::write_global(F5::simulation_t &simulation,
                              bool const have_metafile) const {
  DECLARE_CCTK_PARAMETERS;

  if (verbose or veryverbose) {
    CCTK_INFO("OutputVarAs/write_global");
  }

  F5::unigrid_topology_t topology(simulation);

  int const grouptype = CCTK_GroupTypeFromVarI(m_variable);
  assert(grouptype >= 0);
  assert(grouptype == CCTK_SCALAR or grouptype == CCTK_ARRAY);

  vect<CCTK_REAL, dim> level_origin, level_delta;
  for (int d = 0; d < dim; ++d) {
    level_origin[d] = 0.0;
    level_delta[d] = 1.0;
  }
  F5::Cartesian_coordinate_system_t coordinate_system(topology, level_origin,
                                                      level_delta);

  int const group = CCTK_GroupIndexFromVarI(m_variable);
  assert(group >= 0 and group < CCTK_NumGroups());
  F5::physical_quantity_t physical_quantity(coordinate_system, group);

  int const map = 0;
  int const reflevel = 0;
  int const myproc = CCTK_MyProc(m_cctkGH);
  dh *const dd = Carpet::arrdata.at(group).at(map).dd;
  dh::dboxes const &boxes =
      dd->boxes.at(Carpet::mglevel).at(reflevel).at(myproc);
  bbox<int, dim> const &region = determine_region(boxes);

  if (have_metafile) {
    F5::meta_data_region_t meta_data_region(physical_quantity, region);
    gh *const hh = Carpet::vhh.at(Carpet::map);
    int const proc = hh->processor(Carpet::reflevel, Carpet::component);
    meta_data_region.write(proc);
  }

  F5::data_region_t data_region(physical_quantity, region);

  F5::tensor_component_t tensor_component(data_region, m_variable);
  int const timelevel = 0;
  void const *const varptr = CCTK_VarDataPtrI(m_cctkGH, timelevel, m_variable);
  assert(varptr != 0);
  int const vartype = CCTK_VarTypeI(m_variable);
  assert(vartype >= 0);
  tensor_component.write(varptr, vartype);
}

void f5writer_t::write_one_reflevel(F5::file_t &file,
                                    bool const have_metafile) const {
  DECLARE_CCTK_PARAMETERS;

  if (verbose or veryverbose) {
    CCTK_VInfo(CCTK_THORNSTRING, "OutputVarAs/write_one_reflevel reflevel=%d",
               Carpet::reflevel);
  }

  // int const grouptype = CCTK_GroupTypeFromVarI (m_variable);
  // assert (grouptype >= 0);
  // assert (grouptype == CCTK_GF);

  // A name for the simulation
  string sim_id;
  if (CCTK_IsFunctionAliased("UniqueSimulationID")) {
    sim_id = static_cast<char const *>(UniqueSimulationID(m_cctkGH));
  } else {
    char parfilename[10000];
    int const len = CCTK_ParameterFilename(sizeof parfilename, parfilename);
    assert(len >= 0 and len < sizeof parfilename);
    sim_id = parfilename;
    if (sim_id.length() >= 4 and sim_id.substr(sim_id.length() - 4) == ".par") {
      sim_id = sim_id.substr(0, sim_id.length() - 4);
    }
  }

  // Refinement factor of current refinement level
  vect<hsize_t, dim> const current_refinement =
      vect<int, dim>::ref(m_cctkGH->cctk_levfac);

  // Refinement factor of root refinement level
  vect<hsize_t, dim> const rootlevel_refinement(1);

  // Remember time when the root level was output the last time
  static CCTK_REAL rootlevel_time = 0.0;
  if (all(current_refinement == 1)) {
    rootlevel_time = m_cctkGH->cctk_time;
  }

  F5Path *refinement_path = NULL;
  if (any(current_refinement != 1)) {
    refinement_path = F5Rcreate_relative_vertex_Qrefinement3D(
        file, m_cctkGH->cctk_time, sim_id.c_str(), current_refinement,
        rootlevel_time, rootlevel_refinement);
  }

  // TODO: what is this doing?
  F5T_REFINEMENT3D_POINT = refinement_path->myChart->Point_hid_t;

  if (Carpet::is_level_mode()) {
    for (Carpet::map_iterator map_iter(m_cctkGH, grouptype);
         not map_iter.done(); map_iter.step()) {
      write_one_map(simulation, have_metafile);
    }
  } else {
    write_one_map(simulation, have_metafile);
  }
}

void f5writer_t::write_one_map(F5::simulation_t &simulation,
                               bool const have_metafile) const {
  DECLARE_CCTK_PARAMETERS;

  if (verbose or veryverbose) {
    CCTK_VInfo(CCTK_THORNSTRING, "OutputVarAs/write_one_map map=%d",
               Carpet::map);
  }

  // F5::mesh_refinement_topology_t topology
  //   (simulation, Carpet::map, Carpet::reflevel, Carpet::maxreflevels,
  //    Carpet::spacereflevelfact, Carpet::maxspacereflevelfact);

  // vect<CCTK_REAL, dim> level_origin, level_delta;
  // for (int d=0; d<dim; ++d)
  // {
  //   cGH const * const cctkGH = m_cctkGH;
  //   DECLARE_CCTK_ARGUMENTS;
  //   level_origin[d] = CCTK_ORIGIN_SPACE(d);
  //   level_delta[d]  = CCTK_DELTA_SPACE(d);
  // }
  // F5::Cartesian_coordinate_system_t coordinate_system
  //   (topology, level_origin, level_delta);
  //
  // int const group = CCTK_GroupIndexFromVarI (m_variable);
  // assert (group >= 0 and group < CCTK_NumGroups());
  // F5::physical_quantity_t physical_quantity (coordinate_system, group);

  // TODO: this should depend on the patch number
  char const *const coordinate_system = NULL;

  // Depends on the refinement level
  string const gridname = sim_id + "-" + current_refinement;

  F5Path *const myPath = F5Rcreate_vertex_refinement3D(
      file, m_cctkGH->cctk_time, sim_id.c_str(), refinement, coordinate_system);

  if (not refinement_path) {
    // Make root level appear as unigrid topology as well, just to
    // be nice to tools that do not understand AMR
    F5Rlink_default_vertex_topology(myPath, current_refinement);
  }

  // Global information for this refinement level
  vect<hsize_t, dim> const level_dims = vect<int, dim>::ref(m_cctkGH->cctk_gsh);
  // TODO: switch to double precision
  vect<float, dim> const level_min(CCTK_ORIGIN_SPACE(0), CCTK_ORIGIN_SPACE(1),
                                   CCTK_ORIGIN_SPACE(2));
  vect<float, dim> const level_spacing(CCTK_DELTA_SPACE(0), CCTK_DELTA_SPACE(1),
                                       CCTK_DELTA_SPACE(2));
  vect<float, dim> const level_max =
      level_min + level_spacing * (level_dims - 1);

  // Output geometric information for an entire level; it is
  // uniformely covered in coordinate space
  assert(dim == 3);
  F5Fwrite_linear(myPath, FIBER_HDF5_POSITIONS_STRING, dim, &level_dims[0],
                  F5T_COORD3_FLOAT, &level_min[0], &level_spacing[0]);

  F5Fset_range(myPath, &level_min[0], &level_max[0]);

#error "CONTINUE HERE"

  if (Carpet::is_singlemap_mode()) {
    int const grouptype = CCTK_GroupTypeI(group);
    assert(grouptype >= 0);

    for (Carpet::component_iterator component_iter(m_cctkGH, grouptype);
         not component_iter.done(); component_iter.step()) {
      write_one_component(physical_quantity, have_metafile);
    }
  } else {
    write_one_component(physical_quantity, have_metafile);
  }
}

void f5writer_t::write_one_component(F5::physical_quantity_t &physical_quantity,
                                     bool const have_metafile) const {
  DECLARE_CCTK_PARAMETERS;

  if (verbose or veryverbose) {
    CCTK_VInfo(CCTK_THORNSTRING, "OutputVarAs/write_one_component component=%d",
               Carpet::component);
  }

  gh *const hh = Carpet::vhh.at(Carpet::map);
  bool const is_local = hh->is_local(Carpet::reflevel, Carpet::component);
  if (have_metafile or is_local) {
    dh *const dd = Carpet::vdd.at(Carpet::map);
    bbox<int, dim> const &region = (dd->boxes.at(Carpet::mglevel)
                                        .at(Carpet::reflevel)
                                        .at(Carpet::component)
                                        .exterior);

    if (have_metafile) {
      F5::meta_data_region_t meta_data_region(physical_quantity, region);
      int const proc = hh->processor(Carpet::reflevel, Carpet::component);
      meta_data_region.write(proc);
    }

    if (is_local) {
      F5::data_region_t data_region(physical_quantity, region);

      F5::tensor_component_t tensor_component(data_region, m_variable);
      int const timelevel = 0;
      void const *const varptr =
          CCTK_VarDataPtrI(m_cctkGH, timelevel, m_variable);
      assert(varptr != 0);
      int const vartype = CCTK_VarTypeI(m_variable);
      assert(vartype >= 0);
      tensor_component.write(varptr, vartype);
    }
  }
}

bbox<int, dim> const &
f5writer_t::determine_region(dh::dboxes const &boxes) const {
  DECLARE_CCTK_PARAMETERS;

  // TODO: use superregions instead of regions (? only if the
  // regions are on the same processor?); use HDF5 chunks as well

  bbox<int, dim> dh::dboxes::*boxptr;
  if (CCTK_EQUALS(output_regions, "exterior")) {
    boxptr = &dh::dboxes::exterior;
  } else if (CCTK_EQUALS(output_regions, "owned")) {
    boxptr = &dh::dboxes::owned;
  } else if (CCTK_EQUALS(output_regions, "interior")) {
    boxptr = &dh::dboxes::interior;
  } else {
    assert(0);
  }

  return boxes.*boxptr;
}

} // namespace CarpetIOF5
