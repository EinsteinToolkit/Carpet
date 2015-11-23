#include <cstdlib>
#include <sstream>
#include <string>

#include <hdf5.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "carpet.hh"

#include "data_region.hh"
#include "file.hh"
#include "meta_data_region.hh"
#include "physical_quantity.hh"
#include "simulation.hh"
#include "tensor_component.hh"
#include "timestep.hh"
#include "topology.hh"
#include "writer.hh"

namespace CarpetIOF5 {

writer_t::writer_t(cGH const *const cctkGH, int const variable)
    : m_cctkGH(cctkGH), m_variable(variable) {}

void writer_t::write(F5::file_t &file) const { write_meta(file); }

void writer_t::write_meta(F5::file_t &file) const {
  DECLARE_CCTK_PARAMETERS;

  if (verbose or veryverbose) {
    CCTK_VInfo(CCTK_THORNSTRING, "OutputVarAs/write_meta");
  }

  F5::timestep_t timestep(file, m_cctkGH->cctk_time);

  if (Carpet::is_meta_mode()) {
    for (Carpet::mglevel_iterator mglevel_iter(m_cctkGH);
         not mglevel_iter.done(); mglevel_iter.step()) {
      write_one_mglevel(timestep);
    }
  } else {
    write_one_mglevel(timestep);
  }
}

void writer_t::write_one_mglevel(F5::timestep_t &timestep) const {
  DECLARE_CCTK_PARAMETERS;

  if (verbose or veryverbose) {
    CCTK_VInfo(CCTK_THORNSTRING, "OutputVarAs/write_one_mglevel mglevel=%d",
               Carpet::mglevel);
  }

  // ostringstream namebuf;
  // namebuf << sim_id;
  // if (m_cctkGH->cctk_convlevel != 0)
  // {
  //   namebuf << "-convlevel" << m_cctkGH->cctk_convlevel;
  // }
  // string const name = namebuf.str();
  // F5::simulation_t simulation (timestep, name.c_str());

  int const grouptype = CCTK_GroupTypeFromVarI(m_variable);
  assert(grouptype >= 0);
  switch (grouptype) {
  case CCTK_ARRAY:
  case CCTK_SCALAR: {
    if (Carpet::is_global_mode()) {
      int const rl = 0;
      enter_level_mode(m_cctkGH, rl);
      write_global_reflevel(timestep);
      leave_level_mode(m_cctkGH);
    } else {
      if (Carpet::do_global_mode) {
        write_global_reflevel(timestep);
      }
    }
  } break;
  case CCTK_GF: {
    if (Carpet::is_global_mode()) {
      for (Carpet::reflevel_iterator reflevel_iter(m_cctkGH);
           not reflevel_iter.done(); reflevel_iter.step()) {
        write_one_reflevel(timestep);
      }
    } else {
      write_one_reflevel(simulation);
    }
  } break;
  default:
    assert(0);
  }
}

void writer_t::write_global_reflevel(F5::timestep_t &timestep) const {
  DECLARE_CCTK_PARAMETERS;

  if (verbose or veryverbose) {
    CCTK_INFO("OutputVarAs/write_global_reflevel");
  }

  if (Carpet::is_level_mode()) {
    for (Carpet::map_iterator map_iter(m_cctkGH, grouptype);
         not map_iter.done(); map_iter.step()) {
      write_global_map(timestep);
    }
  } else {
    write_global_map(timestep);
  }
}

void writer_t::write_global_map(F5::timestep_t &timestep) const {
  DECLARE_CCTK_PARAMETERS;

  if (verbose or veryverbose) {
    CCTK_INFO("OutputVarAs/write_global_map");
  }

  write_global_component(timestep);
}

void writer_t::write_global_component(F5::timestep_t &timestep) const {
  DECLARE_CCTK_PARAMETERS;

  if (verbose or veryverbose) {
    CCTK_INFO("OutputVarAs/write_global_component");
  }

  int const group = CCTK_GroupIndexFromVarI(m_variable);
  assert(group >= 0 and group < CCTK_NumGroups());
  cGroup groupdata;
  check(not CCTK_GropuData(group, &groupdata));
  assert(groupdata.grouptype == CCTK_SCALAR or
         groupdata.grouptype == CCTK_ARRAY);

  for (Carpet::component_iterator component_iter(m_cctkGH, grouptype);
       not component_iter.done(); component_iter.step())

    int const myproc = CCTK_MyProc(m_cctkGH);
  if (groupdata.disttype == CCTK_DISTRIB_CONSTANT and myproc != 0) {
    // Output DISTRIB=constant groups only on the root processor
    return;
  }

  // Name the grid after the variable group
  char *const c_name = CCTK_GroupNameFromVarI(m_variable);
  string const name = c_name;
  free(c_name);
  F5::simulation_t simulation(timestep, name.c_str());

  F5::unigrid_topology_t topology(simulation);

  vect<CCTK_REAL, dim> const level_origin(0.0), level_delta(1.0);
  F5::Cartesian_coordinate_system_t coordinate_system(topology, level_origin,
                                                      level_delta);

  F5::physical_quantity_t physical_quantity(coordinate_system, group);

  int const map = 0;
  int const reflevel = 0;
  dh *const dd = Carpet::arrdata.at(group).at(map).dd;
  dh::dboxes const &boxes =
      dd->boxes.at(Carpet::mglevel).at(reflevel).at(myproc);
  bbox<int, dim> const &region = determine_region(boxes);

  F5::file_t &file = timestep.get_file();
  bool const have_metafile =
      file.get_is_metafile() and not file.get_is_datafile();
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

void writer_t::write_global(F5::timestep_t &timestep) const {
  DECLARE_CCTK_PARAMETERS;

  if (verbose or veryverbose) {
    CCTK_INFO("OutputVarAs/write_global");
  }

  int const group = CCTK_GroupIndexFromVarI(m_variable);
  assert(group >= 0 and group < CCTK_NumGroups());
  cGroup groupdata;
  check(not CCTK_GropuData(group, &groupdata));
  assert(groupdata.grouptype == CCTK_SCALAR or
         groupdata.grouptype == CCTK_ARRAY);

  for (Carpet::component_iterator component_iter(m_cctkGH, grouptype);
       not component_iter.done(); component_iter.step())

    int const myproc = CCTK_MyProc(m_cctkGH);
  if (groupdata.disttype == CCTK_DISTRIB_CONSTANT and myproc != 0) {
    // Output DISTRIB=constant groups only on the root processor
    return;
  }

  // Name the grid after the variable group
  char *const c_name = CCTK_GroupNameFromVarI(m_variable);
  string const name = c_name;
  free(c_name);
  F5::simulation_t simulation(timestep, name.c_str());

  F5::unigrid_topology_t topology(simulation);

  vect<CCTK_REAL, dim> const level_origin(0.0), level_delta(1.0);
  F5::Cartesian_coordinate_system_t coordinate_system(topology, level_origin,
                                                      level_delta);

  F5::physical_quantity_t physical_quantity(coordinate_system, group);

  int const map = 0;
  int const reflevel = 0;
  dh *const dd = Carpet::arrdata.at(group).at(map).dd;
  dh::dboxes const &boxes =
      dd->boxes.at(Carpet::mglevel).at(reflevel).at(myproc);
  bbox<int, dim> const &region = determine_region(boxes);

  F5::file_t &file = timestep.get_file();
  bool const have_metafile =
      file.get_is_metafile() and not file.get_is_datafile();
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

void writer_t::write_one_map(F5::simulation_t &simulation) const {
  DECLARE_CCTK_PARAMETERS;

  if (verbose or veryverbose) {
    CCTK_VInfo(CCTK_THORNSTRING, "OutputVarAs/write_one_map map=%d",
               Carpet::map);
  }

  F5::mesh_refinement_topology_t topology(
      simulation, Carpet::map, Carpet::reflevel, Carpet::maxreflevels,
      Carpet::spacereflevelfact, Carpet::maxspacereflevelfact);

  vect<CCTK_REAL, dim> level_origin, level_delta;
  for (int d = 0; d < dim; ++d) {
    cGH const *const cctkGH = m_cctkGH;
    DECLARE_CCTK_ARGUMENTS;
    level_origin[d] = CCTK_ORIGIN_SPACE(d);
    level_delta[d] = CCTK_DELTA_SPACE(d);
  }
  F5::Cartesian_coordinate_system_t coordinate_system(topology, level_origin,
                                                      level_delta);

  int const group = CCTK_GroupIndexFromVarI(m_variable);
  assert(group >= 0 and group < CCTK_NumGroups());
  F5::physical_quantity_t physical_quantity(coordinate_system, group);

  if (Carpet::is_singlemap_mode()) {
    int const grouptype = CCTK_GroupTypeI(group);
    assert(grouptype >= 0);

    for (Carpet::component_iterator component_iter(m_cctkGH, grouptype);
         not component_iter.done(); component_iter.step()) {
      write_one_component(physical_quantity);
    }
  } else {
    write_one_component(physical_quantity);
  }
}

void writer_t::write_one_reflevel(F5::simulation_t &simulation) const {
  DECLARE_CCTK_PARAMETERS;

  if (verbose or veryverbose) {
    CCTK_VInfo(CCTK_THORNSTRING, "OutputVarAs/write_one_reflevel reflevel=%d",
               Carpet::reflevel);
  }

  int const grouptype = CCTK_GroupTypeFromVarI(m_variable);
  assert(grouptype >= 0);
  assert(grouptype == CCTK_GF);

  if (Carpet::is_level_mode()) {
    for (Carpet::map_iterator map_iter(m_cctkGH, grouptype);
         not map_iter.done(); map_iter.step()) {
      write_one_map(simulation);
    }
  } else {
    write_one_map(simulation);
  }
}

void writer_t::write_one_component(
    F5::physical_quantity_t &physical_quantity) const {
  DECLARE_CCTK_PARAMETERS;

  if (verbose or veryverbose) {
    CCTK_VInfo(CCTK_THORNSTRING, "OutputVarAs/write_one_component component=%d",
               Carpet::component);
  }

  F5::file_t &file = (physical_quantity.get_coordinate_system()
                          .get_topology()
                          .get_simulation()
                          .get_timestep()
                          .get_file());
  bool const have_metafile =
      file.get_is_metafile() and not file.get_is_datafile();
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
writer_t::determine_region(dh::dboxes const &boxes) const {
  DECLARE_CCTK_PARAMETERS;

  // TODO: use superregions instead of regions (?  only if the
  // regions are on the same processor?)

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
