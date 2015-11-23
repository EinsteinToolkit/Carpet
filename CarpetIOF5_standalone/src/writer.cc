#include <cstdlib>
#include <sstream>
#include <string>

#include <hdf5.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "carpet.hh"
#include "modes.hh"

#include "data_region.hh"
#include "file.hh"
#include "meta_data_region.hh"
#include "physical_quantity.hh"
#include "simulation.hh"
#include "tensor_component.hh"
#include "timestep.hh"
#include "topology.hh"
#include "utils.hh"
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
    BEGIN_MGLEVEL_LOOP(m_cctkGH) { write_one_mglevel(timestep); }
    END_MGLEVEL_LOOP;
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

  int const grouptype = CCTK_GroupTypeFromVarI(m_variable);
  assert(grouptype >= 0);
  switch (grouptype) {
  case CCTK_ARRAY:
  case CCTK_SCALAR: {
    if (Carpet::is_global_mode()) {
      ENTER_LEVEL_MODE(m_cctkGH, 0) { write_global(timestep); }
      LEAVE_LEVEL_MODE;
    } else {
      if (Carpet::do_global_mode) {
        write_global(timestep);
      }
    }
  } break;
  case CCTK_GF: {
    if (Carpet::is_global_mode()) {
      BEGIN_REFLEVEL_LOOP(m_cctkGH) { write_one_reflevel(timestep); }
      END_REFLEVEL_LOOP;
    } else {
      write_one_reflevel(timestep);
    }
  } break;
  default:
    assert(0);
  }
}

void writer_t::write_global(F5::timestep_t &timestep) const {
  DECLARE_CCTK_PARAMETERS;

  if (verbose or veryverbose) {
    CCTK_VInfo(CCTK_THORNSTRING, "OutputVarAs/write_global");
  }

  int const grouptype = CCTK_GroupTypeFromVarI(m_variable);
  assert(grouptype == CCTK_SCALAR or grouptype == CCTK_ARRAY);

  BEGIN_MAP_LOOP(m_cctkGH, grouptype) {

    // Name the grid after the variable group
    ostringstream namebuf;
    char *const c_name = CCTK_GroupNameFromVarI(m_variable);
    namebuf << "Cactus-" << c_name;
    free(c_name);
    string const name = namebuf.str();
    F5::simulation_t simulation(timestep, name.c_str());

    F5::unigrid_topology_t topology(simulation);

    vect<CCTK_REAL, dim> const level_origin(0.0), level_delta(1.0);
    F5::Cartesian_coordinate_system_t coordinate_system(topology, level_origin,
                                                        level_delta);

    int const group = CCTK_GroupIndexFromVarI(m_variable);
    assert(group >= 0 and group < CCTK_NumGroups());
    F5::physical_quantity_t physical_quantity(coordinate_system, group);

    int const myproc = CCTK_MyProc(m_cctkGH);

    F5::file_t &file = timestep.get_file();
    bool const write_metafile =
        file.get_is_metafile() and not file.get_is_datafile();
    if (write_metafile) {

      BEGIN_COMPONENT_LOOP(m_cctkGH, grouptype) {
        dh *const dd = Carpet::arrdata.at(group).at(Carpet::map).dd;
        dh::light_dboxes const &light_boxes =
            dd->light_boxes.at(Carpet::mglevel).at(Carpet::reflevel).at(myproc);
        bbox<int, dim> const &region = determine_region(light_boxes);
        F5::meta_data_region_t meta_data_region(physical_quantity, region);

        gh *const hh = Carpet::vhh.at(Carpet::map);
        int const proc = hh->processor(Carpet::reflevel, Carpet::component);
        meta_data_region.write(proc);
      }
      END_COMPONENT_LOOP;

    } else // if not write_metafile
    {

      BEGIN_LOCAL_COMPONENT_LOOP(m_cctkGH, grouptype) {
        dh *const dd = Carpet::arrdata.at(group).at(Carpet::map).dd;
        dh::light_dboxes const &light_boxes =
            dd->light_boxes.at(Carpet::mglevel).at(Carpet::reflevel).at(myproc);
        bbox<int, dim> const &region = determine_region(light_boxes);
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
      END_LOCAL_COMPONENT_LOOP;

    } // if not write_metafile
  }
  END_MAP_LOOP;
}

void writer_t::write_one_reflevel(F5::timestep_t &timestep) const {
  DECLARE_CCTK_PARAMETERS;

  if (verbose or veryverbose) {
    CCTK_VInfo(CCTK_THORNSTRING, "OutputVarAs/write_one_reflevel reflevel=%d",
               Carpet::reflevel);
  }

  int const grouptype = CCTK_GroupTypeFromVarI(m_variable);
  assert(grouptype == CCTK_GF);

  if (Carpet::is_level_mode()) {
    BEGIN_MAP_LOOP(m_cctkGH, grouptype) { write_one_map(timestep); }
    END_MAP_LOOP;
  } else {
    write_one_map(timestep);
  }
}

void writer_t::write_one_map(F5::timestep_t &timestep) const {
  DECLARE_CCTK_PARAMETERS;

  if (verbose or veryverbose) {
    CCTK_VInfo(CCTK_THORNSTRING, "OutputVarAs/write_one_map map=%d",
               Carpet::map);
  }

  // Name the grid after the map number
  ostringstream namebuf;
  namebuf << "Carpet";
  if (Carpet::maps > 1) {
    namebuf << "-map" << Carpet::map;
  }
  string const name = namebuf.str();
  F5::simulation_t simulation(timestep, name.c_str());

  F5::mesh_refinement_topology_t topology(
      simulation, Carpet::map, Carpet::maps, Carpet::reflevel,
      Carpet::maxreflevels, Carpet::spacereflevelfact,
      Carpet::maxspacereflevelfact);

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

    BEGIN_COMPONENT_LOOP(m_cctkGH, grouptype) {
      write_one_component(physical_quantity);
    }
    END_COMPONENT_LOOP;
  } else {
    write_one_component(physical_quantity);
  }
}

void writer_t::write_one_component(
    F5::physical_quantity_t &physical_quantity) const {
  DECLARE_CCTK_PARAMETERS;

  if (verbose or veryverbose) {
    CCTK_VInfo(CCTK_THORNSTRING, "OutputVarAs/write_one_component component=%d",
               Carpet::component);
  }

  gh *const hh = Carpet::vhh.at(Carpet::map);
  bool const is_local = hh->is_local(Carpet::reflevel, Carpet::component);

  F5::file_t &file = (physical_quantity.get_coordinate_system()
                          .get_topology()
                          .get_simulation()
                          .get_timestep()
                          .get_file());
  bool const write_metafile =
      file.get_is_metafile() and not file.get_is_datafile();
  if (write_metafile or is_local) {
    dh *const dd = Carpet::vdd.at(Carpet::map);
    bbox<int, dim> const &region = (dd->light_boxes.at(Carpet::mglevel)
                                        .at(Carpet::reflevel)
                                        .at(Carpet::component)
                                        .exterior);

    if (write_metafile) {
      F5::meta_data_region_t meta_data_region(physical_quantity, region);
      int const proc = hh->processor(Carpet::reflevel, Carpet::component);
      meta_data_region.write(proc);
    } else if (is_local) {
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
writer_t::determine_region(dh::light_dboxes const &light_boxes) {
  DECLARE_CCTK_PARAMETERS;

  // TODO: use superregions instead of regions (?  only if the
  // regions are on the same processor?)

  bbox<int, dim> dh::light_dboxes::*boxptr;
  if (CCTK_EQUALS(output_regions, "exterior")) {
    boxptr = &dh::light_dboxes::exterior;
  } else if (CCTK_EQUALS(output_regions, "owned")) {
    boxptr = &dh::light_dboxes::owned;
  } else if (CCTK_EQUALS(output_regions, "interior")) {
    boxptr = &dh::light_dboxes::interior;
  } else {
    assert(0);
  }

  return light_boxes.*boxptr;
}

} // namespace CarpetIOF5
