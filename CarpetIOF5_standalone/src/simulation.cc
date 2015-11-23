#include <cassert>

#include <hdf5.h>

#include "cctk.h"

#include "defs.hh"

#include "simulation.hh"
#include "utils.hh"

namespace CarpetIOF5 {

namespace F5 {

using std::ostringstream;

simulation_t::simulation_t(timestep_t &timestep, char const *const name)
    : m_timestep(timestep), m_name(name) {
  m_hdf5_simulation =
      open_or_create_group(m_timestep.get_hdf5_timestep(), m_name.c_str());
  assert(m_hdf5_simulation >= 0);

  cGH const *const cctkGH = m_timestep.get_file().get_cctkGH();
  write_or_check_attribute(m_hdf5_simulation, "TimeStep",
                           cctkGH->cctk_iteration);

  assert(invariant());
}

simulation_t::~simulation_t() { check(not H5Gclose(m_hdf5_simulation)); }

timestep_t &simulation_t::get_timestep() const { return m_timestep; }

hid_t simulation_t::get_hdf5_simulation() const { return m_hdf5_simulation; }

void simulation_t::get_link_destination(int const proc, string &filename,
                                        string &objectname) const {
  get_timestep().get_link_destination(proc, filename, objectname);
  if (objectname.empty()) {
    objectname = m_name;
  } else {
    objectname += string("/") + m_name;
  }
}

bool simulation_t::invariant() const { return m_hdf5_simulation >= 0; }

} // namespace F5

} // namespace CarpetIOF5
