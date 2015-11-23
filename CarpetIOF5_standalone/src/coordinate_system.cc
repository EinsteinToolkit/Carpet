#include <sstream>
#include <string>

#include <hdf5.h>

#include "cctk.h"

#include "coordinate_system.hh"
#include "utils.hh"

namespace CarpetIOF5 {

namespace F5 {

coordinate_system_t::coordinate_system_t(topology_t &topology)
    : m_topology(topology) {}

coordinate_system_t::~coordinate_system_t() {}

topology_t &coordinate_system_t::get_topology() const { return m_topology; }

hid_t coordinate_system_t::get_hdf5_coordinate_system() const {
  return m_hdf5_coordinate_system;
}

bool coordinate_system_t::invariant() const { return true; }

Cartesian_coordinate_system_t::Cartesian_coordinate_system_t(
    topology_t &topology, vect<CCTK_REAL, dim> const &level_origin,
    vect<CCTK_REAL, dim> const &level_delta)
    : coordinate_system_t(topology), m_level_origin(level_origin),
      m_level_delta(level_delta) {
  assert(all(m_level_delta > (CCTK_REAL)0));

  init();

  assert(invariant());
}

Cartesian_coordinate_system_t::Cartesian_coordinate_system_t(
    topology_t &topology, vect<CCTK_REAL, dim> const &coarse_origin,
    vect<CCTK_REAL, dim> const &coarse_delta,
    vect<int, dim> const &level_offset,
    vect<int, dim> const &level_offset_denominator)
    : coordinate_system_t(topology) {
  assert(all(coarse_delta > (CCTK_REAL)0));
  assert(all(level_offset_denominator > 0));

  // mesh_refinement_topology_t * mesh_refinement_topology
  //   = dynamic_cast<mesh_refinement_topology_t *> (& topology);
  mesh_refinement_topology_t *mesh_refinement_topology =
      (mesh_refinement_topology_t *)(&topology);
  assert(mesh_refinement_topology != 0);
  mesh_refinement_topology->calculate_level_origin_delta(
      coarse_origin, coarse_delta, level_offset, level_offset_denominator,
      m_level_origin, m_level_delta);

  init();

  assert(invariant());
}

Cartesian_coordinate_system_t::~Cartesian_coordinate_system_t() {
  check(not H5Gclose(m_hdf5_coordinate_system));
}

void Cartesian_coordinate_system_t::init() {
  assert(all(m_level_delta > (CCTK_REAL)0));

  ostringstream namebuf;
  namebuf << "Cartesian3D-" << name_from_ivect(m_level_origin) << "-"
          << name_from_ivect(m_level_delta);
  string const namestr = namebuf.str();
  m_name = namestr;
  char const *const name = namestr.c_str();

  m_hdf5_coordinate_system =
      open_or_create_group(m_topology.get_hdf5_topology(), name);
  assert(m_hdf5_coordinate_system >= 0);

  // TODO: don't output coordinates as attributes
  write_or_check_attribute(m_hdf5_coordinate_system, "origin", m_level_origin);
  write_or_check_attribute(m_hdf5_coordinate_system, "delta", m_level_delta);
}

void coordinate_system_t::get_link_destination(int const proc, string &filename,
                                               string &objectname) const {
  get_topology().get_link_destination(proc, filename, objectname);
  if (objectname.empty()) {
    objectname = m_name;
  } else {
    objectname += string("/") + m_name;
  }
}

bool Cartesian_coordinate_system_t::invariant() const {
  return (coordinate_system_t::invariant() and
          all(m_level_delta > (CCTK_REAL)0));
}

} // namespace F5

} // namespace CarpetIOF5
