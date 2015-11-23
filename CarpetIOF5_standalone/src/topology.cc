#include <cassert>
#include <sstream>

#include <hdf5.h>

#include "cctk.h"

#include "defs.hh"
#include "vect.hh"

#include "topology.hh"

namespace CarpetIOF5 {

namespace F5 {

topology_t::topology_t(simulation_t &simulation) : m_simulation(simulation) {}

topology_t::~topology_t() {}

simulation_t &topology_t::get_simulation() const { return m_simulation; }

hid_t topology_t::get_hdf5_topology() const { return m_hdf5_topology; }

bool topology_t::invariant() const { return m_hdf5_topology >= 0; }

unigrid_topology_t::unigrid_topology_t(simulation_t &simulation)
    : topology_t(simulation) {
  char const *const name = "uniform";
  m_name = string(name);
  m_hdf5_topology =
      open_or_create_group(m_simulation.get_hdf5_simulation(), name);
  assert(m_hdf5_topology >= 0);

  assert(invariant());
}

unigrid_topology_t::~unigrid_topology_t() {
  herr_t const herr = H5Gclose(m_hdf5_topology);
  assert(not herr);
}

bool unigrid_topology_t::invariant() const { return topology_t::invariant(); }

mesh_refinement_topology_t::mesh_refinement_topology_t(
    simulation_t &simulation, int const map, int const maps,
    int const refinement_level, int const max_refinement_levels,
    vect<int, dim> const &level_refinement_factor,
    vect<int, dim> const &max_refinement_factor)
    : topology_t(simulation), m_refinement_level(refinement_level),
      m_max_refinement_levels(max_refinement_levels),
      m_level_refinement_factor(level_refinement_factor),
      m_max_refinement_factor(max_refinement_factor) {
  assert(refinement_level >= 0);
  assert(refinement_level < max_refinement_levels);
  assert(all(level_refinement_factor > 0));
  assert(all(level_refinement_factor <= max_refinement_factor));

  ostringstream namebuf;
  namebuf << "Vertices";
  if (maps > 1) {
    namebuf << "-map" << map;
  }
  if (max_refinement_levels > 1) {
    namebuf << "-level" << refinement_level;
  }
  string const namestr = namebuf.str();
  m_name = namestr;
  char const *const name = namestr.c_str();

  m_hdf5_topology =
      open_or_create_group(m_simulation.get_hdf5_simulation(), name);
  assert(m_hdf5_topology >= 0);

  assert(invariant());
}

void mesh_refinement_topology_t::calculate_level_origin_delta(
    vect<CCTK_REAL, dim> const &coarse_origin,
    vect<CCTK_REAL, dim> const &coarse_delta,
    vect<int, dim> const &level_offset,
    vect<int, dim> const &level_offset_denominator,
    vect<CCTK_REAL, dim> &level_origin,
    vect<CCTK_REAL, dim> &level_delta) const {
  assert(all(coarse_delta > (CCTK_REAL)0));
  assert(all(level_offset_denominator > 0));

  level_delta = coarse_delta / vect<CCTK_REAL, dim>(m_level_refinement_factor);
  level_origin =
      (coarse_origin + (vect<CCTK_REAL, dim>(level_offset) /
                        vect<CCTK_REAL, dim>(level_offset_denominator) /
                        vect<CCTK_REAL, dim>(m_max_refinement_factor)));
}

mesh_refinement_topology_t::~mesh_refinement_topology_t() {
  herr_t const herr = H5Gclose(m_hdf5_topology);
  assert(not herr);
}

void topology_t::get_link_destination(int const proc, string &filename,
                                      string &objectname) const {
  get_simulation().get_link_destination(proc, filename, objectname);
  if (objectname.empty()) {
    objectname = m_name;
  } else {
    objectname += string("/") + m_name;
  }
}

bool mesh_refinement_topology_t::invariant() const {
  return (topology_t::invariant() and m_refinement_level >= 0 and
          m_refinement_level < m_max_refinement_levels and
          all(m_level_refinement_factor > 0) and
          all(m_level_refinement_factor <= m_max_refinement_factor));
}

} // namespace F5

} // namespace CarpetIOF5
