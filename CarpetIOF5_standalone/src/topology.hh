#ifndef TOPOLOGY_HH
#define TOPOLOGY_HH

#include <hdf5.h>

#include "cctk.h"

#include "defs.hh"
#include "vect.hh"

#include "simulation.hh"
#include "utils.hh"

namespace CarpetIOF5 {

namespace F5 {

class topology_t {

  topology_t();
  topology_t(topology_t const &);
  topology_t operator=(topology_t const &);

protected:
  simulation_t &m_simulation;

  string m_name;

  hid_t m_hdf5_topology;

  topology_t(simulation_t &simulation);

public:
  virtual ~topology_t();

public:
  simulation_t &get_simulation() const;

  void get_link_destination(int proc, string &filename,
                            string &objectname) const;

  hid_t get_hdf5_topology() const;

  virtual bool invariant() const;
};

class unigrid_topology_t : public topology_t {

  unigrid_topology_t();
  unigrid_topology_t(unigrid_topology_t const &);
  unigrid_topology_t operator=(unigrid_topology_t const &);

public:
  // Create unigrid topology
  unigrid_topology_t(simulation_t &simulation);

  virtual ~unigrid_topology_t();

  virtual bool invariant() const;
};

class mesh_refinement_topology_t : public topology_t {

  int const m_refinement_level;
  int const m_max_refinement_levels;
  vect<int, dim> const m_level_refinement_factor;
  vect<int, dim> const m_max_refinement_factor;

  mesh_refinement_topology_t();
  mesh_refinement_topology_t(mesh_refinement_topology_t const &);
  mesh_refinement_topology_t operator=(mesh_refinement_topology_t const &);

public:
  mesh_refinement_topology_t(simulation_t &simulation, int map, int maps,
                             int refinement_level, int max_refinement_levels,
                             vect<int, dim> const &level_refinement_factors,
                             vect<int, dim> const &max_refinement_factors);

  void
  calculate_level_origin_delta(vect<CCTK_REAL, dim> const &coarse_origin,
                               vect<CCTK_REAL, dim> const &coarse_delta,
                               vect<int, dim> const &level_offset,
                               vect<int, dim> const &level_offset_denominator,
                               vect<CCTK_REAL, dim> &level_origin,
                               vect<CCTK_REAL, dim> &level_delta) const;

  virtual ~mesh_refinement_topology_t();

  void get_link_destination(string &filename, string &objectname) const;

  virtual bool invariant() const;
};

} // namespace F5

} // namespace CarpetIOF5

#endif // #ifndef TOPOLOGY_HH
