#ifndef WRITER_HH
#define WRITER_HH

#include "cctk.h"

#include "bbox.hh"
#include "defs.hh"
#include "dh.hh"
#include "vect.hh"

#include "file.hh"
#include "physical_quantity.hh"
#include "simulation.hh"
#include "timestep.hh"
#include "topology.hh"

namespace CarpetIOF5 {

class writer_t {

  cGH const *const m_cctkGH;
  int const m_variable;

public:
  writer_t(cGH const *cctkGH, int variable);

  void write(F5::file_t &file) const;

private:
  void write_meta(F5::file_t &file) const;

  void write_one_mglevel(F5::timestep_t &timestep) const;

  void write_global(F5::timestep_t &timestep) const;

  void write_one_reflevel(F5::timestep_t &timestep) const;

  void write_one_map(F5::timestep_t &timestep) const;

  void write_one_component(F5::physical_quantity_t &physical_quantity) const;

  static bbox<int, dim> const &
  determine_region(dh::light_dboxes const &light_boxes);
};

} // namespace CarpetIOF5

#endif // #ifndef WRITER_HH
