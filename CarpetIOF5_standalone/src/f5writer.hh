#ifndef F5WRITER_HH
#define F5WRITER_HH

#include "cctk.h"

#include "bbox.hh"
#include "defs.hh"
#include "dh.hh"

#include "file.hh"
#include "physical_quantity.hh"
#include "simulation.hh"
#include "timestep.hh"
#include "topology.hh"

namespace CarpetIOF5 {

class f5writer_t {

  cGH const *const m_cctkGH;
  int const m_variable;

public:
  f5writer_t(cGH const *cctkGH, int variable);

  void write(F5::file_t &file) const;

private:
  void write_meta(F5::file_t &file, bool have_metafile) const;

  void write_one_mglevel(F5::timestep_t &timestep, bool have_metafile) const;

  void write_global(F5::simulation_t &simulation, bool have_metafile) const;

  void write_one_reflevel(F5::simulation_t &simulation,
                          bool have_metafile) const;

  void write_one_map(F5::simulation_t &simulation, bool have_metafile) const;

  void write_one_component(F5::physical_quantity_t &physical_quantity,
                           bool have_metafile) const;

  bbox<int, dim> const &determine_region(dh::dboxes const &boxes) const;
};

} // namespace CarpetIOF5

#endif // #ifndef F5WRITER_HH
