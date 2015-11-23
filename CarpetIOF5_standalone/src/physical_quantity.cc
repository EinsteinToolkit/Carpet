#include <cassert>
#include <cstdlib>

#include <hdf5.h>

#include "cctk.h"

#include "physical_quantity.hh"

namespace CarpetIOF5 {

namespace F5 {

physical_quantity_t::physical_quantity_t(coordinate_system_t &coordinate_system,
                                         int const group)
    : m_coordinate_system(coordinate_system), m_group(group) {
  assert(group >= 0 and group < CCTK_NumGroups());

  char *const name = CCTK_GroupName(group);
  assert(name != 0);
  for (char *p = name; *p; ++p) {
    *p = static_cast<char>(tolower(*p));
  }
  m_name = string(name);

  m_hdf5_physical_quantity = open_or_create_group(
      m_coordinate_system.get_hdf5_coordinate_system(), name);
  assert(m_hdf5_physical_quantity >= 0);

  free(name);

  assert(invariant());
}

physical_quantity_t::~physical_quantity_t() {
  herr_t const herr = H5Gclose(m_hdf5_physical_quantity);
  assert(not herr);
}

coordinate_system_t &physical_quantity_t::get_coordinate_system() const {
  return m_coordinate_system;
}

int physical_quantity_t::get_group() const { return m_group; }

string physical_quantity_t::get_name() const { return m_name; }

hid_t physical_quantity_t::get_hdf5_physical_quantity() const {
  return m_hdf5_physical_quantity;
}

void physical_quantity_t::get_link_destination(int const proc, string &filename,
                                               string &objectname) const {
  get_coordinate_system().get_link_destination(proc, filename, objectname);
  if (objectname.empty()) {
    objectname = m_name;
  } else {
    objectname += string("/") + m_name;
  }
}

bool physical_quantity_t::invariant() const {
  return (m_group >= 0 and m_group < CCTK_NumGroups() and
          m_hdf5_physical_quantity >= 0);
}

} // namespace F5

} // namespace CarpetIOF5
