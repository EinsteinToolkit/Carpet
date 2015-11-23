#include <cassert>
#include <cstdlib>
#include <sstream>
#include <string>

#include <hdf5.h>

#include "cctk.h"

#include "carpet.hh"

#include "data_region.hh"
#include "utils.hh"

namespace CarpetIOF5 {

namespace F5 {

using std::ostringstream;

data_region_t::data_region_t(physical_quantity_t &physical_quantity,
                             bbox<int, dim> const &region)
    : m_physical_quantity(physical_quantity), m_region(region),
      m_name(string("region-") + F5::name_from_ibbox(region)) {
  assert(not region.empty());

  m_hdf5_data_region = open_or_create_group(
      m_physical_quantity.get_hdf5_physical_quantity(), m_name.c_str());
  assert(m_hdf5_data_region >= 0);

  assert(invariant());
}

data_region_t::~data_region_t() {
  herr_t const herr = H5Gclose(m_hdf5_data_region);
  assert(not herr);
}

physical_quantity_t &data_region_t::get_physical_quantity() const {
  return m_physical_quantity;
}

bbox<int, dim> const &data_region_t::get_region() const { return m_region; }

hid_t data_region_t::get_hdf5_data_region() const { return m_hdf5_data_region; }

bool data_region_t::invariant() const {
  return (not m_region.empty() and m_hdf5_data_region >= 0);
}

} // namespace F5

} // namespace CarpetIOF5
