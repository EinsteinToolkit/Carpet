#include <cassert>
#include <cstdlib>
#include <sstream>
#include <string>

#include <hdf5.h>

#include "cctk.h"

#include "carpet.hh"

#include "defs.hh"

#include "meta_data_region.hh"
#include "utils.hh"

namespace CarpetIOF5 {

namespace F5 {

using std::ostringstream;

meta_data_region_t::meta_data_region_t(physical_quantity_t &physical_quantity,
                                       bbox<int, dim> const &region)
    : m_physical_quantity(physical_quantity), m_region(region),
      m_name(string("region-") + F5::name_from_ibbox(region)) {
  assert(not region.empty());

  assert(invariant());
}

meta_data_region_t::~meta_data_region_t() {}

physical_quantity_t &meta_data_region_t::get_physical_quantity() const {
  return m_physical_quantity;
}

void meta_data_region_t::write(int const proc) const {
  DECLARE_CCTK_PARAMETERS;

  string filename;
  string objectname;
  get_physical_quantity().get_link_destination(proc, filename, objectname);

  hid_t const hdf5_physical_quantity =
      m_physical_quantity.get_hdf5_physical_quantity();

  bool link_exists;
  H5E_BEGIN_TRY {
    H5L_info_t linkbuf;
    herr_t const herr = H5Lget_info(hdf5_physical_quantity, m_name.c_str(),
                                    &linkbuf, H5P_DEFAULT);
    link_exists = herr == 0;
  }
  H5E_END_TRY;
  if (veryverbose) {
    CCTK_VInfo(CCTK_THORNSTRING, "Link \"%s\" %s", m_name.c_str(),
               link_exists ? "exists" : "does not exist");
  }

  if (not link_exists) {
    if (veryverbose) {
      CCTK_VInfo(CCTK_THORNSTRING,
                 "H5Lcreate_external "
                 "(filename=\"%s\", objectname=\"%s\", linkname=\"%s\")",
                 filename.c_str(), objectname.c_str(), m_name.c_str());
    }
    check(not H5Lcreate_external(filename.c_str(), objectname.c_str(),
                                 hdf5_physical_quantity, m_name.c_str(),
                                 H5P_DEFAULT, H5P_DEFAULT));
  }
}

bool meta_data_region_t::invariant() const { return not m_region.empty(); }

} // namespace F5

} // namespace CarpetIOF5
