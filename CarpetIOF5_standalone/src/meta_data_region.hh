#ifndef META_DATA_REGION_HH
#define META_DATA_REGION_HH

#include <string>

#include <hdf5.h>

#include "bbox.hh"
#include "defs.hh"

#include "physical_quantity.hh"

namespace CarpetIOF5 {

namespace F5 {

using std::string;

class meta_data_region_t {

  physical_quantity_t &m_physical_quantity;

  bbox<int, dim> const m_region;
  string const m_name;

  hid_t m_properties;
  hid_t m_dataset;
  hid_t m_dataspace;

  meta_data_region_t();
  meta_data_region_t(meta_data_region_t const &);
  meta_data_region_t operator=(meta_data_region_t const &);

public:
  meta_data_region_t(physical_quantity_t &physical_quantity,
                     bbox<int, dim> const &region);

  virtual ~meta_data_region_t();

  physical_quantity_t &get_physical_quantity() const;

  void write(int proc) const;

  virtual bool invariant() const;
};

} // namespace F5

} // namespace CarpetIOF5

#endif // #ifndef META_DATA_REGION_HH
