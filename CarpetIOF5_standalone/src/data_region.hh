#ifndef DATA_REGION_HH
#define DATA_REGION_HH

#include <string>

#include <hdf5.h>

#include "bbox.hh"
#include "defs.hh"

#include "physical_quantity.hh"

namespace CarpetIOF5 {

namespace F5 {

using std::string;

class data_region_t {

  physical_quantity_t &m_physical_quantity;

  bbox<int, dim> const m_region;
  string const m_name;

  hid_t m_hdf5_data_region;

  data_region_t();
  data_region_t(data_region_t const &);
  data_region_t operator=(data_region_t const &);

public:
  data_region_t(physical_quantity_t &physical_quantity,
                bbox<int, dim> const &region);

  virtual ~data_region_t();

  physical_quantity_t &get_physical_quantity() const;

  bbox<int, dim> const &get_region() const;

  hid_t get_hdf5_data_region() const;

  virtual bool invariant() const;
};

} // namespace F5

} // namespace CarpetIOF5

#endif // #ifndef DATA_REGION_HH
