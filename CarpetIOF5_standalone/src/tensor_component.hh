#ifndef TENSOR_COMPONENT_HH
#define TENSOR_COMPONENT_HH

#include <string>

#include <hdf5.h>

#include "data_region.hh"

namespace CarpetIOF5 {

namespace F5 {

class tensor_component_t {

  data_region_t &m_data_region;

  int const m_variable;
  string m_name;

  hid_t m_dataspace;
  hid_t m_properties;
  hid_t m_dataset;

  tensor_component_t();
  tensor_component_t(tensor_component_t const &);
  tensor_component_t operator=(tensor_component_t const &);

public:
  tensor_component_t(data_region_t &data_region, int variable);

  virtual ~tensor_component_t();

  data_region_t &get_data_region() const;

  hid_t get_variable() const;

  void write(void const *data, int cactus_datatype) const;

  virtual bool invariant() const;
};

} // namespace F5

} // namespace CarpetIOF5

#endif // #ifndef TENSOR_COMPONENT_HH
