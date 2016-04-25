#include <cassert>
#include <cstdlib>

// force HDF5 1.8.x installations to use the new API
#define H5Dcreate_vers 2

#include <hdf5.h>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "vect.hh"

#include "tensor_component.hh"

namespace CarpetIOF5 {

namespace F5 {

tensor_component_t::tensor_component_t(data_region_t &data_region,
                                       int const variable)
    : m_data_region(data_region), m_variable(variable) {
  DECLARE_CCTK_PARAMETERS;

  assert(variable >= 0 and variable < CCTK_NumVars());

  char const *const name = CCTK_VarName(variable);
  assert(name != 0);
  m_name = string(name);

  int const vartype = CCTK_VarTypeI(variable);
  assert(vartype >= 0);
  hid_t const hdf5_datatype = hdf5_datatype_from_cactus_datatype(vartype);
  assert(hdf5_datatype >= 0);

  bbox<int, dim> const &region = m_data_region.get_region();

  vect<hsize_t, dim> const dims = (region.shape() / region.stride()).reverse();
  m_dataspace = H5Screate_simple(dim, &dims[0], 0);
  assert(m_dataspace >= 0);

  m_properties = H5Pcreate(H5P_DATASET_CREATE);
  assert(m_properties >= 0);
  vect<int, dim> const user_chunk_size(chunk_size_x, chunk_size_y,
                                       chunk_size_z);
  bool const need_chunks =
      any(user_chunk_size > 0) or compression_level > 0 or write_checksum;
  if (need_chunks) {
    vect<hsize_t, dim> const chunk_size =
        either(user_chunk_size > 0, vect<hsize_t, dim>(user_chunk_size), dims);
    herr_t const herr = H5Pset_chunk(m_properties, dim, &chunk_size[0]);
    assert(not herr);
  }
  if (compression_level > 0) {
    herr_t herr;
    herr = H5Pset_shuffle(m_properties);
    assert(not herr);
    const herr = H5Pset_deflate(m_properties, compression_level);
    assert(not herr);
  }
  if (write_checksum) {
    herr_t const herr = H5Pset_fletcher32(m_properties);
    assert(not herr);
  }

  m_dataset = H5Dcreate(m_data_region.get_hdf5_data_region(), m_name.c_str(),
                        hdf5_datatype, m_dataspace, H5P_DEFAULT, m_properties,
                        H5P_DEFAULT);
  assert(m_dataset >= 0);

  write_or_check_attribute(m_dataset, "iorigin",
                           region.lower() / region.stride());

  assert(invariant());
}

tensor_component_t::~tensor_component_t() {
  herr_t herr;

  herr = H5Dclose(m_dataset);
  assert(not herr);

  herr = H5Sclose(m_dataspace);
  assert(not herr);

  herr = H5Pclose(m_properties);
  assert(not herr);
}

data_region_t &tensor_component_t::get_data_region() const {
  return m_data_region;
}

hid_t tensor_component_t::get_variable() const { return m_variable; }

// TODO: This assumes that the shape of data is the same as the
// shape of the region; this may not be so if not all of the data
// are written out
void tensor_component_t::write(void const *const data,
                               int const cactus_datatype) const {
  hid_t const memory_hdf5_datatype =
      hdf5_datatype_from_cactus_datatype(cactus_datatype);
  assert(memory_hdf5_datatype >= 0);

  bbox<int, dim> const &region = m_data_region.get_region();

  vect<hsize_t, dim> const dims = (region.shape() / region.stride()).reverse();
  hid_t const memory_dataspace = H5Screate_simple(dim, &dims[0], &dims[0]);
  assert(memory_dataspace >= 0);

  hid_t const transfer_properties = H5Pcreate(H5P_DATASET_XFER);
  assert(transfer_properties >= 0);

  herr_t herr;
  herr = H5Dwrite(m_dataset, memory_hdf5_datatype, memory_dataspace,
                  m_dataspace, transfer_properties, data);
  assert(not herr);

  herr = H5Pclose(transfer_properties);
  assert(not herr);

  herr = H5Sclose(memory_dataspace);
  assert(not herr);
}

bool tensor_component_t::invariant() const {
  return (m_variable >= 0 and m_variable < CCTK_NumVars() and
          m_properties >= 0 and m_dataset >= 0 and m_dataspace >= 0);
}

} // namespace F5

} // namespace CarpetIOF5
