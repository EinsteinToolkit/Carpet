#include <cassert>
#include <complex>
#include <cstring>
#include <sstream>
#include <vector>

// force HDF5 1.8.x installations to use the new API
#define H5Acreate_vers 2
#define H5Gcreate_vers 2
#define H5Gopen_vers 2
#define H5Tarray_create_vers 2

#include <hdf5.h>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "bbox.hh"
#include "defs.hh"
#include "vect.hh"

#include "utils.hh"

namespace CarpetIOF5 {

namespace F5 {

using std::complex;
using std::vector;

hid_t hdf5_datatype_from_dummy(signed char const &dummy) {
  return H5T_NATIVE_SCHAR;
}

hid_t hdf5_datatype_from_dummy(short const &dummy) { return H5T_NATIVE_SHORT; }

hid_t hdf5_datatype_from_dummy(int const &dummy) { return H5T_NATIVE_INT; }

hid_t hdf5_datatype_from_dummy(long const &dummy) { return H5T_NATIVE_LONG; }

hid_t hdf5_datatype_from_dummy(long long const &dummy) {
  return H5T_NATIVE_LLONG;
}

hid_t hdf5_datatype_from_dummy(float const &dummy) { return H5T_NATIVE_FLOAT; }

hid_t hdf5_datatype_from_dummy(double const &dummy) {
  return H5T_NATIVE_DOUBLE;
}

hid_t hdf5_datatype_from_dummy(long double const &dummy) {
  return H5T_NATIVE_LDOUBLE;
}

#ifdef HAVE_CCTK_COMPLEX8
hid_t hdf5_datatype_from_dummy(CCTK_COMPLEX8 const &dummy) {
  CCTK_REAL4 real;
  return hdf5_complex_datatype_from_dummy(dummy, real);
}
#endif

#ifdef HAVE_CCTK_COMPLEX16
hid_t hdf5_datatype_from_dummy(CCTK_COMPLEX16 const &dummy) {
  CCTK_REAL8 real;
  return hdf5_complex_datatype_from_dummy(dummy, real);
}
#endif

#ifdef HAVE_CCTK_COMPLEX32
hid_t hdf5_datatype_from_dummy(CCTK_COMPLEX32 const &dummy) {
  CCTK_REAL16 real;
  return hdf5_complex_datatype_from_dummy(dummy, real);
}
#endif

template <typename T, typename R>
hid_t hdf5_complex_datatype_from_dummy(T const &dummy, R const &real) {
  static bool initialised = false;
  static hid_t hdf_complex_datatype;

  if (not initialised) {
    initialised = true;

    hsize_t const cdim[1] = {2};

    hdf_complex_datatype =
        H5Tarray_create(hdf5_datatype_from_dummy(real), 1, cdim);
    assert(hdf_complex_datatype >= 0);
  }

  return hdf_complex_datatype;
}

hid_t hdf5_datatype_from_cactus_datatype(int const cactus_datatype) {
  switch (cactus_datatype) {
  case CCTK_VARIABLE_VOID:
    return H5I_INVALID_HID;
#define CASE(TYPEID, TYPE)                                                     \
  case TYPEID: {                                                               \
    TYPE dummy;                                                                \
    return hdf5_datatype_from_dummy(dummy);                                    \
  } break
    CASE(CCTK_VARIABLE_BYTE, CCTK_BYTE);
    CASE(CCTK_VARIABLE_INT, CCTK_INT);
#ifdef HAVE_CCTK_INT1
    CASE(CCTK_VARIABLE_INT1, CCTK_INT1);
#endif
#ifdef HAVE_CCTK_INT2
    CASE(CCTK_VARIABLE_INT2, CCTK_INT2);
#endif
#ifdef HAVE_CCTK_INT4
    CASE(CCTK_VARIABLE_INT4, CCTK_INT4);
#endif
#ifdef HAVE_CCTK_INT8
    CASE(CCTK_VARIABLE_INT8, CCTK_INT8);
#endif
#ifdef HAVE_CCTK_INT16
    CASE(CCTK_VARIABLE_INT16, CCTK_INT16);
#endif
    CASE(CCTK_VARIABLE_REAL, CCTK_REAL);
#ifdef HAVE_CCTK_REAL4
    CASE(CCTK_VARIABLE_REAL4, CCTK_REAL4);
#endif
#ifdef HAVE_CCTK_REAL8
    CASE(CCTK_VARIABLE_REAL8, CCTK_REAL8);
#endif
#ifdef HAVE_CCTK_REAL16
    CASE(CCTK_VARIABLE_REAL16, CCTK_REAL16);
#endif
    CASE(CCTK_VARIABLE_COMPLEX, CCTK_COMPLEX);
#ifdef HAVE_CCTK_COMPLEX8
    CASE(CCTK_VARIABLE_COMPLEX8, CCTK_COMPLEX8);
#endif
#ifdef HAVE_CCTK_COMPLEX16
    CASE(CCTK_VARIABLE_COMPLEX16, CCTK_COMPLEX16);
#endif
#ifdef HAVE_CCTK_COMPLEX32
    CASE(CCTK_VARIABLE_COMPLEX32, CCTK_COMPLEX32);
#endif
#undef CASE
  case CCTK_VARIABLE_CHAR:
    return H5I_INVALID_HID;
  case CCTK_VARIABLE_STRING:
    return H5I_INVALID_HID;
  case CCTK_VARIABLE_POINTER:
    return H5I_INVALID_HID;
  case CCTK_VARIABLE_FPOINTER:
    return H5I_INVALID_HID;
  }
  return H5I_INVALID_HID;
}

template <typename T, int D> string name_from_ivect(vect<T, D> const &ivect) {
  ostringstream buf;
  for (int d = 0; d < dim; ++d) {
    if (d > 0) {
      buf << ":";
    }
    buf << ivect[d];
  }
  return buf.str();
}

template string name_from_ivect(vect<int, dim> const &ivect);

template string name_from_ivect(vect<CCTK_REAL, dim> const &ivect);

template <typename T, int D> string name_from_ibbox(bbox<T, D> const &ibbox) {
  ostringstream buf;
  buf << name_from_ivect(ibbox.lower()) << "-" << name_from_ivect(ibbox.upper())
      << "-" << name_from_ivect(ibbox.stride());
  return buf.str();
}

template string name_from_ibbox(bbox<int, dim> const &ivect);

hid_t open_or_create_group(hid_t const where, char const *const name) {
  DECLARE_CCTK_PARAMETERS;

  assert(where >= 0);
  assert(name != 0);

  // bool group_exists;
  // H5E_BEGIN_TRY {
  //   H5G_stat_t statbuf;
  //   herr_t const herr = H5Gget_objinfo (where, name, true, & statbuf);
  //   group_exists = herr == 0;
  // } H5E_END_TRY;
  bool group_exists;
  H5E_BEGIN_TRY {
    H5G_info_t groupinfo;
    herr_t const herr =
        H5Gget_info_by_name(where, name, &groupinfo, H5P_DEFAULT);
    group_exists = herr == 0;
  }
  H5E_END_TRY;
  if (veryverbose) {
    CCTK_VInfo(CCTK_THORNSTRING, "Group \"%s\" %s", name,
               group_exists ? "exists" : "does not exist");
  }

  hid_t group;
  if (group_exists) {
    if (veryverbose) {
      CCTK_VInfo(CCTK_THORNSTRING, "H5Gopen (name=\"%s\")", name);
    }
    group = H5Gopen(where, name, H5P_DEFAULT);
  } else {
    if (veryverbose) {
      CCTK_VInfo(CCTK_THORNSTRING, "H5Gcreate (name=\"%s\")", name);
    }
    group = H5Gcreate(where, name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  }

  return group;
}

template <typename T>
void write_or_check_attribute(hid_t const where, char const *const name,
                              T const *const values, int const num_values) {
  assert(where >= 0);
  assert(name != 0);
  assert(num_values >= 0);
  assert(num_values == 0 or values != 0);

  T dummy;
  hid_t const datatype = hdf5_datatype_from_dummy(dummy);

  hid_t attribute;
  H5E_BEGIN_TRY { attribute = H5Aopen_name(where, name); }
  H5E_END_TRY;

  if (attribute < 0) {
    // The attribute does not yet exist; create it
    hsize_t const adim = num_values;
    hid_t const dataspace = H5Screate_simple(1, &adim, &adim);
    assert(dataspace >= 0);
    attribute =
        H5Acreate(where, name, datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    assert(attribute >= 0);
    herr_t herr;
    herr = H5Awrite(attribute, datatype, values);
    assert(not herr);
    herr = H5Aclose(attribute);
    assert(not herr);
    herr = H5Sclose(dataspace);
    assert(not herr);
  } else {
    // The attribute already exists; read and check it
    hid_t const dataspace = H5Aget_space(attribute);
    htri_t const is_simple = H5Sis_simple(dataspace);
    assert(is_simple >= 0);
    assert(is_simple > 0);
    int const ndims = H5Sget_simple_extent_ndims(dataspace);
    assert(ndims == 1);
    hsize_t adim;
    herr_t herr;
    herr = H5Sget_simple_extent_dims(dataspace, &adim, 0);
    assert(adim == hsize_t(num_values));
    vector<T> buf(adim);
    herr = H5Aread(attribute, datatype, &buf.front());
    assert(not herr);
    herr = H5Sclose(dataspace);
    assert(not herr);
    herr = H5Aclose(attribute);
    assert(not herr);
    for (int n = 0; n < num_values; ++n) {
      assert(values[n] == buf[n]);
    }
  }
}

#if SIZEOF_INT != CCTK_INTEGER_PRECISION
template void write_or_check_attribute(hid_t const where,
                                       char const *const name,
                                       int const *const values,
                                       int const num_values);
#endif
template void write_or_check_attribute(hid_t const where,
                                       char const *const name,
                                       CCTK_INT const *const values,
                                       int const num_values);
template void write_or_check_attribute(hid_t const where,
                                       char const *const name,
                                       CCTK_REAL const *const values,
                                       int const num_values);

template <typename T>
void write_or_check_attribute(hid_t const where, char const *const name,
                              T const &value) {
  assert(where >= 0);
  assert(name != 0);

  write_or_check_attribute(where, name, &value, 1);
}

#if SIZEOF_INT != CCTK_INTEGER_PRECISION
template void write_or_check_attribute(hid_t const where,
                                       char const *const name,
                                       int const &value);
#endif
template void write_or_check_attribute(hid_t const where,
                                       char const *const name,
                                       CCTK_INT const &value);
template void write_or_check_attribute(hid_t const where,
                                       char const *const name,
                                       CCTK_REAL const &value);

template <typename T, int D>
void write_or_check_attribute(hid_t where, char const *name,
                              vect<T, D> const &value) {
  assert(where >= 0);
  assert(name != 0);

  write_or_check_attribute(where, name, &value[0], D);
}

#if SIZEOF_INT != CCTK_INTEGER_PRECISION
template void write_or_check_attribute(hid_t where, char const *name,
                                       vect<int, dim> const &value);
#endif
template void write_or_check_attribute(hid_t where, char const *name,
                                       vect<CCTK_INT, dim> const &value);
template void write_or_check_attribute(hid_t where, char const *name,
                                       vect<CCTK_REAL, dim> const &value);

template <>
void write_or_check_attribute(hid_t const where, char const *const name,
                              char const *const &value) {
  assert(where >= 0);
  assert(name != 0);
  assert(value != 0);

  hid_t attribute;
  H5E_BEGIN_TRY { attribute = H5Aopen_name(where, name); }
  H5E_END_TRY;

  if (attribute < 0) {
    // The attribute does not yet exist; create it
    hid_t const datatype = H5Tcopy(H5T_C_S1);
    assert(datatype >= 0);
    herr_t herr;
    int const length = strlen(value) + 1;
    herr = H5Tset_size(datatype, length);
    assert(not herr);
    hsize_t const dim = 1;
    hid_t const dataspace = H5Screate_simple(1, &dim, &dim);
    assert(dataspace >= 0);
    attribute =
        H5Acreate(where, name, datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    assert(attribute >= 0);
    herr = H5Awrite(attribute, datatype, value);
    assert(not herr);
    herr = H5Aclose(attribute);
    assert(not herr);
    herr = H5Sclose(dataspace);
    assert(not herr);
    herr = H5Tclose(datatype);
    assert(not herr);
  } else {
    // The attribute already exists; read and check it
    hid_t datatype = H5Aget_type(attribute);
    assert(datatype >= 0);
    hid_t typeclass = H5Tget_class(datatype);
    assert(typeclass == H5T_STRING);
    int const length = H5Tget_size(datatype);
    assert(length >= 0);
    hid_t const dataspace = H5Aget_space(attribute);
    htri_t const is_simple = H5Sis_simple(dataspace);
    assert(is_simple >= 0);
    assert(is_simple > 0);
    int const ndims = H5Sget_simple_extent_ndims(dataspace);
    assert(ndims == 1);
    hsize_t dim;
    herr_t herr;
    herr = H5Sget_simple_extent_dims(dataspace, &dim, 0);
    assert(dim == 1);
    vector<char> buf(length);
    herr = H5Aread(attribute, datatype, &buf.front());
    assert(not herr);
    herr = H5Sclose(dataspace);
    assert(not herr);
    herr = H5Aclose(attribute);
    assert(not herr);
    herr = H5Tclose(datatype);
    assert(not herr);
    assert(strcmp(&buf.front(), value) == 0);
  }
}

} // namespace F5

} // namespace CarpetIOF5
